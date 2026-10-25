#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "gtl/btree.hpp"
#define HT_SHARD_BC 8

#include <blockingconcurrentqueue.h>

#include <atomic>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <limits>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "asm.h"
#include "bloomfilter.h"
#include "config.h"
#include "inline_vec.h"
#include "lcmer.h"

template <class map, class F>
static void for_each_map(const map &mp, const F &fCallback, int threads) {
    std::vector<std::thread> workers;
    for (int tid = 0; tid < threads; ++tid) {
        workers.emplace_back(
            [&](const map &mp, int tid) {
                for (size_t i = tid; i < mp.subcnt(); i += threads) {
                    mp.with_submap(i, [&](const map::EmbeddedSet &set) {
                        for (const auto &p : set) {
                            fCallback(tid, p);
                        }
                    });
                }
            },
            std::cref(mp), tid);
    }
    for (auto &worker : workers) {
        worker.join();
    }
}

// Calls fn(i) for every i in [0, n) on up to `threads` threads (the calling
// thread is one of them), handing out indices one at a time.
template <class F>
static void parallel_for(size_t n, int threads, F &&fn) {
    std::atomic<size_t> next{0};
    auto work = [&] {
        for (size_t i = next++; i < n; i = next++) fn(i);
    };
    std::vector<std::thread> workers;
    for (size_t t = 1; t < std::min<size_t>(threads, n); ++t)
        workers.emplace_back(work);
    work();
    for (auto &worker : workers) worker.join();
}

// A queue holding at most `capacity` items, with the dequeue interface of
// moodycamel::BlockingConcurrentQueue. enqueue() blocks on a semaphore of free
// slots while the queue is full (no spinning or polling). Dequeues free slots
// in groups of up to kSlotGroup: while the queue stays full, every freed slot
// would otherwise wake a blocked producer through the kernel. A group never
// exceeds the capacity, so a full queue always holds enough items to free one.
// Blocks for `capacity` items from `producers` producer threads are allocated
// up front, but enqueue() can still allocate one: a block returns to the pool
// only once all of its items are dequeued, so partly dequeued blocks can use
// up the pool. Blocks allocated this way are freed when emptied.
template <class T>
class bounded_queue {
   public:
    // Enough to make the wake-up negligible next to a producer's work for the
    // group, small next to the capacities in use.
    static constexpr size_t kSlotGroup = 16;
    bounded_queue(size_t capacity, size_t producers)
        : queue_{capacity, 0, producers},
          free_slots_{static_cast<ssize_t>(capacity), /*maxSpins=*/0},
          group_{std::min(capacity, kSlotGroup)} {}
    void enqueue(T item) {
        free_slots_.wait();
        queue_.enqueue(std::move(item));
    }
    void wait_dequeue(T &item) {
        queue_.wait_dequeue(item);
        free_slot();
    }
    template <class Rep, class Period>
    bool wait_dequeue_timed(T &item,
                            std::chrono::duration<Rep, Period> timeout) {
        if (!queue_.wait_dequeue_timed(item, timeout)) return false;
        free_slot();
        return true;
    }
    bool try_dequeue(T &item) {
        if (!queue_.try_dequeue(item)) return false;
        free_slot();
        return true;
    }

   private:
    void free_slot() {
        if (dequeued_.fetch_add(1) % group_ == group_ - 1) {
            free_slots_.signal(static_cast<ssize_t>(group_));
        }
    }
    moodycamel::BlockingConcurrentQueue<T> queue_;
    moodycamel::LightweightSemaphore free_slots_;
    const size_t group_;
    std::atomic<size_t> dequeued_{0};
};

template <std::integral int_type>
consteval int_type nextpowerof2(int_type v) {
    v--;
    for (size_t i = 1; i < sizeof(int_type) * 8; i <<= 1) {
        v |= v >> i;
    }
    v++;
    return v;
}

using u64 = uint64_t;
using i64 = int64_t;
using u32 = uint32_t;
using DefaultBF = bloom_interface<bloom_ext::BlockedAtomicBloomFilter<>>;
// A bubble removed by GraphCleaning::pop_bubbles, for a planned diploid
// output of the removed allele as an alternate sequence (collected only under
// dbg_index_base::kKeepPoppedAlleles). The fragments of the source's last and
// the sink's first lcmer keep their origin ids through every later merge, so
// they locate the bubble in the final unitigs.
struct popped_bubble {
    u64 source_origin;
    u64 sink_origin;
    bool error;                 // removed nodes read below half a copy
    float removed_depth;        // highest avg_core_count among removed nodes
    float kept_depth;           // lowest avg_core_count on the kept path
    SequenceContainerType alt;  // best removed path, source to sink
};
// The word_arenas of one graph's fragments, one per worker id
// (0..MAXIMUM_WORK_THREADS).
using fragment_arenas = std::array<word_arena, MAXIMUM_WORK_THREADS + 1>;

template <class BFT = DefaultBF>
struct dbg_index_base {
    using BloomFilterType = BFT;
    using mutex_type = std::mutex;
    using fragment_list = inline_vec<SequenceContainerType, 1>;
    using prev_list = inline_vec<u64, 1>;
    struct bidbg_val {
        struct counts_t {
            unsigned core_count;
        };

        bidbg_val(unsigned core_count)
            : redirected_flag(false),
              seeded(false),
              forward_edge_count(0),
              counts{core_count} {}
        bidbg_val(unsigned core_count, u64 prev)
            : redirected_flag(false),
              seeded(false),
              forward_edge_count(0),
              counts{core_count} {
            this->prev.push_back(prev);
        }

        bool redirected_flag;
        // Holds a window of an earlier pass's segment (--multipass): admitted
        // without the Bloom filter, adding no count, and spared by the
        // depth-only cleaning rules. Sits in the padding before
        // forward_edge_count.
        bool seeded;
        // forward_edge_count stays outside the redirect_to union: a
        // neighbouring chain may read it while this node is being redirected.
        // fec is written only on representatives; core_count stays unioned
        // since only the node's own chain touches it.
        unsigned forward_edge_count;
        union {
            counts_t counts;
            u64 redirect_to;
        };

        // The node's sequence: one fragment per absorbed lcmer node, in path
        // order. A fresh node has one, held inline; its words come from the
        // inserting thread's fragment_arena.
        fragment_list seq;

        bool redirected() const { return redirected_flag; }
        // An absorbed node keeps only its redirect target: seq and prev are
        // released, so it holds no heap and lists no predecessors.
        void mark_redirected(u64 to) {
            redirected_flag = true;
            redirect_to = to;
            seq = fragment_list{};
            prev = prev_list{};
        }

        // Predecessor node ids; a node usually has one, held inline.
        prev_list prev;
        void emplace_prev(u64 val) {
            if (val == 0) return;
            auto it = std::find(prev.begin(), prev.end(), val);
            if (it == prev.end()) {
                prev.push_back(val);
            }
        }

        constexpr double avg_core_count() const {
            // Precondition: live node. seq.empty() means the node was absorbed
            // during simplification (mark_redirected); skip via redirected()
            // first.
            assert(!seq.empty());
            return seq.empty() ? 0.0 : (double)counts.core_count / seq.size();
        }
    };
    // One fragment and one predecessor inline keep a node at one 64 B value.
    static_assert(sizeof(bidbg_val) == 64);
    using dbg_map_t = gtl::parallel_flat_hash_map<
        u64, bidbg_val, std::hash<u64>, gtl::priv::hash_default_eq<u64>,
        std::allocator<std::pair<const u64, bidbg_val>>, HT_SHARD_BC,
        mutex_type>;

    struct lci_type {
        u64 id;
        static consteval lci_type empty() { return {0}; }
    };

    dbg_map_t dbg;
    gtl::btree_set<u64> active_set;
    // Holds popped alleles for planned diploid alternate-allele output.
    static constexpr bool kKeepPoppedAlleles = false;
    gtl::vector<popped_bubble> popped_bubbles;
    // Depth above which a node holds more than one locus (set by cleaning;
    // 0 when no depth model was fitted).
    double locus_limit = 0.0;

    // First-sighting filter for read lcmers. Only read insertion
    // (construct_dbg) probes it, and construct_dbg releases it when done.
    std::optional<BloomFilterType> bf;
    // --retro-admission: the first sighting of a window, kept when it is on
    // a read's own strand, with its neighbouring window keys in that read.
    // `prev_of` holds every kept key with the preceding key (0 at the read's
    // start). The following key is 0 at the read's end, and also when that
    // window was kept from the same read, whose sighting then holds the
    // adjacency. Kept windows mostly come in runs from one read, so it is 0
    // for most sightings; `next_of` holds only the nonzero ones. Only read
    // insertion and resolve_retro_admissions use them; construct_dbg frees
    // them.
    struct retro_admissions {
        using table = gtl::parallel_flat_hash_map<
            u64, u64, std::hash<u64>, gtl::priv::hash_default_eq<u64>,
            std::allocator<std::pair<const u64, u64>>, HT_SHARD_BC, mutex_type>;
        table prev_of;
        table next_of;
        // A key's Bloom filter test and its record are made under the key's
        // stripe, so a read that finds the key in the filter also finds the
        // record, and two reads never both take the first sighting. 32
        // stripes per possible worker keep two workers on one stripe rare.
        static constexpr size_t stripes_per_worker = 32;
        static constexpr size_t stripe_count =
            nextpowerof2((u64)MAXIMUM_WORK_THREADS) * stripes_per_worker;
        std::array<mutex_type, stripe_count> stripes;
        mutex_type &stripe(u64 key) { return stripes[key % stripe_count]; }
    };
    std::unique_ptr<retro_admissions> retro;
    // --retro-admission, a window of a read's own strand: the key if the
    // Bloom filter has seen it, else 0 and the sighting is kept.
    u64 sight_own(u64 key, u64 prev, u64 next) {
        // A repeat sighting needs no record and no stripe: bits are only
        // ever set, so a filter that holds the key now also holds it under
        // the stripe.
        if (bf->contains(key)) return key;
        std::lock_guard lock(retro->stripe(key));
        if (bf->insert_if_absent(key)) {
            retro->prev_of.try_emplace(key, prev);
            if (next != 0) retro->next_of.try_emplace(key, next);
            return 0;
        }
        return key;
    }
    // A window of a read's reverse-complement parse: inserts it into the
    // Bloom filter. `admitted`: the filter had already seen the window.
    // `needs_node` (--retro-admission only): the window is now admitted with
    // a kept own-strand sighting, so it needs a node this read's bases can
    // give.
    struct reverse_sighting {
        bool admitted;
        bool needs_node;
    };
    reverse_sighting sight_reverse(u64 key) {
        if (!retro) {
            return {!bf->insert_if_absent(key), false};
        }
        std::lock_guard lock(retro->stripe(key));
        const bool admitted = !bf->insert_if_absent(key);
        return {admitted, admitted && retro->prev_of.contains(key)};
    }
    // After reading, adds each kept sighting of an admitted window as its
    // read would have added it had the window been admitted then
    // (insert_lcmer_helper); a window is admitted iff it has a node. Each
    // adjacency with an admitted neighbour gives the upstream window a count
    // and the downstream one a predecessor. Every adjacency is held by one
    // sighting (see retro_admissions) and nothing else changes the graph
    // meanwhile, so the result does not depend on thread timing. Returns the
    // number of admitted kept sightings.
    u64 resolve_retro_admissions(int threads) {
        const auto add_count = [](auto &val) {
            val.second.counts.core_count += 1;
        };
        std::atomic<u64> admitted{0};
        parallel_for(retro->prev_of.subcnt(), threads, [&](size_t submap) {
            retro->prev_of.with_submap(submap, [&](const auto &set) {
                for (const auto &[key, prev] : set) {
                    if (!dbg.contains(key)) continue;
                    admitted.fetch_add(1, std::memory_order_relaxed);
                    if (dbg.modify_if(prev, add_count)) {
                        dbg.modify_if(key, [&](auto &val) {
                            val.second.emplace_prev(prev);
                        });
                    }
                    // Reading is over and nothing writes next_of any more,
                    // so the lookup needs no lock.
                    u64 next = 0;
                    retro->next_of.if_contains_unsafe(
                        key, [&](const auto &v) { next = v.second; });
                    if (dbg.modify_if(next, [&](auto &val) {
                            val.second.emplace_prev(key);
                        })) {
                        dbg.modify_if(key, add_count);
                    }
                }
            });
        });
        return admitted.load();
    }
    // Erases the nodes made only to hold a place for a kept sighting's
    // adjacency (count 0, no predecessor). The own-strand path and
    // resolve_retro_admissions give a count to every node they list as a
    // predecessor, so none of these is referenced by them. --link-reverse may
    // list one as a predecessor without a count; that edge goes with the
    // node, as cleaning would cut such a tip anyway. Seeding comes after
    // construction, so no node is seeded yet. Returns the number erased.
    u64 erase_unlinked_placeholders(int threads) {
        std::atomic<u64> erased{0};
        parallel_for(dbg.subcnt(), threads, [&](size_t submap) {
            dbg.with_submap_m(submap, [&](auto &set) {
                for (auto it = set.begin(); it != set.end();) {
                    const bidbg_val &v = it->second;
                    if (v.counts.core_count == 0 && v.prev.empty()) {
                        it = set.erase(it);
                        erased.fetch_add(1, std::memory_order_relaxed);
                    } else {
                        ++it;
                    }
                }
            });
        });
        return erased.load();
    }

    dbg_index_base(double false_positive_rate, uint64_t projected_count,
                   uint64_t seed)
        : bf{std::in_place, false_positive_rate, projected_count, seed} {}

    // The words of node fragments, from the arena of the inserting worker
    // (`tid` 0..MAXIMUM_WORK_THREADS, one thread per id at a time). Nodes
    // borrow them; they are freed with the index or by clear_parallel,
    // unless take_fragment_arenas hands them to a new owner first.
    word_arena &fragment_arena(int tid) { return arenas_[tid]; }
    // Moves the arenas out (unitig_graph takes them with the fragments).
    fragment_arenas take_fragment_arenas() { return std::move(arenas_); }

   private:
    fragment_arenas arenas_;

   public:
    static constexpr u64 lcmer_id(const LCmer &l) { return l.value(); }

    // The node an lcmer id ends at after following redirects, or empty when
    // the chain ends at an erased node. Only path simplification redirects:
    // it merges a chain into the chain's most upstream node, which is live
    // and not merged in the same sweep, and a redirected node is never merged
    // or revived again. Each redirect thus points to a node redirected in a
    // later sweep or never, so every chain ends, after at most one hop per
    // sweep.
    lci_type resolve(u64 id) const {
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;
        auto &dbg_nc = const_cast<dbg_map_t &>(dbg);
        auto it = dbg_nc.find(id, dbg_nc.hash(id), no_lock);
        while (it != dbg_nc.end() && it->second.redirected()) {
            id = it->second.redirect_to;
            it = dbg_nc.find(id, dbg_nc.hash(id), no_lock);
        }
        if (it == dbg_nc.end()) return lci_type::empty();
        return {id};
    }
    // Sets every node's forward_edge_count to the number of nodes listing it
    // in prev. Submaps are split across threads; nothing is inserted or
    // erased meanwhile, so lookups take no lock, and counts are bumped
    // atomically (a node can be the successor of nodes in any submap).
    void count_forward_edges(int threads) {
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;
        parallel_for(dbg.subcnt(), threads, [&](size_t submap) {
            dbg.with_submap(submap, [&](const auto &set) {
                for (const auto &node : set) {
                    for (u64 edge : node.second.prev) {
                        auto it = dbg.find(edge, dbg.hash(edge), no_lock);
                        if (it == dbg.end()) continue;
                        std::atomic_ref<unsigned>(it->second.forward_edge_count)
                            .fetch_add(1, std::memory_order_relaxed);
                    }
                }
            });
        });
    }
    // Destroys every node, one submap per task on up to `threads` threads:
    // a large graph's per-node allocations take long to free on one thread.
    // The index stays usable (empty).
    void clear_parallel(int threads) {
        parallel_for(dbg.subcnt(), threads,
                     [&](size_t submap) { dbg.clear(submap); });
        for (auto &a : arenas_) a.release();
        gtl::btree_set<u64>().swap(active_set);
        gtl::vector<popped_bubble>().swap(popped_bubbles);
    }
    u64 emplace_nobloomfilter(const LCmer &l, u64 &idx, int tid) {
        (void)idx;
        (void)tid;
        return lcmer_id(l);
    }
    u64 emplace(const LCmer &l, u64 &idx, int tid) {
        (void)idx;
        (void)tid;
        const auto new_hash = std::hash<LCmer>{}(l);
        if (bf->insert_if_absent(new_hash)) {
            return 0;
        }
        return lcmer_id(l);
    }
};

using dbg_index = dbg_index_base<DefaultBF>;

// The cleaned graph in compact form. After cleaning nearly every node-table
// entry is a redirect stub, kept only so that an lcmer id resolves to its
// final unitig (dbg_index::resolve). unitig_graph holds the live unitigs in a
// dense array and the lcmer ids in a sorted index that answers every id as
// resolve() did, so the node table is freed once cleaning ends.
struct unitig_graph {
    using lci_type = dbg_index::lci_type;
    struct unitig {
        u64 id;
        unsigned core_count;
        dbg_index::fragment_list seq;  // words borrowed from arenas_
        dbg_index::prev_list prev;     // ids of live unitigs

        double avg_core_count() const {
            if (seq.empty()) return 0.0;
            return (double)core_count / seq.size();
        }
    };

    // Live unitigs in node-table iteration order; those of node-table submap
    // s are [submap_begin[s], submap_begin[s + 1]).
    gtl::vector<unitig> unitigs;
    gtl::vector<size_t> submap_begin;
    // Cleaning results read by later stages (see dbg_index_base); active_set
    // may still name erased nodes, so look ids up with find().
    gtl::btree_set<u64> active_set;
    gtl::vector<popped_bubble> popped_bubbles;
    double locus_limit = 0.0;
    size_t node_table_size = 0;  // entries of the node table, stubs included

    // Takes the live unitigs and the cleaning results out of `cleaned` and
    // clears its node table; the caller then frees `cleaned`.
    unitig_graph(dbg_index &cleaned, int threads);

    // The live unitig `id`, or nullptr.
    const unitig *find(u64 id) const {
        auto it = index_of_.find(id);
        if (it == index_of_.end()) return nullptr;
        return &unitigs[it->second];
    }
    // What dbg_index::resolve(key) returned on the cleaned node table.
    lci_type resolve(u64 key) const {
        auto it = std::lower_bound(keys_.begin(), keys_.end(), key);
        if (it == keys_.end() || *it != key) return lci_type::empty();
        return {unitigs[targets_[it - keys_.begin()]].id};
    }
    lci_type get(const LCmer &l) const {
        return resolve(dbg_index::lcmer_id(l));
    }

    // fn(tid, unitig) for every unitig on `threads` threads. Thread tid takes
    // the unitigs of submaps tid, tid + threads, ...: the split for_each_map
    // made over the node table, so per-thread results combine as before.
    template <class F>
    void for_each_unitig(int threads, const F &fn) const {
        std::vector<std::thread> workers;
        for (int tid = 0; tid < threads; ++tid) {
            workers.emplace_back([&, tid] {
                for (size_t s = tid; s + 1 < submap_begin.size(); s += threads)
                    for (size_t u = submap_begin[s]; u < submap_begin[s + 1];
                         ++u)
                        fn(tid, unitigs[u]);
            });
        }
        for (auto &worker : workers) worker.join();
    }

    // Frees everything; the fragments, many small allocations, are freed one
    // unitig per task on up to `threads` threads.
    void clear_parallel(int threads) {
        parallel_for(unitigs.size(), threads, [&](size_t u) {
            unitigs[u].seq = dbg_index::fragment_list{};
        });
        *this = unitig_graph{};
    }

   private:
    unitig_graph() = default;
    // Steps of the constructor: the unitig slots, ids and (unsorted) lcmer
    // ids; the lcmer-id index; then the unitigs' data.
    void place_unitigs(const dbg_index::dbg_map_t &dbg, int threads);
    void index_keys(const dbg_index &cleaned, int threads);
    void take_unitig_data(dbg_index::dbg_map_t &dbg, int threads);

    gtl::flat_hash_map<u64, u32> index_of_;  // live id -> position in unitigs
    // Every lcmer id that resolves to a node, sorted, and per id that node's
    // position in unitigs.
    gtl::vector<u64> keys_;
    gtl::vector<u32> targets_;
    // The arenas holding the unitigs' fragment words, taken over from the
    // node table: they live and die with this graph.
    fragment_arenas arenas_;
};

inline unitig_graph::unitig_graph(dbg_index &cleaned, int threads)
    : active_set(std::move(cleaned.active_set)),
      popped_bubbles(std::move(cleaned.popped_bubbles)),
      locus_limit(cleaned.locus_limit),
      node_table_size(cleaned.dbg.size()) {
    place_unitigs(cleaned.dbg, threads);
    index_keys(cleaned, threads);
    take_unitig_data(cleaned.dbg, threads);
    arenas_ = cleaned.take_fragment_arenas();
    cleaned.clear_parallel(threads);
}

inline void unitig_graph::place_unitigs(const dbg_index::dbg_map_t &dbg,
                                        int threads) {
    const size_t submaps = dbg.subcnt();
    // Live and total entries per submap, then their prefix sums.
    submap_begin.assign(submaps + 1, 0);
    gtl::vector<size_t> key_begin(submaps + 1, 0);
    parallel_for(submaps, threads, [&](size_t s) {
        dbg.with_submap(s, [&](const auto &set) {
            for (const auto &node : set)
                if (!node.second.redirected()) ++submap_begin[s + 1];
            key_begin[s + 1] = set.size();
        });
    });
    for (size_t s = 0; s < submaps; ++s) {
        submap_begin[s + 1] += submap_begin[s];
        key_begin[s + 1] += key_begin[s];
    }
    unitigs.resize(submap_begin.back());
    keys_.resize(key_begin.back());
    parallel_for(submaps, threads, [&](size_t s) {
        dbg.with_submap(s, [&](const auto &set) {
            size_t u = submap_begin[s];
            size_t k = key_begin[s];
            for (const auto &node : set) {
                keys_[k++] = node.first;
                if (!node.second.redirected()) unitigs[u++].id = node.first;
            }
        });
    });
    index_of_.reserve(unitigs.size());
    for (size_t u = 0; u < unitigs.size(); ++u)
        index_of_.emplace(unitigs[u].id, static_cast<u32>(u));
}

// Resolves every lcmer id through the node table; kNoNode marks the ids that
// resolve to no node until they are dropped.
inline void unitig_graph::index_keys(const dbg_index &cleaned, int threads) {
    constexpr u32 kNoNode = std::numeric_limits<u32>::max();
    assert(unitigs.size() < kNoNode);
    std::sort(keys_.begin(), keys_.end());
    targets_.resize(keys_.size());
    constexpr size_t kKeysPerTask = size_t{1} << 16;  // amortizes the hand-out
    parallel_for((keys_.size() + kKeysPerTask - 1) / kKeysPerTask, threads,
                 [&](size_t task) {
                     const size_t end =
                         std::min(keys_.size(), (task + 1) * kKeysPerTask);
                     for (size_t i = task * kKeysPerTask; i < end; ++i) {
                         auto it = index_of_.find(cleaned.resolve(keys_[i]).id);
                         targets_[i] = kNoNode;
                         if (it != index_of_.end()) targets_[i] = it->second;
                     }
                 });
    size_t kept = 0;
    for (size_t i = 0; i < keys_.size(); ++i) {
        if (targets_[i] == kNoNode) continue;
        keys_[kept] = keys_[i];
        targets_[kept++] = targets_[i];
    }
    keys_.resize(kept);
    targets_.resize(kept);
    keys_.shrink_to_fit();
    targets_.shrink_to_fit();
}

inline void unitig_graph::take_unitig_data(dbg_index::dbg_map_t &dbg,
                                           int threads) {
    const size_t submaps = dbg.subcnt();
    parallel_for(submaps, threads, [&](size_t s) {
        dbg.with_submap_m(s, [&](auto &set) {
            size_t u = submap_begin[s];
            for (auto &node : set) {
                if (node.second.redirected()) continue;
                unitig &dst = unitigs[u++];
                dst.core_count = node.second.counts.core_count;
                dst.seq = std::move(node.second.seq);
                dst.prev = std::move(node.second.prev);
            }
        });
    });
}

struct GraphCleaning {
    using bidbg_val = dbg_index::bidbg_val;
    using lci_type = dbg_index::lci_type;
    using dbg_map_t = dbg_index::dbg_map_t;
    using redir_map = gtl::parallel_flat_hash_map_d<u64, u64, 4, std::mutex>;
    static double median_of_sorted(const gtl::vector<double> &vals) {
        if (vals.empty()) return 0.0;
        const size_t mid = vals.size() / 2;
        if (vals.size() % 2 == 1) return vals[mid];
        return (vals[mid - 1] + vals[mid]) / 2.0;
    }

    // Depth (avg_core_count) of one copy: the lowest peak of the depth
    // histogram of long nodes, weighted by cores, that holds at least 10% of
    // the long-node cores within +-25% of it. Long nodes average many
    // lcmers, so their depths are sharp. Returns 0 without long nodes.
    static double one_copy_depth(const dbg_map_t &dbg,
                                 const gtl::btree_set<u64> &active_set) {
        constexpr size_t kMinCores = 32;
        constexpr int kBins = 512;
        constexpr double kLogLow = -3.0, kLogStep = 0.025;
        gtl::vector<double> hist(kBins, 0.0), smooth(kBins, 0.0);
        double total = 0.0;
        for (u64 node : active_set) {
            auto it = dbg.find(node);
            if (it == dbg.end() || it->second.redirected()) continue;
            const size_t len = it->second.seq.size();
            if (len < kMinCores) continue;
            const double d = it->second.avg_core_count();
            if (d <= 0.0) continue;
            const int bin =
                std::clamp(static_cast<int>((std::log(d) - kLogLow) / kLogStep),
                           0, kBins - 1);
            hist[bin] += len;
            total += len;
        }
        if (total == 0.0) return 0.0;
        const auto bin_depth = [](int bin) {
            return std::exp(kLogLow + (bin + 0.5) * kLogStep);
        };
        for (int bin = 0; bin < kBins; ++bin) {
            for (int k = -3; k <= 3; ++k) {
                if (bin + k >= 0 && bin + k < kBins)
                    smooth[bin] += hist[bin + k] * (4 - std::abs(k));
            }
        }
        for (int bin = 1; bin + 1 < kBins; ++bin) {
            if (smooth[bin] < smooth[bin - 1] ||
                smooth[bin] <= smooth[bin + 1]) {
                continue;  // not a peak
            }
            const double peak = bin_depth(bin);
            double near = 0.0;
            for (int other = 0; other < kBins; ++other) {
                const double d = bin_depth(other);
                if (d >= 0.75 * peak && d <= 1.25 * peak) near += hist[other];
            }
            if (near >= 0.1 * total) return peak;
        }
        return 0.0;
    }
    static auto attempt_simplify_one_path_lockless(dbg_map_t &dbg,
                                                   redir_map &redirect_deleted,
                                                   u64 id,
                                                   auto starting_from) -> bool {
        // Find a node we cannot simplify
        // Iterate backwards and merge until we cannot simplify
        // o->o->o->o->o->o -> o->o
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;

        auto current = starting_from;
        if (current == dbg.end()) {
            return false;
        }
        unsigned first_fec = current->second.forward_edge_count;
        bool simplified = false;
        auto prev_hash = dbg.hash(current->second.prev[0]);
        auto prev_it = dbg.find(current->second.prev[0], prev_hash, no_lock);
        gtl::vector<dbg_map_t::iterator> merged;
        unsigned cc_sum = 0;
        while (prev_it != dbg.end() && current->second.prev.size() == 1 &&
               prev_it->first != current->first &&  // never step onto a
                                                    // self-loop
               prev_it->second.forward_edge_count == 1) {
            if (prev_it != dbg.end()) {
                cc_sum += current->second.counts.core_count;
                merged.push_back(current);

                simplified = true;
            } else {  // Prev is erased
            }
            current = prev_it;
            if (prev_it->second.prev.size() > 0) {
                prev_hash = dbg.hash(current->second.prev[0]);
                prev_it = dbg.find(current->second.prev[0], prev_hash, no_lock);
            } else {
                break;
            }
        }
        if (!simplified) {
            // Degenerate: nothing merged (current == seed). Skip the seed
            // write: it is value-identical but would race a neighbour chain's
            // read of forward_edge_count.
            return false;
        }
        // simplified implies the walk moved, i.e. current != seed. Assign
        // under the submap lock: operator[] returns a reference after the
        // lock is released, and a concurrent insert may rehash the submap.
        redirect_deleted.insert_or_assign(id, current->first);

        current->second.forward_edge_count = first_fec;
        current->second.counts.core_count += cc_sum;
        size_t merged_frags = current->second.seq.size();
        for (const auto &m : merged) merged_frags += m->second.seq.size();
        current->second.seq.reserve(merged_frags);
        for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
            for (auto &frag : (*it)->second.seq) {
                current->second.seq.push_back(std::move(frag));
            }
            current->second.seeded |= (*it)->second.seeded;
            (*it)->second.mark_redirected(current->first);
        }

        return true;
    };

    static auto simplify_simple_parallel_mark(dbg_index::dbg_map_t &dbg,
                                              int threads) {
        gtl::vector<gtl::vector<std::pair<u64, dbg_index::dbg_map_t::iterator>>>
            to_simplify_from(threads);
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;

        for_each_map(
            dbg,
            [&to_simplify_from, &no_lock, &dbg](
                int tid, const dbg_index::dbg_map_t::value_type &item) {
                if (item.second.forward_edge_count != 1 &&
                    item.second.prev.size() == 1) {
                    auto ih = dbg.hash(item.first);
                    auto iit = dbg.find(item.first, ih, no_lock);
                    to_simplify_from[tid].emplace_back(item.first, iit);
                } else if (item.second.prev.size() > 1) {
                    for (auto &prev_id : item.second.prev) {
                        if (prev_id == item.first) {
                            continue;  // never seed a self-loop
                        }
                        auto prev_hash = dbg.hash(prev_id);

                        auto prev_it = dbg.find(prev_id, prev_hash, no_lock);
                        // Seed only fec==1 predecessors; fec!=1 predecessors
                        // are already type-1 seeds (above). This keeps type-1
                        // and type-2 seeds disjoint so each chain is owned by
                        // one thread.
                        if (prev_it != dbg.end() &&
                            prev_it->second.prev.size() == 1 &&
                            prev_it->second.forward_edge_count == 1) {
                            to_simplify_from[tid].emplace_back(prev_id,
                                                               prev_it);
                        }
                    }
                }
            },
            threads);
        return to_simplify_from;
    }
    static void simplify_simple_parallel_sweep(
        dbg_index::dbg_map_t &dbg,
        const gtl::vector<
            gtl::vector<std::pair<u64, dbg_index::dbg_map_t::iterator>>>
            &to_simplify,
        redir_map &redirect_deleted) {
        // Chains are node-disjoint (see the seed guard in
        // simplify_simple_parallel_mark), so each thread mutates only its own
        // chain's nodes without synchronization. redirect_deleted is shared:
        // it is written only through its locking insert_or_assign.
        gtl::vector<std::thread> threads;
        for (size_t i = 0; i < to_simplify.size(); ++i) {
            threads.emplace_back(
                [&](int i) {
                    for (const auto &[node, it] : to_simplify[i]) {
                        attempt_simplify_one_path_lockless(
                            dbg, redirect_deleted, node, it);
                    }
                },
                i);
        }

        for (auto &thread : threads) {
            thread.join();
        }
    }
    // Prune absorbed nodes: a merged node is redirected() and a consumed
    // node is gone from dbg.
    static void prune_active_set(const dbg_map_t &dbg, auto &active_set) {
        for (auto it = active_set.begin(); it != active_set.end();) {
            auto nit = dbg.find(*it);
            if (nit == dbg.end() || nit->second.redirected()) {
                it = active_set.erase(it);
            } else {
                ++it;
            }
        }
    }
    // active_set from scratch: every dbg node, or only the live
    // (non-redirected) ones. Built from sorted keys in one pass instead of
    // one btree insert (or erase) per node.
    static void build_active_set(const dbg_map_t &dbg,
                                 gtl::btree_set<u64> &active_set,
                                 bool live_only) {
        gtl::vector<u64> keys;
        keys.reserve(live_only ? 0 : dbg.size());
        for (const auto &node : dbg) {
            if (!live_only || !node.second.redirected()) {
                keys.push_back(node.first);
            }
        }
        std::sort(keys.begin(), keys.end());
        active_set = gtl::btree_set<u64>(keys.begin(), keys.end());
    }
    // Rewrite prev edges through redirect_deleted and drop the self-edges
    // this creates. Only live nodes can change: a redirected node's prev is
    // empty. Every live node is in active_set (nodes only ever leave it), so
    // iterating it instead of all of dbg (mostly redirected husks) gives the
    // same result; each node edits only itself.
    static void redirect_prev_edges(dbg_map_t &dbg,
                                    const redir_map &redirect_deleted,
                                    const auto &active_set) {
        for (u64 id : active_set) {
            auto nit = dbg.find(id);
            if (nit == dbg.end()) continue;
            auto &node = *nit;
            for (size_t idx = 0; idx < node.second.prev.size(); ++idx) {
                auto &p = node.second.prev[idx];
                auto it = redirect_deleted.find(p);
                if (it != redirect_deleted.end()) {
                    p = it->second;
                }
                if (p == node.first) {
                    // Cycle collapse: the redirect rewrite made this prev entry
                    // a self-reference. Drop the spurious self-edge
                    // (emplace_prev dedups, so at most one exists).
                    if (node.second.forward_edge_count > 0) {
                        node.second.forward_edge_count--;
                    }
                    node.second.prev.erase(node.second.prev.begin() + idx);
                    --idx;
                }
            }
        }
    }
    // Unified simple-path simplification: discover node-disjoint chains in
    // parallel (mark), merge them lock-free (sweep), prune active_set and
    // redirect the remaining edges. Used for every cleaning iteration;
    // threads==1 just means a single partition.
    static void simplify_simple_paths(dbg_index::dbg_map_t &dbg,
                                      redir_map &redirect_deleted,
                                      auto &active_set, int threads) {
        auto to_simplify = simplify_simple_parallel_mark(dbg, threads);
        simplify_simple_parallel_sweep(dbg, to_simplify, redirect_deleted);
        prune_active_set(dbg, active_set);
        redirect_prev_edges(dbg, redirect_deleted, active_set);
    }

    // Returns the number of nodes removed.
    static size_t cut_tips(dbg_map_t &dbg, int min_core_depth,
                           size_t max_tip_length, auto &active_set) {
        const auto prune_prev = [&dbg](auto current) {
            auto &pvec = current->second.prev;
            pvec.erase(std::remove_if(pvec.begin(), pvec.end(),
                                      [&dbg](u64 pid) {
                                          auto it = dbg.find(pid);
                                          return it == dbg.end() ||
                                                 it->second.redirected();
                                      }),
                       pvec.end());
            std::sort(pvec.begin(), pvec.end());
        };
        size_t n_cut = 0;

        for (u64 node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end()) {
                continue;
            }
            if (current->second.redirected()) {
                continue;
            }
            // Seeded nodes take the branch of a supported node: they are
            // never removed for depth, but still drop low-depth tips.
            if (current->second.avg_core_count() < min_core_depth &&
                !current->second.seeded) {
                // Forlorn
                if (current->second.forward_edge_count == 0 &&
                    current->second.prev.size() == 0) {
                    dbg.erase(current);
                    ++n_cut;
                } else if (current->second.forward_edge_count == 0 &&
                           current->second.prev.size() == 1) {
                    auto prev_it = dbg.find(current->second.prev[0]);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        prev_it->second.forward_edge_count--;
                    }
                    dbg.erase(current);
                    ++n_cut;
                }
            } else {
                gtl::vector<u64> prev_ids(current->second.prev.begin(),
                                          current->second.prev.end());
                std::sort(prev_ids.begin(), prev_ids.end());
                for (u64 pid : prev_ids) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() < min_core_depth &&
                            !prev_it->second.seeded &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);
                            ++n_cut;
                        }
                    }
                }

                current = dbg.find(node);
                if (current == dbg.end() || current->second.redirected()) {
                    continue;
                }
                prune_prev(current);

                prev_ids.assign(current->second.prev.begin(),
                                current->second.prev.end());
                for (u64 pid : prev_ids) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() <
                                current->second.avg_core_count() &&
                            !prev_it->second.seeded &&
                            prev_it->second.seq.size() < max_tip_length &&
                            current->second.seq.size() > 2 * max_tip_length &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);
                            ++n_cut;
                        }
                    }
                }

                current = dbg.find(node);
                if (current == dbg.end() || current->second.redirected()) {
                    continue;
                }
                prune_prev(current);
            }
        }
        return n_cut;
    }

    // Live forward adjacency: successor lists derived from prev edges,
    // sorted per node for deterministic arm walking. Iterates the live
    // active_set, not the whole graph. `live(id)` is the live node `id` or
    // nullptr.
    template <class Live>
    static gtl::flat_hash_map<u64, gtl::vector<u64>> build_live_fwd(
        const gtl::btree_set<u64> &active_set, const Live &live) {
        gtl::flat_hash_map<u64, gtl::vector<u64>> fwd;
        fwd.reserve(active_set.size());
        for (u64 id : active_set) {
            const auto *node = live(id);
            if (node == nullptr) continue;
            for (u64 p : node->prev) {
                if (live(p) != nullptr) fwd[p].push_back(id);
            }
        }
        for (auto &[_, succs] : fwd) {
            std::sort(succs.begin(), succs.end());
        }
        return fwd;
    }
    static gtl::flat_hash_map<u64, gtl::vector<u64>> build_live_fwd(
        const dbg_map_t &dbg, const gtl::btree_set<u64> &active_set) {
        return build_live_fwd(active_set, [&](u64 id) -> const bidbg_val * {
            auto it = dbg.find(id);
            if (it == dbg.end() || it->second.redirected()) return nullptr;
            return &it->second;
        });
    }
    static gtl::flat_hash_map<u64, gtl::vector<u64>> build_live_fwd(
        const unitig_graph &graph) {
        return build_live_fwd(graph.active_set,
                              [&](u64 id) { return graph.find(id); });
    }

    // Bubble search from fork `sid` (Onodera et al. 2013): a node is visited
    // once all its predecessors are visited; the bubble closes when exactly
    // one reached node is left, the sink. Visited nodes other than `sid` go
    // to `inner`, each after its predecessors. Returns the sink, or 0 on a
    // dead end, a stale edge, a cycle, a self-loop, or when the inner nodes
    // exceed max_nodes nodes or max_len cores.
    static u64 find_bubble_sink(
        const dbg_map_t &dbg,
        const gtl::flat_hash_map<u64, gtl::vector<u64>> &fwd, u64 sid,
        size_t max_nodes, size_t max_len, gtl::vector<u64> &inner) {
        inner.clear();
        gtl::flat_hash_set<u64> visited;
        gtl::flat_hash_set<u64> reached;  // reached but not visited yet
        gtl::vector<u64> ready{sid};
        size_t len = 0;
        while (!ready.empty()) {
            const u64 cur = ready.back();
            ready.pop_back();
            auto cit = dbg.find(cur);
            if (cit == dbg.end() || cit->second.redirected()) return 0;
            const auto &cprev = cit->second.prev;
            if (std::find(cprev.begin(), cprev.end(), cur) != cprev.end()) {
                return 0;  // self-loop: a tandem repeat unit
            }
            visited.insert(cur);
            reached.erase(cur);
            if (cur != sid) {
                inner.push_back(cur);
                len += cit->second.seq.size();
                if (inner.size() > max_nodes || len > max_len) return 0;
            }
            auto nxt_it = fwd.find(cur);
            if (nxt_it == fwd.end() || nxt_it->second.empty()) {
                return 0;  // dead end inside the bubble
            }
            for (u64 next : nxt_it->second) {
                if (next == sid) return 0;  // cycle back onto the fork
                auto nit = dbg.find(next);
                if (nit == dbg.end() || nit->second.redirected()) {
                    return 0;  // stale edge from an earlier mutation
                }
                reached.insert(next);
                bool all_prev_visited = true;
                for (u64 p : nit->second.prev) {
                    auto pit = dbg.find(p);
                    if (pit == dbg.end() || pit->second.redirected()) {
                        continue;  // stale edge to a removed node
                    }
                    if (!visited.contains(p)) {
                        all_prev_visited = false;
                        break;
                    }
                }
                if (all_prev_visited) ready.push_back(next);
            }
            if (ready.size() == 1 && reached.size() == 1) {
                const u64 tid = ready[0];
                const auto &tprev = dbg.find(tid)->second.prev;
                if (std::find(tprev.begin(), tprev.end(), tid) != tprev.end()) {
                    return 0;  // self-loop on the sink
                }
                if (auto tnext = fwd.find(tid);
                    tnext != fwd.end() &&
                    std::find(tnext->second.begin(), tnext->second.end(),
                              sid) != tnext->second.end()) {
                    return 0;  // sink leads straight back to the fork
                }
                return inner.empty() ? 0 : tid;
            }
        }
        return 0;
    }

    // Heaviest sid -> tid path through a bubble found by find_bubble_sink
    // (`inner` lists every node after its predecessors). A node weighs its
    // lcmer support (core_count), or nothing when it is in `weightless`.
    // Returns the inner nodes of the path in order.
    static gtl::vector<u64> heaviest_bubble_path(
        const dbg_map_t &dbg, u64 sid, u64 tid, const gtl::vector<u64> &inner,
        const gtl::flat_hash_set<u64> &weightless) {
        // node -> (weight of the best path ending here, its predecessor)
        gtl::flat_hash_map<u64, std::pair<double, u64>> best;
        best[sid] = {0.0, 0};
        const auto best_into = [&](u64 nid, double weight) {
            std::pair<double, u64> out{-1.0, 0};
            for (u64 p : dbg.find(nid)->second.prev) {
                auto bit = best.find(p);
                if (bit == best.end()) continue;  // not in the bubble
                if (bit->second.first + weight > out.first) {
                    out = {bit->second.first + weight, p};
                }
            }
            return out;
        };
        for (u64 nid : inner) {
            const double weight =
                weightless.contains(nid)
                    ? 0.0
                    : (double)dbg.find(nid)->second.counts.core_count;
            best[nid] = best_into(nid, weight);
        }
        gtl::vector<u64> path;
        for (u64 nid = best_into(tid, 0.0).second; nid != sid;
             nid = best[nid].second) {
            path.push_back(nid);
        }
        std::reverse(path.begin(), path.end());
        return path;
    }

    // Copy of a bubble's removed allele with its depths and anchors.
    static popped_bubble make_popped_bubble(const dbg_map_t &dbg, u64 sid,
                                            u64 tid,
                                            const gtl::vector<u64> &alt_path,
                                            bool error, double removed_depth,
                                            double kept_depth) {
        gtl::vector<SequenceContainerType> frags;
        size_t len = 0;
        for (u64 nid : alt_path) {
            for (const auto &frag : dbg.find(nid)->second.seq) {
                frags.push_back(SequenceContainerType::copy_of(frag));
                len += frag.len;
            }
        }
        const auto &sseq = dbg.find(sid)->second.seq;
        const auto &tseq = dbg.find(tid)->second.seq;
        return popped_bubble{sseq.empty() ? 0 : sseq.back().origin_node,
                             tseq.empty() ? 0 : tseq.front().origin_node,
                             error,
                             (float)removed_depth,
                             (float)kept_depth,
                             frags.empty()
                                 ? SequenceContainerType()
                                 : SequenceContainerType(frags, 0, len)};
    }

    // Depth above which a node holds more copies than one locus of `ploidy`
    // copies, measured once before cleaning on the fork nodes (nodes with
    // two or more successors; most sit at one locus): median m and
    // cv = 1.4826 * MAD / m of the fork depths within half to one and a half
    // times ploidy * one_copy. A node beyond m * (1 + 3 cv) holds more copies
    // with about 1 in 1000 chance of being a locus node by noise. Measured
    // once because popping leaves mostly repeat forks behind, which would
    // shift the estimate up.
    static double locus_depth_limit(const dbg_map_t &dbg,
                                    const gtl::btree_set<u64> &active_set,
                                    double one_copy, int ploidy) {
        const double expected = ploidy * one_copy;
        const auto fwd = build_live_fwd(dbg, active_set);
        gtl::vector<double> depths;
        for (const auto &[id, succs] : fwd) {
            if (succs.size() < 2) continue;
            const double d = dbg.find(id)->second.avg_core_count();
            if (d >= 0.5 * expected && d <= 1.5 * expected) depths.push_back(d);
        }
        if (depths.size() < 100) return (ploidy + 0.5) * one_copy;
        std::sort(depths.begin(), depths.end());
        const double m = median_of_sorted(depths);
        gtl::vector<double> dev;
        dev.reserve(depths.size());
        for (double d : depths) dev.push_back(std::abs(d - m));
        std::sort(dev.begin(), dev.end());
        const double cv = 1.4826 * median_of_sorted(dev) / m;
        return m * (1.0 + 3.0 * cv);
    }

    // Bubble popping by copy number. `one_copy` is the depth of one copy (see
    // one_copy_depth), `locus_limit` the depth of one locus (see
    // locus_depth_limit). A bubble whose source, sink and inner nodes all
    // stay below locus_limit lies at one locus: its arms are alleles or read
    // errors, and it is popped, keeping its heaviest path. A bubble with any
    // node beyond it is a repeat whose arms are different copies; popping it
    // would collapse them, so it is kept. The removed side's depth only marks
    // the record as a read error (below half a copy). Removed nodes may hold
    // at most bubble_length_bound times the kept path's cores. Bubbles are
    // decided smallest first, so nested ones resolve from the inside out.
    // Each pop is appended to `popped` (under kKeepPoppedAlleles) before its
    // nodes go. max_bubble_len (inner cores) 0 disables popping.
    static size_t pop_bubbles(dbg_map_t &dbg, double one_copy,
                              double locus_limit, double bubble_length_bound,
                              size_t max_bubble_len, auto &active_set,
                              gtl::vector<popped_bubble> &popped) {
        if (max_bubble_len == 0 || one_copy <= 0.0) return 0;
        constexpr size_t kMaxBubbleNodes = 256;
        auto fwd = build_live_fwd(dbg, active_set);

        // Forks that open a bubble, smallest bubble first.
        gtl::vector<std::pair<size_t, u64>> forks;
        gtl::vector<u64> inner;
        for (u64 sid : active_set) {
            auto ffit = fwd.find(sid);
            if (ffit == fwd.end() || ffit->second.size() < 2) continue;
            if (find_bubble_sink(dbg, fwd, sid, kMaxBubbleNodes, max_bubble_len,
                                 inner) == 0) {
                continue;
            }
            size_t len = 0;
            for (u64 nid : inner) len += dbg.find(nid)->second.seq.size();
            forks.push_back({len, sid});
        }
        std::sort(forks.begin(), forks.end());

        gtl::vector<u64> to_rm;
        size_t n_popped = 0;
        for (const auto &[_, sid] : forks) {
            // Earlier pops change the graph, so the bubble is searched again.
            auto sit = dbg.find(sid);
            if (sit == dbg.end() || sit->second.redirected()) continue;
            auto ffit = fwd.find(sid);
            if (ffit == fwd.end() || ffit->second.size() < 2) continue;
            const u64 tid = find_bubble_sink(dbg, fwd, sid, kMaxBubbleNodes,
                                             max_bubble_len, inner);
            if (tid == 0) continue;

            const gtl::vector<u64> keep =
                heaviest_bubble_path(dbg, sid, tid, inner, {});
            const gtl::flat_hash_set<u64> kept(keep.begin(), keep.end());

            double removed_depth = 0.0, inner_high = 0.0;
            double kept_depth =
                std::min(sit->second.avg_core_count(),
                         dbg.find(tid)->second.avg_core_count());
            size_t removed_len = 0, kept_len = 0;
            for (u64 nid : inner) {
                const auto &node = dbg.find(nid)->second;
                const double d = node.avg_core_count();
                inner_high = std::max(inner_high, d);
                if (kept.contains(nid)) {
                    kept_len += node.seq.size();
                } else {
                    removed_depth = std::max(removed_depth, d);
                    removed_len += node.seq.size();
                }
            }
            for (u64 nid : keep) {
                kept_depth = std::min(kept_depth,
                                      dbg.find(nid)->second.avg_core_count());
            }
            if (removed_len == 0) continue;  // only skip edges
            if ((double)removed_len >
                bubble_length_bound * (double)std::max<size_t>(kept_len, 1)) {
                continue;  // removed side much longer than the kept path
            }
            const double end_depth =
                std::max(sit->second.avg_core_count(),
                         dbg.find(tid)->second.avg_core_count());
            if (end_depth >= locus_limit || inner_high >= locus_limit) {
                continue;  // a repeat: the arms are different copies
            }
            const bool error = removed_depth < 0.5 * one_copy;

            // The removed allele is the heaviest path over removed nodes.
            if constexpr (dbg_index::kKeepPoppedAlleles)
                popped.push_back(make_popped_bubble(
                    dbg, sid, tid,
                    heaviest_bubble_path(dbg, sid, tid, inner, kept), error,
                    removed_depth, kept_depth));

            // Delete removed nodes; along the kept path, drop every other
            // edge inside the bubble. Inner nodes have no edges outside it.
            gtl::flat_hash_set<u64> in_bubble(inner.begin(), inner.end());
            in_bubble.insert(sid);
            in_bubble.insert(tid);
            for (u64 nid : inner) {
                if (kept.contains(nid)) continue;
                dbg.erase(nid);
                fwd.erase(nid);
                to_rm.push_back(nid);
            }
            u64 from = sid;
            for (size_t i = 0; i <= keep.size(); ++i) {
                const u64 to = i < keep.size() ? keep[i] : tid;
                auto &pvec = dbg.find(to)->second.prev;
                pvec.erase(std::remove_if(pvec.begin(), pvec.end(),
                                          [&](u64 p) {
                                              return p != from &&
                                                     in_bubble.contains(p);
                                          }),
                           pvec.end());
                auto &next = fwd[from];
                next.erase(std::remove_if(next.begin(), next.end(),
                                          [&](u64 n) {
                                              return n != to &&
                                                     in_bubble.contains(n);
                                          }),
                           next.end());
                dbg.find(from)->second.forward_edge_count =
                    (unsigned)next.size();
                from = to;
            }
            ++n_popped;
        }
        for (u64 nid : to_rm) {
            active_set.erase(nid);
        }
        return n_popped;
    }

    // Simple chain leaving `start` away from `from`, in the forward
    // (successor) or backward (predecessor) direction, that dead-ends within
    // max_len cores. Fills `chain` and returns true for such a tip; returns
    // false when the chain continues, merges, forks, loops or is longer.
    static bool walk_tip(const dbg_map_t &dbg,
                         const gtl::flat_hash_map<u64, gtl::vector<u64>> &fwd,
                         u64 from, u64 start, bool forward, size_t max_len,
                         gtl::vector<u64> &chain) {
        chain.clear();
        size_t len = 0;
        for (u64 cur = start;;) {
            auto cit = dbg.find(cur);
            if (cit == dbg.end() || cit->second.redirected()) return false;
            // Edges into `cur` along the walk must come from the chain only.
            gtl::vector<u64> in, out;
            for (u64 p : cit->second.prev) {
                auto pit = dbg.find(p);
                if (pit != dbg.end() && !pit->second.redirected())
                    (forward ? in : out).push_back(p);
            }
            if (auto f = fwd.find(cur); f != fwd.end()) {
                for (u64 n : f->second) {
                    auto nit = dbg.find(n);
                    if (nit != dbg.end() && !nit->second.redirected())
                        (forward ? out : in).push_back(n);
                }
            }
            if (in.size() != 1) return false;  // merge into the chain
            chain.push_back(cur);
            len += cit->second.seq.size();
            if (len > max_len) return false;
            if (out.empty()) return true;  // dead end: a tip
            if (out.size() != 1 || out[0] == from || out[0] == start ||
                out[0] == cur) {
                return false;  // fork, loop or cycle
            }
            cur = out[0];
        }
    }

    // Fork tips by copy number. At a node with two or more branches on one
    // side, a branch that dead-ends within max_tip_len cores while another
    // branch continues is an allele or a read error when the fork and the
    // tip stay below locus_limit (see pop_bubbles); it is cut and recorded
    // (alt = the tip, anchored at the fork). Tips at deeper forks are repeat
    // copies and kept. The record marks tips below half a copy as errors.
    static size_t cut_fork_tips(dbg_map_t &dbg, double one_copy,
                                double locus_limit, size_t max_tip_len,
                                auto &active_set,
                                gtl::vector<popped_bubble> &popped) {
        if (one_copy <= 0.0) return 0;
        auto fwd = build_live_fwd(dbg, active_set);
        gtl::vector<u64> to_rm;
        gtl::vector<u64> chain;
        size_t n_cut = 0;
        for (u64 xid : active_set) {
            for (int forward = 1; forward >= 0; --forward) {
                auto xit = dbg.find(xid);
                if (xit == dbg.end() || xit->second.redirected()) break;
                if (xit->second.avg_core_count() >= locus_limit) break;
                gtl::vector<u64> branches;
                if (forward) {
                    if (auto f = fwd.find(xid); f != fwd.end())
                        branches = f->second;
                } else {
                    branches.assign(xit->second.prev.begin(),
                                    xit->second.prev.end());
                }
                branches.erase(
                    std::remove_if(branches.begin(), branches.end(),
                                   [&](u64 b) {
                                       auto bit = dbg.find(b);
                                       return bit == dbg.end() ||
                                              bit->second.redirected();
                                   }),
                    branches.end());
                if (branches.size() < 2) continue;
                gtl::vector<gtl::vector<u64>> tips;
                size_t continuing = 0;
                for (u64 b : branches) {
                    if (walk_tip(dbg, fwd, xid, b, forward, max_tip_len,
                                 chain)) {
                        tips.push_back(chain);
                    } else {
                        ++continuing;
                    }
                }
                if (continuing == 0) continue;
                for (const auto &tip : tips) {
                    double tip_depth = 0.0;
                    bool seeded = false;
                    for (u64 nid : tip) {
                        const auto &node = dbg.find(nid)->second;
                        tip_depth = std::max(tip_depth, node.avg_core_count());
                        seeded |= node.seeded;
                    }
                    if (tip_depth >= locus_limit) continue;  // repeat copy
                    if (seeded) continue;  // an earlier pass's segment
                    if constexpr (dbg_index::kKeepPoppedAlleles) {
                        // Tip sequence in forward order.
                        gtl::vector<SequenceContainerType> frags;
                        size_t len = 0;
                        for (size_t k = 0; k < tip.size(); ++k) {
                            const u64 nid =
                                forward ? tip[k] : tip[tip.size() - 1 - k];
                            for (const auto &frag : dbg.find(nid)->second.seq) {
                                frags.push_back(
                                    SequenceContainerType::copy_of(frag));
                                len += frag.len;
                            }
                        }
                        const auto &xseq = dbg.find(xid)->second.seq;
                        const u64 anchor = xseq.empty() ? 0
                                           : forward ? xseq.back().origin_node
                                                     : xseq.front().origin_node;
                        popped.push_back(popped_bubble{
                            forward ? anchor : 0, forward ? 0 : anchor,
                            tip_depth < 0.5 * one_copy, (float)tip_depth,
                            (float)xit->second.avg_core_count(),
                            frags.empty()
                                ? SequenceContainerType()
                                : SequenceContainerType(frags, 0, len)});
                    }
                    // Detach the tip from the fork, then delete it.
                    const u64 first = tip.front();
                    if (forward) {
                        auto &next = fwd[xid];
                        next.erase(std::remove(next.begin(), next.end(), first),
                                   next.end());
                        dbg.find(xid)->second.forward_edge_count =
                            (unsigned)next.size();
                    } else {
                        auto &pvec = dbg.find(xid)->second.prev;
                        pvec.erase(std::remove(pvec.begin(), pvec.end(), first),
                                   pvec.end());
                    }
                    for (u64 nid : tip) {
                        dbg.erase(nid);
                        fwd.erase(nid);
                        to_rm.push_back(nid);
                    }
                    ++n_cut;
                }
            }
        }
        for (u64 nid : to_rm) {
            active_set.erase(nid);
        }
        return n_cut;
    }

    static void cut_low_support_edges(dbg_map_t &dbg, auto &active_set) {
        for (auto &node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end() || current->second.redirected()) {
                // Skip consumed husks: a redirected node's seq is empty and
                // its active union member is redirect_to, not counts.
                continue;
            }
            if (current->second.avg_core_count() < 2 &&
                !current->second.seeded && current->second.prev.size() == 1 &&
                current->second.forward_edge_count == 1) {
                auto prev_it = dbg.find(current->second.prev[0]);
                if (prev_it != dbg.end())
                    prev_it->second.forward_edge_count--;
                dbg.erase(current);
            }
        }
    }
    // Drops prev edges to erased or redirected nodes. Every live node is in
    // active_set, so iterating it visits the same nodes as all of dbg.
    static void fix_edges(dbg_map_t &dbg, const auto &active_set) {
        for (u64 id : active_set) {
            auto nit = dbg.find(id);
            if (nit == dbg.end() || nit->second.redirected()) continue;
            auto &node = *nit;
            auto &pvec = node.second.prev;
            for (auto it = pvec.begin(); it != pvec.end(); ++it) {
                if (dbg.find(*it) == dbg.end() ||
                    dbg.find(*it)->second.redirected()) {
                    *it = std::numeric_limits<u64>::max();
                }
            }
            pvec.erase(std::remove_if(
                           pvec.begin(), pvec.end(),
                           [](u64 val) {
                               return val == std::numeric_limits<u64>::max();
                           }),
                       pvec.end());
        }
    }
    static void run(dbg_index &dbgi, const asm_config &acfg) {
        auto cfg = acfg.graph_cleaning;
        INIT_TIMERS(GraphCleaning);
        redir_map redirect_deleted;
        decltype(dbgi.active_set) &active_set = dbgi.active_set;

        if (cfg.iterations > -1) {
            // Initial simplification over all of dbg. active_set is not read
            // before the first prune, which would leave exactly the live
            // nodes, so it is built from those after the merge.
            auto to_simplify =
                simplify_simple_parallel_mark(dbgi.dbg, acfg.threads.work);
            simplify_simple_parallel_sweep(dbgi.dbg, to_simplify,
                                           redirect_deleted);
            build_active_set(dbgi.dbg, active_set, true);
            redirect_prev_edges(dbgi.dbg, redirect_deleted, active_set);
            IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                               "Initial Simplification");
            cut_low_support_edges(dbgi.dbg, active_set);
            const double one_copy = one_copy_depth(dbgi.dbg, active_set);
            const double locus_limit =
                one_copy > 0.0 ? locus_depth_limit(dbgi.dbg, active_set,
                                                   one_copy, cfg.ploidy)
                               : 0.0;
            dbgi.locus_limit = locus_limit;
            if (!acfg.silent)
                std::fprintf(stderr,
                             "[graph-clean] one-copy depth %.2f, ploidy %d, "
                             "locus depth limit %.2f\n",
                             one_copy, cfg.ploidy, locus_limit);
            // Isolated tips below half a copy are read errors; 3 without a
            // depth.
            const int tip_depth =
                one_copy > 0.0 ? std::max(1, (int)std::lround(0.5 * one_copy))
                               : 3;
            // An allele tip runs from its fork to the allele's lcmer hole,
            // about one lcmer window; 64 cores leave room for several. Only a
            // diploid sample has allele tips: on a haploid one the cut tips
            // were mostly satellite/repeat copies and added misassemblies.
            constexpr size_t kMaxForkTipCores = 64;
            const bool fork_tips = cfg.ploidy >= 2;
            // A round that removes and merges nothing leaves the graph as it
            // found it, apart from sorting prev lists and pruning erased
            // nodes from active_set, so the next round would change nothing.
            int rounds = 1;  // the initial simplification
            for (int i = 1; i < cfg.iterations; ++i) {
                ++rounds;
                fix_edges(dbgi.dbg, active_set);
                const size_t cut = cut_tips(dbgi.dbg, tip_depth,
                                            cfg.max_tip_length, active_set);
                const size_t popped = pop_bubbles(
                    dbgi.dbg, one_copy, locus_limit,
                    cfg.min_bubble_length_difference, cfg.max_bubble_length,
                    active_set, dbgi.popped_bubbles);
                const size_t tips =
                    fork_tips ? cut_fork_tips(dbgi.dbg, one_copy, locus_limit,
                                              kMaxForkTipCores, active_set,
                                              dbgi.popped_bubbles)
                              : 0;
                if (!acfg.silent)
                    std::fprintf(stderr,
                                 "[graph-clean] popped %zu bubbles, cut %zu "
                                 "fork tips\n",
                                 popped, tips);
                redirect_deleted.clear();
                // Iterations run simplify serially (single partition): the
                // parallel sweep has a cross-chain forward_edge_count race.
                simplify_simple_paths(dbgi.dbg, redirect_deleted, active_set,
                                      1);
                IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                                   "Iteration %d done.", i);
                if (cut + popped + tips + redirect_deleted.size() == 0) break;
            }
            fix_edges(dbgi.dbg, active_set);
            if (!acfg.silent)
                std::fprintf(stderr, "[graph-clean] %d cleaning rounds\n",
                             rounds);
        } else {
            build_active_set(dbgi.dbg, active_set, false);
        }
    }
};

__attribute__((format(printf, 2, 3))) static inline void vec_appendf(
    std::vector<char> &v, const char *fmt, ...) {
    va_list args, args_copy;
    va_start(args, fmt);
    va_copy(args_copy, args);
    int n = vsnprintf(nullptr, 0, fmt, args);
    va_end(args);
    if (n > 0) {
        size_t old = v.size();
        v.resize(old + static_cast<size_t>(n) + 1);
        vsnprintf(v.data() + old, static_cast<size_t>(n) + 1, fmt, args_copy);
        v.pop_back();  // drop the null terminator
    }
    va_end(args_copy);
}

inline void print_dbg_seqs(FILE *out, const dbg_index &dbgi, int nthreads = 0) {
    if (nthreads <= 0)
        nthreads = std::max(1u, std::thread::hardware_concurrency()) - 1;
    if (nthreads < 1) nthreads = 1;

    // Work item: a pointer into the dbg map (stable while we only read).
    struct WorkItem {
        uint64_t lval;
        const typename dbg_index::bidbg_val *edges;
    };

    // Items in flight per worker before the producer blocks.
    constexpr size_t kItemsPerWorker = 1024;
    bounded_queue<WorkItem> work_queue{
        kItemsPerWorker * static_cast<size_t>(nthreads), /*producers=*/1};
    std::atomic<bool> work_done{false};
    std::atomic<int> active_workers{nthreads};
    // ── Worker thread pool ──────────────────────────────────────────
    // Each worker dequeues nodes, formats the S-line + L-lines into a
    // vector<char> and writes it to `out`.
    // Format one node (S-line + L-lines) into line_buf.
    auto format_node = [&](std::vector<char> &line_buf, const WorkItem &item) {
        const auto &edges = *item.edges;
        const uint64_t lval = item.lval;
        line_buf.clear();

        // Pre-compute total sequence length.
        size_t sofar = 0;
        for (const auto &s : edges.seq) {
            sofar += s.len;
        }

        // S-line header.
        vec_appendf(line_buf, "S\tC%" PRIu64 "\t", lval);

        // Convert all sequence chunks directly into line_buf.
        {
            size_t seq_start = line_buf.size();
            line_buf.resize(seq_start + sofar + 1);  // +1 for trailing \0
            size_t offset = 0;
            for (const auto &s : edges.seq) {
                s.convert_to_ascii(line_buf.data() + seq_start + offset);
                offset += s.len;
            }
            line_buf.resize(seq_start + sofar);  // trim null terminator
        }

        // S-line trailer.
        vec_appendf(line_buf,
                    "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n", sofar,
                    (int)(edges.avg_core_count() * sofar),
                    edges.counts.core_count, edges.prev.size(),
                    edges.forward_edge_count);

        // L-lines for incoming edges.
        for (uint64_t nl : edges.prev) {
            if (dbgi.dbg.contains(nl))
                vec_appendf(line_buf,
                            "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", nl,
                            lval, 0);
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        workers.emplace_back([&] {
            WorkItem item;
            std::vector<char> line_buf;
            line_buf.reserve(1 << 16);

            while (true) {
                if (!work_queue.wait_dequeue_timed(item, WAIT_QUEUE_TIMEOUT)) {
                    if (work_done.load()) break;
                    continue;
                }
                format_node(line_buf, item);
                fwrite(line_buf.data(), line_buf.size(), 1, out);
            }

            // Drain remaining work after done flag is set.
            while (work_queue.try_dequeue(item)) {
                format_node(line_buf, item);
                fwrite(line_buf.data(), line_buf.size(), 1, out);
            }
            active_workers.fetch_sub(1);
        });
    }

    // ── Producer: main thread enqueues work ─────────────────────────
    for (const auto &[lval, edges] : dbgi.dbg) {
        if (edges.counts.core_count == 0 && edges.prev.size() == 0) continue;
        if (edges.seq.size() == 0) continue;

        work_queue.enqueue(WorkItem{lval, &edges});
    }
    work_done.store(true);

    for (auto &w : workers) w.join();
    // All workers have exited and decremented active_workers to 0.
}

#endif  // _GRAPH_H_
