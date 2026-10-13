#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "gtl/btree.hpp"
#define PREV_WITH_VEC
#define HT_SHARD_BC 8

#include <blockingconcurrentqueue.h>

#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <limits>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>
#include <cinttypes>

#include "asm.h"
#include "bloomfilter.h"
#include "config.h"
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
// A bubble removed by GraphCleaning::pop_bubbles, kept so the removed allele
// can be written as an alternate sequence after layout. The fragments of the
// source's last and the sink's first lcmer keep their origin ids through
// every later merge, so they locate the bubble in the final unitigs.
struct popped_bubble {
    u64 source_origin;
    u64 sink_origin;
    bool error;                 // removed nodes read below half a copy
    float removed_depth;        // highest avg_core_count among removed nodes
    float kept_depth;           // lowest avg_core_count on the kept path
    SequenceContainerType alt;  // best removed path, source to sink
};
template <class BFT = DefaultBF>
struct dbg_index_base {
    using BloomFilterType = BFT;
    static constexpr u64 lcindex_increment =
        nextpowerof2((u64)MAXIMUM_WORK_THREADS);
    using mutex_type = std::mutex;
    struct bidbg_val {
        struct counts_t {
            unsigned core_count;
        };

        bidbg_val(unsigned core_count)
            : redirected_flag(false), forward_edge_count(0), counts{core_count},
              seq{}, prev{} {}
        bidbg_val(unsigned core_count, u64 prev)
            : redirected_flag(false), forward_edge_count(0), counts{core_count},
              seq{}, prev{prev} {}

        bool redirected_flag;
        // forward_edge_count stays outside the redirect_to union: a
        // neighbouring chain may read it while this node is being redirected.
        // fec is written only on representatives; core_count stays unioned
        // since only the node's own chain touches it.
        unsigned forward_edge_count;
        union {
            counts_t counts;
            u64 redirect_to;
        };

        gtl::vector<SequenceContainerType> seq;

        bool redirected() const {
            return redirected_flag || (prev.size() == 1 && prev[0] == 0);
        }
        void mark_redirected(u64 to) {
            redirected_flag = true;
            redirect_to = to;
            gtl::vector<SequenceContainerType>().swap(seq);
            gtl::vector<u64>().swap(prev);
            prev.push_back(0);  // 0 is an illegal node id; marks redirected.
        }

        gtl::vector<u64> prev;
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
    using dbg_map_t = gtl::parallel_flat_hash_map<
        u64, bidbg_val, std::hash<u64>, gtl::priv::hash_default_eq<u64>,
        std::allocator<std::pair<const u64, bidbg_val>>, HT_SHARD_BC,
        mutex_type>;

    struct lci_type {
        u64 id;
        // u64 moved : 1;
        static constexpr size_t max_id = std::numeric_limits<u64>::max();
        static consteval lci_type empty() { return {0}; }
    };

    dbg_map_t dbg;
    gtl::btree_set<u64> active_set;
    gtl::vector<popped_bubble> popped_bubbles;

    // merge_info_map_t merge_map;
    BloomFilterType bf;

    dbg_index_base(double false_positive_rate, uint64_t projected_count,
                   uint64_t seed)
        : bf{false_positive_rate, projected_count, seed} {}

    // dbg_index(const dbg_index &other)
    static constexpr u64 lcmer_id(const LCmer &l) { return l.value(); }

    lci_type get(const LCmer &l) const {
        u64 id = lcmer_id(l);
        gtl::LockableBaseImpl<gtl::NullMutex>::DoNothing no_lock;
        auto &dbg_nc = const_cast<dbg_map_t &>(dbg);
        auto h = dbg_nc.hash(id);
        auto it = dbg_nc.find(id, h, no_lock);
        if (it == dbg_nc.end()) {
            return lci_type::empty();
        }
        for (int hop = 0;
             hop < 128 && it != dbg_nc.end() && it->second.redirected();
             ++hop) {
            const u64 next = it->second.redirect_to;
            if (next == 0 || next == id) {
                return lci_type::empty();
            }
            id = next;
            h = dbg_nc.hash(id);
            it = dbg_nc.find(id, h, no_lock);
        }
        return it == dbg_nc.end() ? lci_type::empty() : lci_type{id};
    }
    void count_forward_edges() {
        for (auto &node : dbg) {
            for (auto &edge : node.second.prev) {
                auto it = dbg.find(edge);
                if (it == dbg.end()) {
                    continue;
                }
                it->second.forward_edge_count++;
            }
        }
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
        if (!bf.contains(new_hash)) {
            bf.insert(new_hash);
            return 0;
        }
        return lcmer_id(l);
    }
    gtl::vector<u64> shortest_path_bfs(u64 source, u64 target) const {
        if (source == target) return gtl::vector<u64>{source};
        if (source == 0 || target == 0) return {};

        std::queue<u64> q;
        gtl::flat_hash_map<u64, u64> parent;
        gtl::flat_hash_map<u64, uint8_t> visited;

        visited[source] = 1;
        q.push(source);

        while (!q.empty()) {
            const u64 cur = q.front();
            q.pop();

            auto it = dbg.find(cur);
            if (it == dbg.end()) {
                continue;
            }
            for (u64 nb : it->second.prev) {
                if (visited.find(nb) != visited.end()) continue;
                auto neit = dbg.find(nb);
                if (neit == dbg.end()) continue;

                visited.emplace(nb, 1);
                parent.emplace(nb, cur);

                if (nb == target) {
                    gtl::vector<u64> path;
                    u64 v = target;
                    path.push_back(v);
                    while (v != source) {
                        auto pit = parent.find(v);
                        if (pit == parent.end()) return {};
                        v = pit->second;
                        path.push_back(v);
                    }
                    std::reverse(path.begin(), path.end());
                    return path;
                }

                q.push(nb);
            }
        }
        return {};
    }
};

using dbg_index = dbg_index_base<DefaultBF>;
using dbg_index_light = dbg_index_base<MockBF>;

struct FlatDBG {
    struct Node {
        gtl::vector<u64> prev;
        gtl::vector<u64> next;
        gtl::vector<SequenceContainerType> seq;
    };
    gtl::vector<Node> nodes;
};
struct GraphCleaning {
    using bidbg_val = dbg_index::bidbg_val;
    using lci_type = dbg_index::lci_type;
    using dbg_map_t = dbg_index::dbg_map_t;
    using redir_map = gtl::parallel_flat_hash_map_d<u64, u64, 4, std::mutex>;
    struct avg_core_count_stats {
        double baseline = 0.0;
    };
    static double median_of_sorted(const gtl::vector<double> &vals) {
        if (vals.empty()) return 0.0;
        const size_t mid = vals.size() / 2;
        if (vals.size() % 2 == 1) return vals[mid];
        return (vals[mid - 1] + vals[mid]) / 2.0;
    }

    static avg_core_count_stats collect_active_avg_core_count_stats(
        const dbg_map_t &dbg, const gtl::btree_set<u64> &active_set) {
        avg_core_count_stats stats;
        constexpr size_t kMinDepthSample = 16;
        struct depth_len_sample {
            double depth;
            size_t len_proxy;
            bool topology_stable;
        };
        gtl::vector<depth_len_sample> samples;
        samples.reserve(active_set.size());
        gtl::vector<size_t> lengths;
        lengths.reserve(active_set.size());

        for (u64 node : active_set) {
            auto it = dbg.find(node);
            if (it == dbg.end() || it->second.redirected()) continue;
            if (it->second.seq.empty()) continue;
            const double d = it->second.avg_core_count();
            if (d <= 0.0) continue;

            const size_t len_proxy = it->second.seq.size();
            const bool topology_stable =
                (it->second.prev.size() <= 1 &&
                 it->second.forward_edge_count <= 1);
            samples.push_back({d, len_proxy, topology_stable});
            lengths.push_back(len_proxy);
        }
        if (samples.empty()) return stats;

        std::sort(lengths.begin(), lengths.end());
        const size_t len_q75 = lengths[(lengths.size() - 1) * 3 / 4];

        auto collect_depths = [&](auto pred) {
            gtl::vector<double> out;
            out.reserve(samples.size());
            for (const auto &s : samples) {
                if (pred(s)) out.push_back(s.depth);
            }
            return out;
        };

        gtl::vector<double> depths = collect_depths([&](const auto &s) {
            return s.topology_stable && s.len_proxy >= len_q75;
        });
        if (depths.size() < kMinDepthSample) {
            depths = collect_depths(
                [&](const auto &s) { return s.len_proxy >= len_q75; });
        }
        if (depths.size() < kMinDepthSample) {
            depths = collect_depths([](const auto &) { return true; });
        }

        std::sort(depths.begin(), depths.end());
        const auto qidx = [n = depths.size()](double q) -> size_t {
            const size_t idx =
                static_cast<size_t>(q * static_cast<double>(n - 1));
            return std::min(idx, n - 1);
        };
        const double raw_q1 = depths[qidx(0.25)];
        const double raw_q3 = depths[qidx(0.75)];
        const double iqr = raw_q3 - raw_q1;
        const double lower = raw_q1 - 1.5 * iqr;
        const double upper = raw_q3 + 1.5 * iqr;

        gtl::vector<double> filtered;
        filtered.reserve(depths.size());
        for (double d : depths) {
            if (d >= lower && d <= upper) filtered.push_back(d);
        }
        const auto &pool = filtered.empty() ? depths : filtered;
        stats.baseline = median_of_sorted(pool);
        return stats;
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
    struct loop_copy_estimate {
        double base_cn = 1.0;
        int base_cn_rounded = 1;
        int other_cn_rounded = 1;
        double global_ref_depth = 0.0;
        double raw_loop_cn = 1.0;
        double raw_curr_cn = 2.0;
        int loop_cn_rounded = 1;
        int curr_cn_rounded = 2;
        int loop_turns = 1;
        double turn_event_gain = 0.0;
        double turn_event_threshold = 0.0;
        int turn_forced_single = 0;
        double score = std::numeric_limits<double>::infinity();
    };

    static loop_copy_estimate infer_loop_copy_numbers(
        double other_in_depth, double loop_depth, double curr_depth,
        const avg_core_count_stats &depth_stats) {
        assert(depth_stats.baseline > 0.0);
        constexpr int kMaxSearchTurns = 64;
        constexpr int kMaxSearchBaseCn = 64;
        constexpr double kBaseCnPriorWeight = 0.35;
        constexpr double kTurnGainBaseThreshold = 0.25;

        const double baseline = depth_stats.baseline;
        const double other_base_cn = std::max(1.0, other_in_depth / baseline);

        loop_copy_estimate est;
        est.global_ref_depth = baseline;
        est.base_cn = other_base_cn;
        est.base_cn_rounded =
            std::max(1, static_cast<int>(std::llround(est.base_cn)));
        est.other_cn_rounded = est.base_cn_rounded;
        est.raw_loop_cn = loop_depth / baseline;
        est.raw_curr_cn = curr_depth / baseline;
        est.loop_cn_rounded =
            std::max(1, static_cast<int>(std::llround(est.raw_loop_cn)));
        est.curr_cn_rounded =
            std::max(1, static_cast<int>(std::llround(est.raw_curr_cn)));

        double best_err = std::numeric_limits<double>::infinity();
        int best_t = 1;
        int best_b = est.base_cn_rounded;
        double best_err_t1 = std::numeric_limits<double>::infinity();
        int best_b_t1 = est.base_cn_rounded;
        const auto fit_err = [&](int b, int t) {
            const double exp_loop = static_cast<double>(t * b);
            const double exp_curr = static_cast<double>((t + 1) * b);
            return std::abs(est.raw_loop_cn - exp_loop) +
                   std::abs(est.raw_curr_cn - exp_curr) +
                   kBaseCnPriorWeight *
                       std::abs(est.base_cn - static_cast<double>(b));
        };
        for (int b = 1; b <= kMaxSearchBaseCn; ++b) {
            for (int t = 1; t <= kMaxSearchTurns; ++t) {
                const double err = fit_err(b, t);
                if (err < best_err ||
                    (err == best_err &&
                     (t < best_t || (t == best_t && b < best_b)))) {
                    best_err = err;
                    best_t = t;
                    best_b = b;
                }
                if (t == 1 && err < best_err_t1) {
                    best_err_t1 = err;
                    best_b_t1 = b;
                }
            }
        }
        const double turn_gain = best_err_t1 - best_err;
        const double turn_threshold =
            kTurnGainBaseThreshold * static_cast<double>(est.other_cn_rounded);
        if (best_t > 1 && turn_gain < turn_threshold) {
            best_t = 1;
            best_b = best_b_t1;
            best_err = best_err_t1;
            est.turn_forced_single = 1;
        }

        est.base_cn_rounded = best_b;
        est.loop_turns = best_t;
        est.turn_event_gain = turn_gain;
        est.turn_event_threshold = turn_threshold;
        est.score = best_err;
        return est;
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
        for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
            // for (auto &seq : (*it)->second.seq)
            gtl::vector<SequenceContainerType> &seqs = (*it)->second.seq;
            for (auto sit = seqs.begin(); sit != seqs.end(); ++sit) {
                current->second.seq.emplace_back(std::move(*sit));
            }
            seqs.clear();
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
                // printf("TTTT:%d\n", tid);
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
    // [0] and 0 is never a node id. Every live node is in active_set (nodes
    // only ever leave it), so iterating it instead of all of dbg (mostly
    // redirected husks) gives the same result; each node edits only itself.
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

    static void cut_tips(dbg_map_t &dbg, int min_core_depth,
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

        // for (auto &dn : lci) {
        for (u64 node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end()) {
                continue;
            }
            if (current->second.redirected()) {
                continue;
            }
            if (current->second.avg_core_count() < min_core_depth) {
                // Forlorn
                if (current->second.forward_edge_count == 0 &&
                    current->second.prev.size() == 0) {
                    dbg.erase(current);  // MARK
                } else if (current->second.forward_edge_count == 0 &&
                           current->second.prev.size() == 1) {
                    auto prev_it = dbg.find(current->second.prev[0]);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        prev_it->second.forward_edge_count--;
                    }
                    dbg.erase(current);  // MARK
                }
            } else {
                gtl::vector<u64> prev_ids(current->second.prev.begin(),
                                          current->second.prev.end());
                std::sort(prev_ids.begin(), prev_ids.end());
                for (u64 pid : prev_ids) {
                    auto prev_it = dbg.find(pid);
                    if (prev_it != dbg.end() && !prev_it->second.redirected()) {
                        if (prev_it->second.avg_core_count() < min_core_depth &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);  // MARK
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
                            prev_it->second.seq.size() < max_tip_length &&
                            current->second.seq.size() > 2 * max_tip_length &&
                            prev_it->second.prev.size() == 0 &&
                            prev_it->second.forward_edge_count == 1) {
                            dbg.erase(prev_it);  // MARK
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
    }

    // Live forward adjacency: successor lists derived from prev edges,
    // sorted per node for deterministic arm walking. Iterates the live
    // active_set, not the whole dbg map.
    static gtl::flat_hash_map<u64, gtl::vector<u64>> build_live_fwd(
        const dbg_map_t &dbg, const gtl::btree_set<u64> &active_set) {
        gtl::flat_hash_map<u64, gtl::vector<u64>> fwd;
        fwd.reserve(active_set.size());
        for (u64 id : active_set) {
            auto it = dbg.find(id);
            if (it == dbg.end() || it->second.redirected()) continue;
            for (u64 p : it->second.prev) {
                auto pit = dbg.find(p);
                if (pit != dbg.end() && !pit->second.redirected()) {
                    fwd[p].push_back(id);
                }
            }
        }
        for (auto &[_, succs] : fwd) {
            std::sort(succs.begin(), succs.end());
        }
        return fwd;
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
    // Each pop is appended to `popped` before its nodes go. max_bubble_len
    // (inner cores) 0 disables popping.
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
            popped.push_back(make_popped_bubble(
                dbg, sid, tid, heaviest_bubble_path(dbg, sid, tid, inner, kept),
                error, removed_depth, kept_depth));

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
                    branches = xit->second.prev;
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
                    gtl::vector<SequenceContainerType> frags;
                    size_t len = 0;
                    for (u64 nid : tip) {
                        const auto &node = dbg.find(nid)->second;
                        tip_depth = std::max(tip_depth, node.avg_core_count());
                    }
                    if (tip_depth >= locus_limit) continue;  // repeat copy
                    // Tip sequence in forward order.
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
                                       : forward    ? xseq.back().origin_node
                                                    : xseq.front().origin_node;
                    popped.push_back(popped_bubble{
                        forward ? anchor : 0, forward ? 0 : anchor,
                        tip_depth < 0.5 * one_copy, (float)tip_depth,
                        (float)xit->second.avg_core_count(),
                        frags.empty() ? SequenceContainerType()
                                      : SequenceContainerType(frags, 0, len)});
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
                current->second.prev.size() == 1 &&
                current->second.forward_edge_count == 1) {
                auto prev_it = dbg.find(current->second.prev[0]);
                if (prev_it != dbg.end())
                    prev_it->second.forward_edge_count--;
                dbg.erase(current);
            }
        }
    }
    static void open_simple_loopsv2(dbg_map_t &dbg,
                                    gtl::btree_set<u64> &active_set,
                                    const avg_core_count_stats &depth_stats) {
        gtl::vector<u64> to_rm;
        for (auto &node : active_set) {
            auto current = dbg.find(node);
            if (current == dbg.end() || current->second.redirected()) {
                continue;
            }
            if (current->second.forward_edge_count == 2 &&
                current->second.prev.size() == 2) {
                // 2 in 2 out
                //
                //      o----->[-Cur-]----->o
                //             |     ^
                //             v     |
                //             [-Pre-]
                //
                // for (auto &pnode : current->second.prev) {
                int merged_node_index_in_prev = [&]() -> int {
                    for (size_t i = 0; i < current->second.prev.size(); ++i) {
                        auto pnode = current->second.prev[i];
                        const auto &pit = dbg.find(pnode);
                        if (pit == dbg.end() || pit->second.redirected()) {
                            continue;
                        }
                        if (pit->second.forward_edge_count != 1 ||
                            pit->second.prev.size() != 1) {  // Prev is simple
                            continue;
                        }
                        if (pit->second.prev[0] == node) {  // Loop confirmed.
                            // Try to find other Since we confirmed 2 node
                            // exist, we can flip it for other
                            auto other_prev = current->second.prev[1 - i];
                            double curr_depth =
                                current->second.avg_core_count();
                            double pr_depth = pit->second.avg_core_count();
                            assert(depth_stats.baseline > 0.0);
                            double other_depth = depth_stats.baseline;
                            if (auto other_it = dbg.find(other_prev);
                                other_it != dbg.end() &&
                                !other_it->second.redirected()) {
                                const double d =
                                    other_it->second.avg_core_count();
                                if (d > 0.0) other_depth = d;
                            }
                            // Lets combine the sequences here.
                            SequenceContainerType cur_seq{
                                current->second.seq.begin(),
                                current->second.seq.end()};
                            SequenceContainerType pre_seq{
                                pit->second.seq.begin(), pit->second.seq.end()};

                            const auto copy_est = infer_loop_copy_numbers(
                                other_depth, pr_depth, curr_depth, depth_stats);
                            const int extra_turns =
                                std::max(0, copy_est.loop_turns - 1);

                            for (int turn = 0; turn < extra_turns; ++turn) {
                                current->second.seq.push_back(
                                    SequenceContainerType::copy_of(pre_seq));
                                current->second.seq.push_back(
                                    SequenceContainerType::copy_of(cur_seq));
                            }
                            current->second.seq.push_back(std::move(pre_seq));
                            current->second.seq.push_back(std::move(cur_seq));

                            // Get rid of the looping node
                            return i;
                        }
                    }
                    return -1;
                }();
                if (merged_node_index_in_prev == -1) {
                    continue;
                }
                to_rm.push_back(
                    current->second.prev[merged_node_index_in_prev]);
                // active_set.erase(current->second.prev[merged_node_index_in_prev]);
                current->second.forward_edge_count = 1;
                if (merged_node_index_in_prev == 0) {
                    current->second.prev[0] = current->second.prev[1];
                }
                current->second.prev.pop_back();
            }
        }
        for (auto node : to_rm) {
            active_set.erase(node);
            dbg.erase(node);
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
            for (int i = 1; i < cfg.iterations; ++i) {
                fix_edges(dbgi.dbg, active_set);
                cut_tips(dbgi.dbg, tip_depth, cfg.max_tip_length, active_set);
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
                // erase_if(dbgi.lci, [](const auto &node) {  // MARK
                //     return node.second.id == lci_type::max_id;
                // });
                // Iterations run simplify serially (single partition): the
                // parallel sweep has a cross-chain forward_edge_count race.
                simplify_simple_paths(dbgi.dbg, redirect_deleted, active_set, 1);
                IF_TIME_CHECKPOINT(!acfg.silent, GraphCleaning, stderr,
                                   "Iteration %d done.", i);
            }
            fix_edges(dbgi.dbg, active_set);
        } else {
            build_active_set(dbgi.dbg, active_set, false);
        }
    }
};

inline void print_raw_dbg(FILE *out, const dbg_index &dbgi) {
    for (const auto &[lval, edges] : dbgi.dbg) {
        fprintf(out, "S\tC%" PRIu64 "\t%s\tLN:i:%d\tRC:i:%d", lval, "*", 100,
                edges.counts.core_count);
    }

    for (const auto &[lval, edges] : dbgi.dbg) {
        for (const auto &p : edges.prev) {
            fprintf(out, "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", p, lval,
                    1);
        }
    }
}

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

    using work_q = moodycamel::BlockingConcurrentQueue<WorkItem>;
    work_q work_queue;
    std::atomic<bool> work_done{false};
    std::atomic<int> active_workers{nthreads};
    // ── Worker thread pool ──────────────────────────────────────────
    // Each worker dequeues nodes, formats the S-line + L-lines into a
    // vector<char>, then enqueues it.
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
                const auto &edges = *item.edges;
                const uint64_t lval = item.lval;
                line_buf.clear();

                size_t sofar = 0;
                for (const auto &s : edges.seq) {
                    sofar += s.len;
                }

                vec_appendf(line_buf, "S\tC%" PRIu64 "\t", lval);

                {
                    size_t seq_start = line_buf.size();
                    line_buf.resize(seq_start + sofar + 1);
                    size_t offset = 0;
                    for (const auto &s : edges.seq) {
                        s.convert_to_ascii(line_buf.data() + seq_start +
                                           offset);
                        offset += s.len;
                    }
                    line_buf.resize(seq_start + sofar);
                }

                vec_appendf(line_buf,
                            "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n",
                            sofar, (int)(edges.avg_core_count() * sofar),
                            edges.counts.core_count, edges.prev.size(),
                            edges.forward_edge_count);

                for (uint64_t nl : edges.prev) {
                    if (dbgi.dbg.contains(nl))
                        vec_appendf(line_buf,
                                    "L\tC%" PRIu64 "\t+\tC%" PRIu64
                                    "\t+\t%dM\n",
                                    nl, lval, 0);
                }
                fwrite(line_buf.data(), line_buf.size(), 1, out);

                line_buf.clear();
            }
            active_workers.fetch_sub(1);
        });
    }

    // ── Producer: main thread enqueues work ─────────────────────────
    for (const auto &[lval, edges] : dbgi.dbg) {
        if (edges.counts.core_count == 0 && edges.prev.size() == 0) continue;
        if (edges.seq.size() == 0) continue;

        while (!work_queue.try_enqueue(WorkItem{lval, &edges})) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(10));
        }
    }
    work_done.store(true);

    for (auto &w : workers) w.join();
    // All workers have exited and decremented active_workers to 0;
    // the writer thread will notice and drain the output queue.
}

inline void print_dbg_dot(const dbg_index &dbgi, int i) {
    std::string filename = std::string{"graph"} + std::to_string(i) + ".dot";
    FILE *dot_file = fopen(filename.c_str(), "w");
    if (dot_file) {
        // print_raw_dbg(dot_file, dbgi);
        fprintf(dot_file, "digraph G {\n");
        // Legend
        fprintf(dot_file, "\tlegend [shape=none, label=<\n");
        fprintf(
            dot_file,
            "\t\t<table border=\"0\" cellborder=\"1\" cellspacing=\"0\">\n");
        fprintf(dot_file,
                "\t\t\t<tr><td bgcolor=\"lightgrey\">Node ID</td><td "
                "bgcolor=\"lightgrey\">Forward</td><td "
                "bgcolor=\"lightgrey\">Prev</td><td "
                "bgcolor=\"lightgrey\">Merged count</td><td "
                "bgcolor=\"lightgrey\">Core Count </td></tr>\n");
        fprintf(dot_file, "\t\t</table>\n");
        fprintf(dot_file, "\t>];\n");
        for (auto &node : dbgi.dbg) {
            fprintf(dot_file,
                    "\t%" PRIu64 " [label=\"%" PRIu64 "->%u,%zu,%zu,%d\"];\n",
                    node.first, node.first,
                    node.second.forward_edge_count,
                    node.second.prev.size(), node.second.seq.size(),
                    node.second.counts.core_count);
            for (auto &prev : node.second.prev) {
                fprintf(dot_file, "\t%" PRIu64 " -> %" PRIu64 ";\n", prev,
                        node.first);
            }
        }
        fprintf(dot_file, "}\n");
        fclose(dot_file);
    }
}

#endif  // _GRAPH_H_
