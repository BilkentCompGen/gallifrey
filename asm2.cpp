
// #include <concurrentqueue.h>
#include <blockingconcurrentqueue.h>
#include <fcntl.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <inttypes.h>
#include <stdio.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gtl/phmap.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>

#include "asm.h"
#include "bgzf.h"
#include "bloomfilter.h"
#include "config.h"
#include "ext/doctest.h"
#include "lcmer.h"
#include "lcptools/lcptools_ho.h"
#include "nibble2base.hpp"
#include "savebuffer_bh.h"
#include "syncmer.h"

using u64 = uint64_t;
using u32 = uint32_t;

using std::string;
enum class input_type { SAM, BAM, FASTX, FASTXGZ, FASTXBZ2, UNKNOWN };

template <std::integral int_type>
consteval int_type nextpowerof2(int_type v) {
    v--;
    for (size_t i = 1; i < sizeof(int_type) * 8; i <<= 1) {
        v |= v >> i;
    }
    v++;
    return v;
}

#define PREV_WITH_VEC
#define HT_SHARD_BC 8
struct dbg_index {
    static constexpr u64 lcindex_increment =
        nextpowerof2((u64)MAXIMUM_WORK_THREADS);
    using mutex_type = std::mutex;
    struct bidbg_val {
        unsigned cc;
        unsigned fec;
        gtl::vector<bam_nibble_data> seq;
#ifdef PREV_WITH_SET
        gtl::flat_hash_set<u64> prev;
        void emplace_prev(u64 val) { prev.emplace(val); }
#else
#ifdef PREV_WITH_VEC
        gtl::vector<u64> prev;

        void emplace_prev(u64 val) {
            auto it = std::find(prev.begin(), prev.end(), val);
            if (it == prev.end()) {
                prev.push_back(val);
            }
        }
#endif
#endif
        constexpr double avg_core_count() const {
            return (double)cc / seq.size();
        }
    };
    using dbg_map_t = gtl::parallel_flat_hash_map<
        u64, bidbg_val, std::hash<u64>, gtl::priv::hash_default_eq<u64>,
        std::allocator<std::pair<const u64, bidbg_val>>, HT_SHARD_BC,
        mutex_type>;

    struct lci_type {
        u64 id;
        //    u64 pos; // TODO remove this
        // lci_type() : id(0), pos(-1) {}
        // lci_type(u64 _id, u64 _pos) : id(_id), pos(_pos) {}
        // lci_type(u64 _id) : id(_id), pos(-1) {}
    };
    using lcmer_index_t = gtl::parallel_flat_hash_map<
        lcmer<LCMER_L>, lci_type, std::hash<lcmer<LCMER_L>>,
        gtl::priv::hash_default_eq<lcmer<LCMER_L>>,
        std::allocator<std::pair<const lcmer<LCMER_L>, lci_type>>, HT_SHARD_BC,
        mutex_type>;

    dbg_map_t dbg;
    lcmer_index_t lci;
    // merge_info_map_t merge_map;
    bloom_interface<bloom_ext::RegisterBlockedBloomFilter<1>, null_mutex> bf;

    dbg_index(double false_positive_rate, uint64_t projected_count,
              uint64_t seed)
        : bf{false_positive_rate, projected_count, seed} {}
    // dbg_index(const dbg_index &other)
    u64 get(const lcmer<LCMER_L> &l) const {
        auto it = lci.find(l);
        return it == lci.end() ? 0 : it->second.id;
        // if (lci.contains(l)) {
        //     return lci.at(l);
        // }
        // // lci.if_contains(l,[&ret](const lcmer_index_t::value_type &v){ ret
        // = v.second;}); return 0;
    }
    void count_forward_edges() {
        for (auto &node : dbg) {
            for (auto &edge : node.second.prev) {
                auto it = dbg.find(edge);
                if (it == dbg.end()) {
                    continue;
                }
                it->second.fec++;
            }
        }
    }
    void simplify_graph(int max_iterations = 1, size_t max_tip_length = 5,
                        size_t bubble_diff_eps = 2, size_t min_core_depth = 3,
                        double bubble_depth_epsilon = 2) {
        gtl::flat_hash_map<u64, u64> redirect_deleted;
        auto try_simplify_simple_path = [&](u64 id, auto current) -> bool {
            // Find a node we cannot simplify
            // Iterate backwards and merge until we cannot simplify
            // o->o->o->o->o->o -> o->o

            if (current == dbg.end()) {
                return false;
            }
            unsigned first_fec = current->second.fec;
            bool simplified = false;
            auto prev_it = dbg.find(current->second.prev[0]);
            gtl::vector<bam_nibble_data> merged;
            unsigned cc_sum = 0;
            while (prev_it != dbg.end() && current->second.prev.size() == 1 &&
                   prev_it->second.fec == 1) {
                if (prev_it != dbg.end()) {
                    for (auto it = current->second.seq.begin();
                         it != current->second.seq.end(); ++it) {
                        merged.emplace_back(std::move(*it));
                    }
                    cc_sum += current->second.cc;
                    dbg.erase(current->first);
                    simplified = true;
                } else {  // Prev is erased
                }
                current = prev_it;
                if (prev_it->second.prev.size() > 0) {
                    prev_it = dbg.find(current->second.prev[0]);
                } else {
                    break;
                }
            }
            if (current->first != id) {
                redirect_deleted[id] = current->first;
            }

            current->second.fec = first_fec;
            current->second.cc += cc_sum;
            for (auto it = merged.rbegin(); it != merged.rend(); ++it) {
                current->second.seq.emplace_back(std::move(*it));
            }

            return simplified;
        };
        for (int iter = 0; iter < max_iterations; ++iter) {
            // Remove simple paths

            for (auto &dn : lci) {
                auto current = dbg.find(dn.second.id);
                if (current == dbg.end()) {
                    dn.second.id = -1;
                    continue;
                }
                // Simple paths
                if ((current->second.fec != 1 &&
                     current->second.prev.size() == 1)) {
                    try_simplify_simple_path(dn.second.id, current);
                } else if (current->second.prev.size() > 1) {
                    for (auto &prev_id : current->second.prev) {
                        auto prev_it = dbg.find(prev_id);
                        if (prev_it != dbg.end() &&
                            prev_it->second.prev.size() == 1) {
                            try_simplify_simple_path(prev_id, prev_it);
                        }
                    }
                }

                // auto prev_it = dbg.find(node.prev[0]);
                // if (prev_it != dbg.end()) {
                //     auto &prev_node = prev_it->second;
                //     dn.second.pos = prev_node.merged.size();
                //     prev_node.merged.push_back(dn.second.id);
                //     prev_node.cc = node.cc;
                //     dbg.erase(dn.second.id);
                // } else {  // Prev is erased
                // }
                // merge_map.emplace(dn.second, node.prev[0], prev_node)
            }
            for (auto &node : dbg) {
                size_t idx = 0;
                int self_idx = -1;
                for (auto &p : node.second.prev) {
                    auto it = redirect_deleted.find(p);
                    if (it != redirect_deleted.end()) {
                        p = it->second;
                    } else if (p == node.first) {  // Self edge
                        node.second.fec--;
                        self_idx = idx;
                    }
                    ++idx;
                }
                if (self_idx != -1) {
                    node.second.prev.erase(node.second.prev.begin() + self_idx);
                }
            }

            // Tips
            for (auto &dn : lci) {
                auto current = dbg.find(dn.second.id);
                if (current == dbg.end()) {
                    dn.second.id = -1;
                    continue;
                }
                if (current->second.avg_core_count() < min_core_depth) {
                    // Forlorn
                    if (current->second.fec == 0 &&
                        current->second.prev.size() == 0) {
                        dbg.erase(current);
                    } else if (current->second.fec == 0 &&
                               current->second.prev.size() == 1) {
                        auto prev_it = dbg.find(current->second.prev[0]);
                        if (prev_it != dbg.end()) {
                            prev_it->second.fec--;
                        }
                        dbg.erase(current);
                    }
                    // else if (current->second.fec == 1 &&
                    //            current->second.prev.size() == 0) {
                    //     dbg.erase(current);
                    // }
                } else {
                    for (u64 &pid : current->second.prev) {
                        auto prev_it = dbg.find(pid);
                        if (prev_it != dbg.end()) {
                            if (prev_it->second.avg_core_count() <
                                    min_core_depth &&
                                prev_it->second.prev.size() == 0 &&
                                prev_it->second.fec == 1) {
                                dbg.erase(prev_it);
                                pid = std::numeric_limits<u64>::max();
                            }
                        }
                    }
                    auto &pvec = current->second.prev;
                    pvec.erase(std::remove_if(
                                   pvec.begin(), pvec.end(),
                                   [](u64 val) {
                                       return val ==
                                              std::numeric_limits<u64>::max();
                                   }),
                               pvec.end());
                    for (u64 &pid : current->second.prev) {
                        auto prev_it = dbg.find(pid);
                        if (prev_it != dbg.end()) {
                            if (prev_it->second.avg_core_count() <
                                    current->second.avg_core_count() &&
                                prev_it->second.seq.size() < max_tip_length &&
                                current->second.seq.size() >
                                    2 * max_tip_length &&
                                prev_it->second.prev.size() == 0 &&
                                prev_it->second.fec == 1) {
                                dbg.erase(prev_it);
                                pid = std::numeric_limits<u64>::max();
                            }
                        }
                    }
                    pvec.erase(std::remove_if(
                                   pvec.begin(), pvec.end(),
                                   [](u64 val) {
                                       return val ==
                                              std::numeric_limits<u64>::max();
                                   }),
                               pvec.end());
                }
                // We have non-zero forward edges, and more than one prev
                // edges
                // if (current->second.fec > 0 &&
                //     current->second.prev.size() > 0) {
                //     // Check if one of the prev edges is not terminating
                //     bool non_terminating_exists = false;
                //     size_t min_edge_len = INT_MAX;
                //     size_t max_edge_len = INT_MIN;
                //     for (auto &prev_id : current->second.prev) {
                //         auto prev_it = dbg.find(prev_id);
                //         if (prev_it == dbg.end()) {
                //             continue;
                //         }
                //         if (prev_it->second.prev.size() > 0) {
                //             non_terminating_exists = true;
                //         }
                //         min_edge_len =
                //             std::min(min_edge_len,
                //             prev_it->second.seq.size());
                //         max_edge_len =
                //             std::max(max_edge_len,
                //             prev_it->second.seq.size());
                //     }
                //     // Remove every tip if non_terminating_exists

                //     if (non_terminating_exists) {
                //         for (auto &prev_id : current->second.prev) {
                //             auto prev_it = dbg.find(prev_id);
                //             if (prev_it == dbg.end()) {
                //                 // delete from prev
                //                 current->second.prev.erase(
                //                     std::remove(current->second.prev.begin(),
                //                                 current->second.prev.end(),
                //                                 prev_id),
                //                     current->second.prev.end());
                //                 continue;
                //             }
                //             if (prev_it->second.prev.size() == 0 &&
                //                 prev_it->second.fec == 1 &&
                //                 prev_it->second.seq.size() < max_tip_length)
                //                 {
                //                 // Erase the node from the
                //                 current->second.prev
                //                 current->second.prev.erase(
                //                     std::remove(current->second.prev.begin(),
                //                                 current->second.prev.end(),
                //                                 prev_id),
                //                     current->second.prev.end());
                //                 dbg.erase(prev_it);
                //             }
                //         }
                //     } else if (max_edge_len - min_edge_len > max_tip_length)
                //     {
                //         for (auto &prev_id : current->second.prev) {
                //             auto prev_it = dbg.find(prev_id);
                //             if (prev_it == dbg.end()) {
                //                 continue;
                //             }
                //             if (prev_it->second.prev.size() == 0 &&
                //                 prev_it->second.fec == 1 &&
                //                 prev_it->second.seq.size() <
                //                     max_tip_length - max_tip_length) {
                //                 // Erase the node from the
                //                 current->second.prev

                //                 current->second.prev.erase(
                //                     std::remove(current->second.prev.begin(),
                //                                 current->second.prev.end(),
                //                                 prev_id),
                //                     current->second.prev.end());
                //                 dbg.erase(prev_it);
                //             }
                //         }
                //     }
                // }
            }

            // Buubles
            for (auto &dn : lci) {
                auto current = dbg.find(dn.second.id);
                if (current == dbg.end()) {
                    dn.second.id = -1;
                    continue;
                }

                if (current->second.prev.size() ==
                    2) {  // Simple case with 2 should expand later to many case
                    auto p1 = dbg.find(current->second.prev[0]);
                    auto p2 = dbg.find(current->second.prev[1]);
                    if (p1 != dbg.end() && p2 != dbg.end()) {
                        if (p1->second.prev.size() == 1 &&
                            p2->second.prev.size() == 1) {
                            if (p1->second.prev[0] == p2->second.prev[0]) {
                                // BUBBLE FOUND
                                if (p2->second.avg_core_count() /
                                        p1->second.avg_core_count() <
                                    bubble_depth_epsilon) {
                                    // Remove p2
                                    dbg.erase(p2);
                                    current->second.prev.pop_back();
                                    auto nn = dbg.find(p1->second.prev[0]);
                                    if (nn != dbg.end()) {
                                        nn->second.fec--;
                                    }
                                } else if (p1->second.avg_core_count() /
                                               p2->second.avg_core_count() <
                                           bubble_depth_epsilon) {
                                    // Remove p1
                                    dbg.erase(p1);
                                    current->second.prev[0] =
                                        current->second.prev[1];
                                    current->second.prev.pop_back();
                                    auto nn = dbg.find(p2->second.prev[0]);
                                    if (nn != dbg.end()) {
                                        nn->second.fec--;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Death by thousand cuts
            for (auto &dn : lci) {
                auto current = dbg.find(dn.second.id);
                if (current == dbg.end()) {
                    continue;
                }
                if (current->second.avg_core_count() < 2 &&
                    current->second.prev.size() == 1 &&
                    current->second.fec == 1) {
                    auto prev_it = dbg.find(current->second.prev[0]);
                    if (prev_it != dbg.end()) prev_it->second.fec--;
                    dbg.erase(current);
                }
            }
            for (auto &node : dbg) {
                auto &pvec = node.second.prev;
                for (auto it = pvec.begin(); it != pvec.end(); ++it) {
                    if (dbg.find(*it) == dbg.end()) {
                        *it = std::numeric_limits<u64>::max();
                    }
                }
                pvec.erase(
                    std::remove_if(pvec.begin(), pvec.end(),
                                   [](u64 val) {
                                       return val ==
                                              std::numeric_limits<u64>::max();
                                   }),
                    pvec.end());
            }
            redirect_deleted.clear();
        }

        erase_if(lci, [](const auto &node) {
            return node.second.id == std::numeric_limits<u64>::max();
        });
        // // Remove bubbles and tips
        // for (auto &node : dbg) {
        //     if (node.second.prev.size() == 2) {  // can be simplified.
        //         u64 p1 = node.second.prev[0];
        //         u64 p2 = node.second.prev[1];
        //         if (dbg[p1].cc == 1 && dbg[p2].cc == 1 &&
        //             dbg[p1].prev.size() == 1 && dbg[p2].prev.size() == 1 &&
        //             dbg[p1].prev[0] == dbg[p2].prev[0]) {
        //             ///     o
        //             ///   o< >o -> o->o->o
        //             ///     o
        //         } else if (node.second.prev.size() == 0) {
        //             // Remove tips
        //             /// ..o->o->o..       o->o->o
        //             ///       \     ->
        //             ///        o
        //         }
        //     }
        // }
    }
    u64 emplace_nobloomfilter(const lcmer<LCMER_L> &l, u64 &idx, int tid) {
        auto vv = lci.try_emplace(l, idx | tid);
        if (vv.second) {  // Means new item on the table
            idx += lcindex_increment;
            // threadcount.
        }
        return vv.first->second.id;
    }
    u64 emplace(const lcmer<LCMER_L> &l, u64 &idx, int tid) {
        auto new_hash = lci.hash(l);
        if (!bf.contains(new_hash)) {
            bf.insert(new_hash);
            return 0;
        }
        auto vv = lci.try_emplace_with_hash(new_hash, l, idx | tid);
        if (vv.second) {  // Means new item on the table
            idx += lcindex_increment;
            // threadcount.
        }
        return vv.first->second.id;
    }
};

struct assembler {
    string input;
    string output;
    decltype(asm_config::threads) threads;
    int lcp_level;
    double density;
    decltype(asm_config::bloom) bloom;
    uint64_t seed;
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType =  syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
    decltype(asm_config::syncmer) syncmer_params;
    SyncmerType syncmer_functor;
    assembler(const asm_config &cfg)
        : input{cfg.input},
          output{cfg.output},
          threads{cfg.threads},
          lcp_level{cfg.lcp_level},
          bloom{cfg.bloom},
          seed{cfg.seed},
          syncmer_params{cfg.syncmer},
          syncmer_functor{syncmer_params, static_cast<unsigned int>(seed)} {}


    static inline auto insert_lcmer_helper(const gtl::vector<uint64_t> &ll,
                                           const uint32_t &i, dbg_index &dbgi)
        -> bool {
        if (ll[i] == 0) {
            return false;
        }
        if (ll[i - 1] == 0) {
            if (ll[i + 1] == 0) {
                return false;
            } else {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](dbg_index::dbg_map_t::value_type &val) {
                        val.second.cc++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{1, 0, {}, {}});
                    });
            }
        } else {
            if (ll[i + 1] == 0) {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](dbg_index::dbg_map_t::value_type &val) {
                        val.second.emplace_prev(ll[i - 1]);
                        // val.second.cc++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{0, 0, {}, {ll[i - 1]}});
                    });
            } else {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](dbg_index::dbg_map_t::value_type &val) {
                        val.second.emplace_prev(ll[i - 1]);
                        val.second.cc++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{1, 0, {}, {ll[i - 1]}});
                    });
            }
        }
    }
    static void process_cores(const lcp_read &lr, dbg_index &dbgi,
                              int thread_id, uint64_t &lcmer_index) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - LCMER_L + 3);
        ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            ll.push_back(dbgi.emplace(l, lcmer_index, thread_id));
        }
        ll.push_back(0);

        for (size_t i = 1; i < ll.size() - 1; ++i) {
            if (insert_lcmer_helper(ll, i, dbgi)) {
                dbgi.dbg[ll[i]].seq.emplace_back(lr.data,
                                                 lr.pos[i + LCMER_L - 1].second,
                                                 lr.pos[i + LCMER_L].second);
            }
        }
    }
    static void process_cores_for_basespace(
        const lcp_read &lr, dbg_index &dbgi, [[maybe_unused]] int thread_id,
        [[maybe_unused]] uint64_t &lcmer_index) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        uint64_t val = dbgi.get(l);

        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            uint64_t from = lr.pos[i - 1].second;
            uint64_t to = lr.pos[i].second;

            dbgi.dbg.modify_if(val, [&](dbg_index::dbg_map_t::value_type &val) {
                // Lets copy the one-hot representation instead of decoding it
                // here.
                //  we need to copy from lr.pos[i-1].second to lr.pos[i].second.
                //  We need to pay attention to odd values (
                //   we need to keep some slack and make sure it is communicated
                //   to downstream
                //  )

                if (val.second.seq.size() == 0) {
                    if (to - from > 0) {
                        val.second.seq.emplace_back(lr.data, from, to);
                    }
                }
            });
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            val = dbgi.get(l);
        }
    }

    static std::pair<int, int> get_thread_dist(int threads) {
        if (threads == 0) {
            return {0, 0};
        }
        if (threads < 16) {
            int sam_threads = threads / 2;
            return {sam_threads, threads - sam_threads};
        }
        return {8, threads - 8};
    }
    using cqt = moodycamel::ConcurrentQueue<Chunk>;
    using dqt = moodycamel::BlockingConcurrentQueue<bam1_t *>;

    // static auto process_data_from_chunks(int thread_id, int lcp_level,
    //                                      string filename, cqt &chunk_queue,
    //                                      dbg_index &dbgi,
    //                                      std::atomic<bool> &done,
    //                                      int &processed_read_count,
    //                                      int &all_reads) {
    //     uint64_t lcmer_index = 0;
    // }
    static auto process_data(
        int thread_id, int lcp_level, dqt &data_queue, dqt &recycle_queue,
        dbg_index &dbgi, std::atomic<bool> &done, int &processed_read_count,
        int &all_reads, decltype(assembler::process_cores) Process_Func,
        SyncmerType syncmer_functor) {
        uint64_t lcmer_index = 0;
        bam1_t *data_p = NULL;
        lps str = {0, 0, 0, 0};
        charsi chs;
        lcp_read lr;
        while (true) {
            if (!data_queue.wait_dequeue_timed(data_p, WAIT_QUEUE_TIMEOUT)) {
                // if (!data_queue.try_dequeue(data_p)) {
                if (done.load()) {
                    break;
                }
                // std::this_thread::sleep_for(std::chrono::nanoseconds(100));
                continue;
            }
            all_reads++;
            lr.clear();
            parse_cores_from_bam_record(
                bam_get_seq(data_p), data_p->core.l_qseq, lcp_level, &str,
                syncmer_functor, chs, lr);
            if (lr.cores.size() < LCMER_L + 2) {
                // printf("Skipping read %lld, cores.size()=%zu\n", here,
                // lr.cores.size());
                goto pd_cleanup;
            }
            ++processed_read_count;
            Process_Func(lr, dbgi, thread_id, lcmer_index);
        pd_cleanup:
            if (!recycle_queue.try_dequeue(data_p)) {
                bam_destroy1(data_p);
                data_p = NULL;
            }
        }
        free_lps(&str);
    }

    auto generate_chunks(std::string input, cqt &chunk_queue) {
        BlockIterator iter(input);
        Block blk;
        // fetch chunks until no more blocks
        while (iter.next(blk)) {
            // first block in chunk
            int64_t first_vstart = blk.vstart();
            int64_t last_vend = blk.vend();

            size_t count = 1;
            while (count < BAMREAD_CHUNK_SIZE && iter.next(blk)) {
                last_vend = blk.vend();
                ++count;
            }
            while (!chunk_queue.try_enqueue({first_vstart, last_vend}));
        }
    }

    auto generate_data(const std::string &input, dqt &data_queue,
                       dqt &recycle_queue, int threads) {
        samFile *fp = sam_open(input.c_str(), "r");

        if (!fp) {
            std::cerr << "Error: Could not open " << input << std::endl;
            return;
        }
        hts_set_threads(fp, threads);
        // htsThreadPool pool = {NULL, 0};
        // if (threads > 1) {
        //     pool.pool = hts_tpool_init(threads);
        //     // hts_set_thread_pool(fp, &pool);
        //     hts_set_cache_size(fp, 1024 * 1024 * 1024);
        //     hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool);
        // }

        bam_hdr_t *hdr = sam_hdr_read(fp);
        bam1_t *rec = nullptr;
        while (true) {
            // Try to get a recycled record object
            if (rec == nullptr && !recycle_queue.try_dequeue(rec)) {
                // If the recycle queue is empty, allocate a new one as a
                // fallback
                rec = bam_init1();
            }

            if (sam_read1(fp, hdr, rec) >= 0) {
                while (!data_queue.try_enqueue(rec)) {
                    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
                }
                rec = nullptr;
            } else {
                // End of file or error. Destroy the unused record and exit.
                bam_destroy1(rec);
                break;
            }
        }
        bam_hdr_destroy(hdr);
        sam_close(fp);
        // if (pool.pool) {
        //     hts_tpool_destroy(pool.pool);
        // }
    }
    auto construct_dbg() {
        dqt data_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        std::atomic<bool> done{false};

        dbg_index dbgi{bloom.false_positive_rate, bloom.projected_count, seed};

        gtl::vector<std::thread> workers;

        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back(process_data, tid, lcp_level,
                                 std::ref(data_queue), std::ref(recycle_queue),
                                 std::ref(dbgi), std::ref(done),
                                 std::ref(valid_read_counts[tid - 1]),
                                 std::ref(all_read_counts[tid - 1]),
                                 process_cores, std::cref(syncmer_functor));
        }

        generate_data(input, std::ref(data_queue), std::ref(recycle_queue),
                      threads.io);

        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }

        int valid_sum = 0;
        int all_sum = 0;
        for (size_t i = 0; i < valid_read_counts.size(); ++i) {
            valid_sum += valid_read_counts[i];
            all_sum += all_read_counts[i];
        }

        fprintf(stderr, "%d processed\n%d total reads\n", valid_sum, all_sum);
        return dbgi;
    }

    auto convert2basespace(dbg_index &dbgi) -> int {
        dqt data_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads.work) * 1024, 1, 1};
        std::atomic<bool> done{false};

        gtl::vector<std::thread> workers;
        gtl::vector<int> valid_read_counts(threads.work, 0);
        gtl::vector<int> all_read_counts(threads.work, 0);

        for (int tid = 1; tid <= threads.work; ++tid) {
            workers.emplace_back(
                process_data, tid, lcp_level, std::ref(data_queue),
                std::ref(recycle_queue), std::ref(dbgi), std::ref(done),
                std::ref(valid_read_counts[tid - 1]),
                std::ref(all_read_counts[tid - 1]), process_cores_for_basespace,
                std::cref(syncmer_functor));
        }
        generate_data(input, std::ref(data_queue), std::ref(recycle_queue),
                      threads.io);

        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }

        int valid_sum = 0;
        int all_sum = 0;
        for (size_t i = 0; i < valid_read_counts.size(); ++i) {
            valid_sum += valid_read_counts[i];
            all_sum += all_read_counts[i];
        }

        fprintf(stderr, "%d processed\n%d total reads\n", valid_sum, all_sum);
        return 0;
    }
};

// TODO do not use extensions to identify input type.
input_type identify_input(const std::string &input) {
    std::string_view sv(input);
    if (sv.ends_with(".bam")) {
        return input_type::BAM;
    } else if (sv.ends_with(".sam")) {
        return input_type::SAM;
    } else if (sv.ends_with(".fastq") || sv.ends_with(".fq") ||
               sv.ends_with(".fasta") || sv.ends_with(".fa") ||
               sv.ends_with(".fna")) {
        return input_type::FASTX;
    } else if (sv.ends_with(".fastq.gz") || sv.ends_with(".fq.gz") ||
               sv.ends_with(".fasta.gz") || sv.ends_with(".fa.gz") ||
               sv.ends_with(".fna.gz")) {
        return input_type::FASTXGZ;
    } else if (sv.ends_with(".fastq.bz2") || sv.ends_with(".fq.bz2") ||
               sv.ends_with(".fasta.bz2") || sv.ends_with(".fa.bz2") ||
               sv.ends_with(".fna.bz2")) {
        return input_type::FASTXBZ2;
    }
    return input_type::UNKNOWN;
}

void print_raw_dbg(FILE *out, const dbg_index &dbgi) {
    for (const auto &[lval, edges] : dbgi.dbg) {
        fprintf(out, "S\tC%" PRIu64 "\t%s\tLN:i:%d\tRC:i:%d", lval, "*", 100,
                edges.cc);
    }

    for (const auto &[lval, edges] : dbgi.dbg) {
        for (const auto &p : edges.prev) {
            fprintf(out, "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", p, lval,
                    1);
        }
    }
}

void print_dbg_seqs(FILE *out, const dbg_index &dbgi) {
    charsi buf;
    for (const auto &[lval, edges] : dbgi.dbg) {
        if (edges.cc == 0 && edges.prev.size() == 0) {
            continue;
        }

        //        __builtin_debugtrap();
        if (edges.seq.size() == 0) {
            fprintf(out,
                    "S\tC%" PRIu64
                    "\t%s\tLN:i:%d\tRC:i:%u\tBE:i:%zu\tFE:i:%u\\n",
                    lval, "", 0, edges.cc, edges.prev.size(), edges.fec);
        } else {
            //            int i = onehotto2bitDNA2(s, buffer);
            fprintf(out, "S\tC%" PRIu64 "\t", lval);
            size_t sofar = 0;
            //            for (const bam_nibble_data &s : edges.seq) {
            // reversed
            for (auto it = edges.seq.begin(); it != edges.seq.end(); ++it) {
                const bam_nibble_data &s = *it;
                buf.require(s.len + 1);
                s.convert_to_ascii((char *)buf.data);
                fprintf(out, "%s", buf.data + (buf.data[0] == 'N' ? 1 : 0));
                sofar += s.len;
            }
            fprintf(out, "\tLN:i:%zu\tRC:i:%d\tCC:i:%u\tBE:i:%zu\tFE:i:%u\n",
                    sofar, (int)(edges.avg_core_count() * sofar), edges.cc,
                    edges.prev.size(), edges.fec);
        }

        for (uint64_t nl : edges.prev) {
            if (dbgi.dbg.contains(nl))
                fprintf(out, "L\tC%" PRIu64 "\t+\tC%" PRIu64 "\t+\t%dM\n", nl,
                        lval, 0);
        }
    }
}

void print_dbg_dot(const dbg_index &dbgi, int i) {
    string filename = string{"graph"} + std::to_string(i) + ".dot";
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
                    node.first, node.first, node.second.fec,
                    node.second.prev.size(), node.second.seq.size(),
                    node.second.cc);
            for (auto &prev : node.second.prev) {
                fprintf(dot_file, "\t%" PRIu64 " -> %" PRIu64 ";\n", prev,
                        node.first);
            }
        }
        fprintf(dot_file, "}\n");
        fclose(dot_file);
    }
}

int run_assembly(const asm_config &cfg) {
    cfg.print_cfg(stderr );
    INIT_TIMERS(HIFI_ASSEMBLY);
    auto inty = identify_input(cfg.input);

    if (inty == input_type::UNKNOWN) {
        fprintf(stderr, "Unknown input type %s\n", cfg.input.c_str());
        return -1;
    }

    assembler as(cfg);
    auto dbg = as.construct_dbg();
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG construction done.");

    dbg.count_forward_edges();
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Count forward edges done.");

    // as.convert2basespace(dbg);
    // TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Converting to basespace done.");

    dbg.simplify_graph(
        cfg.graph_cleaning.iterations, cfg.graph_cleaning.max_tip_length,
        cfg.graph_cleaning.min_bubble_length_difference, 3,
        cfg.graph_cleaning.min_bubble_depth_difference_percentage);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Simplify graph done.");

    FILE *out = fopen(cfg.output.c_str(), "w");
    print_dbg_seqs(out, dbg);
    // print_raw_dbg(out, dbg);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG printing done.");
    fclose(out);

    return 0;
}

void tests_mock_fill_seq(dbg_index &dbgi) {
    for (auto &node : dbgi.dbg) {
        node.second.seq.emplace_back();
    }
}

TEST_CASE("Graph simplification simple") {
    dbg_index dbgi{0.95, 10, 0};
    gtl::vector<u32> mock_cores;
    int MOCK_CORE_COUNT = 40;
    for (int i = 0; i < MOCK_CORE_COUNT; i++) {
        mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
    }
    lcmer<LCMER_L> l(mock_cores.begin());
    u64 idx = 0;

    gtl::vector<u64> ll(1, 0);
    ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
    ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    for (int i = LCMER_L + 1; i < MOCK_CORE_COUNT - 1; ++i) {
        uint32_t in_core = mock_cores[i];
        l.update_lcmer(in_core);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    }
    ll.push_back(0);

    for (size_t i = 1; i <= ll.size() - 1; i++) {
        assembler::insert_lcmer_helper(ll, i, dbgi);
    }

    SUBCASE("Check ll vec") {
        for (size_t i = 1; i < ll.size() - 1; ++i) {
            CHECK(ll[i] == (i - 1) * 128);
        }
    }
    SUBCASE("Pre simplify dbg") {
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
            CHECK(node.second.cc < 2);
            CHECK(node.second.prev.size() < 2);
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.lci.size() == MOCK_CORE_COUNT - 2 - LCMER_L + 1);
    }

    tests_mock_fill_seq(dbgi);
    for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
        n.second.cc *= 10;
    }
    dbgi.count_forward_edges();
    dbgi.simplify_graph(1);
    SUBCASE("Post simplify dbg") {
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.dbg.size() == 1);
        CHECK(dbgi.dbg.begin()->second.seq.size() == 30);
    }

    // // Write Dot file
    // FILE *dot_file = fopen("graph.dot", "w");

    // if (dot_file) {
    //     // print_raw_dbg(dot_file, dbgi);
    //     fprintf(dot_file, "digraph G {\n");
    //     for (auto &node : dbgi.dbg) {
    //         fprintf(dot_file, "\t%" PRIu64 " [label=\"%" PRIu64 "\"];\n",
    //         node.first,
    //                 node.first);
    //         for (auto &prev : node.second.prev) {
    //             fprintf(dot_file, "\t%" PRIu64 " -> %" PRIu64 ";\n", prev,
    //             node.first);
    //         }
    //     }
    //     fprintf(dot_file, "}\n");
    //     fclose(dot_file);
    // }
}

TEST_CASE("Graph simplification Bubble") {
    dbg_index dbgi{0.95, 10, 0};
    gtl::vector<u32> mock_cores;
    int MOCK_CORE_COUNT = 40;
    for (int i = 0; i < MOCK_CORE_COUNT; i++) {
        mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
    }
    lcmer<LCMER_L> l(mock_cores.begin());
    u64 idx = 0;

    gtl::vector<u64> ll(1, 0);
    ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
    ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
        uint32_t in_core = mock_cores[i];
        l.update_lcmer(in_core);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    }
    ll.push_back(0);
    for (size_t i = 1; i <= ll.size() - 1; i++) {
        assembler::insert_lcmer_helper(ll, i, dbgi);
    }

    ll.clear();
    mock_cores[MOCK_CORE_COUNT / 2] = 123123123;
    l = lcmer<LCMER_L>(mock_cores.begin());

    ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
        uint32_t in_core = mock_cores[i];
        l.update_lcmer(in_core);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    }
    ll.push_back(0);
    for (size_t i = 1; i <= ll.size() - 1; i++) {
        assembler::insert_lcmer_helper(ll, i, dbgi);
    }

    dbgi.count_forward_edges();
    for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
        n.second.cc *= 10;
    }
    tests_mock_fill_seq(dbgi);
    print_dbg_dot(dbgi, 0);

    dbgi.simplify_graph(1);
    print_dbg_dot(dbgi, 1);

    SUBCASE("Post simplify dbg") {
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.dbg.size() == 3);
    }
    dbgi.simplify_graph(1);

    SUBCASE("Post simplify dbg iter 2") {
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.dbg.size() == 1);
    }

    print_dbg_dot(dbgi, 2);
}

TEST_CASE("Graph simplification Tips (Outgoing)") {
    dbg_index dbgi{0.95, 10, 0};
    gtl::vector<u32> mock_cores;
    int MOCK_CORE_COUNT = 30;
    for (int i = 0; i < MOCK_CORE_COUNT; i++) {
        mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
    }
    lcmer<LCMER_L> l(mock_cores.begin());
    u64 idx = 0;

    gtl::vector<u64> ll(1, 0);
    ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
    ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
        uint32_t in_core = mock_cores[i];
        l.update_lcmer(in_core);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    }
    ll.push_back(0);
    for (size_t i = 1; i <= ll.size() - 1; i++) {
        assembler::insert_lcmer_helper(ll, i, dbgi);
    }
    for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
        n.second.cc *= 10;
    }
    ll.clear();
    mock_cores[12] = 123123123;
    l = lcmer<LCMER_L>(mock_cores.begin());

    ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
        uint32_t in_core = mock_cores[i];
        l.update_lcmer(in_core);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    }
    ll.push_back(0);
    for (size_t i = 1; i <= ll.size() / 2 - 1; i++) {
        assembler::insert_lcmer_helper(ll, i, dbgi);
    }

    dbgi.count_forward_edges();

    tests_mock_fill_seq(dbgi);
    print_dbg_dot(dbgi, 3);

    dbgi.simplify_graph(1);
    SUBCASE("Post simplify dbg") {
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.dbg.size() == 2);
    }
    print_dbg_dot(dbgi, 4);

    dbgi.simplify_graph(1);
    SUBCASE("Post simplify dbg") {
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.dbg.size() == 1);
    }
    print_dbg_dot(dbgi, 5);
}
TEST_CASE("Graph simplification Tips (Incoming)") {
    dbg_index dbgi{0.95, 10, 0};
    gtl::vector<u32> mock_cores;
    lcmer<LCMER_L> l;
    gtl::vector<u64> ll;

    u64 idx = 0;
    int MOCK_CORE_COUNT = 30;
    for (int i = 0; i < MOCK_CORE_COUNT; i++) {
        mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
    }
    for (int DEPTH = 0; DEPTH < 10; ++DEPTH) {
        l = lcmer<LCMER_L>(mock_cores.begin());

        ll.push_back(0);

        ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
        for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
            uint32_t in_core = mock_cores[i];
            l.update_lcmer(in_core);
            ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
        }
        ll.push_back(0);
        for (size_t i = 1; i <= ll.size() - 1; i++) {
            assembler::insert_lcmer_helper(ll, i, dbgi);
        }

        ll.clear();
    }
    ll.push_back(0);
    mock_cores[8] = 123123123;
    l = lcmer<LCMER_L>(mock_cores.begin() + 7);

    ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    for (int i = LCMER_L + 7; i < MOCK_CORE_COUNT - 1; ++i) {
        uint32_t in_core = mock_cores[i];
        l.update_lcmer(in_core);
        ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
    }
    ll.push_back(0);
    for (size_t i = 1; i <= ll.size() / 2 - 1; i++) {
        assembler::insert_lcmer_helper(ll, i, dbgi);
    }

    dbgi.count_forward_edges();

    tests_mock_fill_seq(dbgi);
    print_dbg_dot(dbgi, 6);

    SUBCASE("Post simplify dbg") {
        dbgi.simplify_graph(2);
        for (auto &node : dbgi.dbg) {
            INFO(node.first, "\t.\t.\t", node.second.seq.size());
        }
        CAPTURE(dbgi.dbg.size());
        CHECK(dbgi.dbg.size() == 2);
    }

    print_dbg_dot(dbgi, 7);
}
