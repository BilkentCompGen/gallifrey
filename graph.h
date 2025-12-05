#ifndef _GRAPH_H_
#define _GRAPH_H_

#define PREV_WITH_VEC
#define HT_SHARD_BC 8

#include <cstdint>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <mutex>

#include "bloomfilter.h"
#include "config.h"
#include "lcmer.h"
#include "nibble2base.hpp"
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
using u32 = uint32_t;

// using SequenceContainerType = bam_nibble_data;
using SequenceContainerType = CrumbledSubSeq;

struct dbg_index {
    static constexpr u64 lcindex_increment =
        nextpowerof2((u64)MAXIMUM_WORK_THREADS);
    using mutex_type = std::mutex;
    struct bidbg_val {
        unsigned cc;
        unsigned fec;
        gtl::vector<SequenceContainerType> seq;
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
            gtl::vector<SequenceContainerType> merged;
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

                if (current->second.prev.size() >=
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
                                    current->second.prev.erase(std::next(
                                        current->second.prev.begin()));
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
                                    current->second.prev.erase(std::next(
                                        current->second.prev.begin()));

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
            // Simple Loops

            for (auto &dn : lci) {
                auto current = dbg.find(dn.second.id);
                if (current == dbg.end()) {
                    continue;
                }
                if (current->second.fec == 1 &&
                    current->second.prev.size() == 1) {
                    auto prev_it = dbg.find(current->second.prev[0]);
                    if (prev_it != dbg.end()) {
                        auto &prepre = prev_it->second.prev;
                        for (auto pp = prepre.begin(); pp != prepre.end();
                             ++pp) {
                            if (*pp == dn.second.id) {  // Loop here
                                *pp = std::numeric_limits<u64>::max();
                                if (current->second.cc >=
                                    2 * prev_it->second
                                            .cc) {  // TODO: this is not right.
                                                    // Need to duplicate the
                                                    // sequence
                                    auto &tomerge = current->second.seq;
                                    for (auto it = tomerge.begin();
                                         it != tomerge.end(); ++it) {
                                        prev_it->second.seq.emplace_back(
                                            std::move(*it));
                                    }
                                }

                                prev_it->second.fec--;
                                dbg.erase(current);
                            }
                        }
                    }
                }
            }
            // Fix edges
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

inline void print_raw_dbg(FILE *out, const dbg_index &dbgi) {
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

inline void print_dbg_seqs(FILE *out, const dbg_index &dbgi) {
    charsi buf;
    for (const auto &[lval, edges] : dbgi.dbg) {
        if (edges.cc == 0 && edges.prev.size() == 0) {
            continue;
        }

        //        __builtin_debugtrap();
        if (edges.seq.size() == 0) {
            fprintf(out,
                    "S\tC%" PRIu64
                    "\t%s\tLN:i:%d\tRC:i:%u\tBE:i:%zu\tFE:i:%u\n",
                    lval, "", 0, edges.cc, edges.prev.size(), edges.fec);
        } else {
            //            int i = onehotto2bitDNA2(s, buffer);
            fprintf(out, "S\tC%" PRIu64 "\t", lval);
            size_t sofar = 0;
            //            for (const SequenceContainerType &s : edges.seq) {
            // reversed
            for (auto it = edges.seq.begin(); it != edges.seq.end(); ++it) {
                const SequenceContainerType &s = *it;
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

#endif  // _GRAPH_H_
