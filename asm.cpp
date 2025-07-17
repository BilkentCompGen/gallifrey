#include "asm.h"

#include <_stdio.h>
#include <concurrentqueue.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>
#include <htslib/hts.h>


#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <tuple>

#include "bgzf.h"
#include "config.h"
#include "ext/bloom_filter.hpp"
#include "gtl/phmap_fwd_decl.hpp"
#include "lcmer.h"

// ......... 8 overlaps 16 bytes. since we overlap we use each 8 times. which
// means each overlap takes 16 bytes instead.

using mutex_type = std::mutex;
int hts_useek(htsFile *fp, off_t uoffset, int where)
{
    if (fp->is_bgzf)
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
        return (hseek(fp->fp.hfile, uoffset, SEEK_SET) >= 0)? 0 : -1;
}
off_t hts_utell(htsFile *fp)
{
    if (fp->is_bgzf)
        return bgzf_utell(fp->fp.bgzf);
    else
        return htell(fp->fp.hfile);
}
struct bidbg_val {
    gtl::flat_hash_set<uint64_t> prev;
    gtl::flat_hash_set<uint64_t> next;
    int cc;
    gtl::vector<std::pair<uint64_t, int>> merged;
};
using bidbg = gtl::parallel_flat_hash_map<
    uint64_t, bidbg_val, std::hash<uint64_t>,
    gtl::priv::hash_default_eq<uint64_t>,
    std::allocator<std::pair<const uint64_t, bidbg_val>>, 12, mutex_type>;
using lcmer_index_map_t = gtl::parallel_node_hash_map<
    lcmer<LCMER_L>, uint64_t, std::hash<lcmer<LCMER_L>>,
    gtl::priv::hash_default_eq<lcmer<LCMER_L>>,
    std::allocator<std::pair<const lcmer<LCMER_L>, uint64_t>>, 12, mutex_type>;

using reverse_lcmer_idx_t = gtl::parallel_node_hash_map<
    uint64_t, const lcmer<LCMER_L>, gtl::priv::hash_default_hash<uint64_t>,
    gtl::priv::hash_default_eq<uint64_t>,
    std::allocator<std::pair<const uint64_t, const lcmer<LCMER_L> *>>, 4,
    mutex_type>;
struct lcmer_index_t {
    lcmer_index_map_t map;
    reverse_lcmer_idx_t rev;
    bloom_filter bf;
    lcmer_index_t(const bloom_parameters &bparam) : bf{bparam} {}

    uint64_t emplace(const lcmer<LCMER_L> &l, uint64_t &idx, int tid) {
        if (!bf.contains(l)) {
            bf.insert(l);
            return 0;
        }
        auto vv = map.try_emplace(l, idx | tid);
        if (vv.second) {
            // rev.try_emplace(vv.first->second, vv.first->first);
            idx += 128;
        }
        // map.lazy_emplace_l(
        //     l, [&](lcmer_index_map_t::value_type &val) { ret = val.second; },
        //     [&](const lcmer_index_map_t::constructor &ctr) {
        //         ctr(l, idx | tid);
        //         ret = idx | tid;
        //         idx += 128;
        //     });
        // rev.try_emplace(ret, &map.find(l)->first);
        return vv.first->second;
    }
};

using bgzf_block = int64_t;
struct coreaddress_in_bgzf {
    bgzf_block seek;
    uint32_t idx;
    uint32_t from;
    uint32_t to;
};

using overlap_idx_t = gtl::parallel_flat_hash_map<
    uint64_t, int32_t, gtl::priv::hash_default_hash<uint64_t>,
    gtl::priv::hash_default_eq<uint64_t>,
    std::allocator<std::pair<const uint64_t, int32_t>>, 6, mutex_type>;

using core2bgzf_index_t = gtl::parallel_flat_hash_map<
    uint32_t, coreaddress_in_bgzf, gtl::priv::hash_default_hash<uint32_t>,
    gtl::priv::hash_default_eq<uint32_t>,
    std::allocator<std::pair<const uint32_t, coreaddress_in_bgzf>>, 6,
    mutex_type>;

struct BGZFCoreIndex {
    std::string filename;
    core2bgzf_index_t idx;
    auto emplace(uint32_t core, const coreaddress_in_bgzf &a) {
        return idx.lazy_emplace_l(
            core,
            [&a](core2bgzf_index_t::value_type &val) {
                if (val.second.idx > a.idx) {
                    val.second = a;
                }
            },
            [&](const core2bgzf_index_t::constructor &ctor) { ctor(core, a); });
    }
};

struct edge {
    uint32_t to;
    int32_t next = -1;
    uint64_t id;
    int count;
};

bool validate_sam_read(int ret, const BAMReader &reader, const Chunk &chunk) {
    if (ret == -1) {
        fprintf(stderr, "End of file %s at %lld %lld\n",
                reader.filename.c_str(), chunk.vstart, chunk.vend);
        return false;
    } else if (ret < -1) {
        fprintf(stderr,
                "sam_read1 failed on %s %lld %lld possibly due to a "
                "read "
                "being longer than block size in the previous block\n",
                reader.filename.c_str(), chunk.vstart, chunk.vend);
        return false;
    }
    return true;
}

/**
 * @brief Generates chunks of SAM/BAM blocks and enqueues them into a concurrent
 * queue.
 *
 * This function instantiates a BAMReader for the given input filename, iterates
 * over its blocks using a BlockIterator, and groups consecutive blocks up to a
 * specified chunk size (BAMREAD_CHUNK_SIZE). It then attempts to enqueue each
 * generated chunk into a moodycamel::ConcurrentQueue of Chunk objects.
 *
 * @param input The path to the input SAM/BAM file.
 * @param chunk_queue A reference to a moodycamel::ConcurrentQueue where chunks
 * will be enqueued.
 */
auto generate_chunks(std::string input,
                     moodycamel::ConcurrentQueue<Chunk> &chunk_queue) {
    BAMReader reader{input};
    BlockIterator iter(reader.fp);
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
/**
 * @brief Processes cores by generating LCmers and calculates overlaps stores
 * core locations on the bam file.
 *
 * This function takes a reference to an lcp_read object representing read cores
 * and their positions, along with various indices (LCMer index, overlap
 * indices, and debugging structure) and a BGZF core index reference. It
 * constructs an LCmer from the provided read cores (skipping the first core),
 * registers it in the LCmer index, and iterates through the remaining cores to
 * update overlap information as well as maintain core location information.
 *
 * The function uses a two-phase approach:
 *   - In the first phase (from index LCMER_L+1 onward), it registers cores and
 * updates the overlap information based on previous core positions.
 *   - In the second phase, it further updates LCmer indices for subsequent
 * cores and records bidirectional connections between consecutive LCmers in the
 * dbg.
 *
 * Each core is registered with its corresponding BGZF block information
 * (position and index), which aids in later retrieval and analysis.
 *
 * @param lr The read object containing core values and their positions.
 * @param lci Reference to the LCmer index structure where LCmers are emplaced.
 * @param thread_id Identifier for the current processing thread.
 * @param lcmer_index Reference to a 64-bit index used for LCmer registration.
 * @param ovlps A reference to the overlap indices container for storing
 * computed overlaps.
 * @param dbg Reference to the dbg
 * @param core_index Reference to the BGZF core index structure used for
 * registering cores.
 * @param bgzf_here The current BGZF block information (position or identifier)
 * for the core.
 * @param read_index Identifier for the current read in processing.
 */
auto process_cores(const lcp_read &lr, lcmer_index_t &lci, int thread_id,
                   uint64_t &lcmer_index, overlap_idx_t &ovlps, bidbg &dbg,
                   BGZFCoreIndex &core_index, int64_t bgzf_here,
                   unsigned read_index) {
    lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
    uint64_t two_cores = lr.cores[0];
    core_index.emplace(
        lr.cores[1], coreaddress_in_bgzf{bgzf_here, read_index, lr.pos[1].first,
                                         lr.pos[1].second});

    for (int i = 1; i < LCMER_L + 1; ++i) {
        two_cores = (two_cores << 32) | lr.cores[i];
        ovlps.try_emplace(two_cores, lr.pos[i - 1].second - lr.pos[i].first);
        core_index.emplace(lr.cores[i], coreaddress_in_bgzf{
                                            bgzf_here, read_index,
                                            lr.pos[i].first, lr.pos[i].second});
    }
    gtl::vector<uint64_t> ll(1, 0);
    ll.reserve(lr.cores.size() - LCMER_L + 3);
    ll.push_back(lci.emplace(l, lcmer_index, thread_id));
    for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
        two_cores = (two_cores << 32) | lr.cores[i];
        ovlps.try_emplace(two_cores, lr.pos[i - 1].second - lr.pos[i].first);
        uint32_t in_core = lr.cores[i];
        uint32_t out_core = l.update_lcmer(in_core);
        ll.push_back(lci.emplace(l, lcmer_index, thread_id));
        core_index.emplace(lr.cores[i], coreaddress_in_bgzf{
                                            bgzf_here, read_index,
                                            lr.pos[i].first, lr.pos[i].second});
    }
    ll.push_back(0);

    for (size_t i = 1; i < ll.size() - 1; ++i) {
        if (ll[i] == 0) {
            continue;
        }
        if (ll[i - 1] == 0) {
            if (ll[i + 1] == 0) {
                continue;
            } else {
                dbg.lazy_emplace_l(
                    ll[i],
                    [&](bidbg::value_type &val) {
                        val.second.next.emplace(ll[i + 1]);
                        val.second.cc++;
                    },
                    [&](const bidbg::constructor &ctr) {
                        ctr(ll[i], bidbg_val{{}, {ll[i + 1]}, 1});
                    });
            }
        } else {
            if (ll[i + 1] == 0) {
                dbg.lazy_emplace_l(
                    ll[i],
                    [&](bidbg::value_type &val) {
                        val.second.prev.emplace(ll[i - 1]);
                        val.second.cc++;
                    },
                    [&](const bidbg::constructor &ctr) {
                        ctr(ll[i], bidbg_val{{ll[i - 1]}, {}, 1});
                    });
            } else {
                dbg.lazy_emplace_l(
                    ll[i],
                    [&](bidbg::value_type &val) {
                        val.second.prev.emplace(ll[i - 1]);
                        val.second.next.emplace(ll[i + 1]);
                        val.second.cc++;
                    },
                    [&](const bidbg::constructor &ctr) {
                        ctr(ll[i], bidbg_val{{ll[i - 1]}, {ll[i + 1]}, 1});
                    });
            }
        }
    }
}

auto process_chunks(const std::string &input, int thread_id, int lcp_level,
                    moodycamel::ConcurrentQueue<Chunk> &chunk_queue, bidbg &dbg,
                    overlap_idx_t &ovlps, BGZFCoreIndex &core_loc_idx,
                    lcmer_index_t &lci, std::atomic<bool> &done,
                    int &processed_read_count, int &all_reads) {
    uint64_t lcmer_index = 0;
    Chunk chunk;
    BAMReader reader{input};
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(),
                                                         bam_destroy1);
    sam_hdr_t *hdr = NULL;

    while (chunk_queue.size_approx() > 0 || !done.load()) {
        if (!chunk_queue.try_dequeue(chunk)) {
            continue;
        }
        if (chunk.vstart > 0) {
            if (bgzf_seek(reader.fp->fp.bgzf, chunk.vstart, SEEK_SET) < 0) {
                std::cerr << "bgzf_seek failed on chunk vstart=" << chunk.vstart
                          << "\n";
                std::exit(EXIT_FAILURE);
            }
        }
        int src = 0;
        while (true) {
            int64_t here = bgzf_tell(reader.fp->fp.bgzf);
            if (here >= chunk.vend) {
                break;
            }
            auto telled = hts_utell(reader.fp);
            int ret = sam_read1(reader.fp, hdr, rec.get());
            ++src;
            if (!validate_sam_read(ret, reader, chunk)) {
                break;
            }
            lcp_read lr = parse_cores_from_bam_record(
                bam_get_seq(rec.get()), rec->core.l_qseq, lcp_level);
            if (lr.cores.size() < LCMER_L + 2) {
                // printf("Skipping read %lld, cores.size()=%zu\n", here,
                // lr.cores.size());
                continue;
            }
            ++processed_read_count;
            process_cores(lr, lci, thread_id, lcmer_index, ovlps, dbg,
                          core_loc_idx, telled, src);
        }
        all_reads += src;
        // printf("%zu\n", );
    }
}

auto construct_dbg_v4(std::string input, int threads, int lcp_level,
                      const bloom_parameters &bparam) {
    moodycamel::ConcurrentQueue<Chunk> chunk_queue{100000, 1, 1};
    std::atomic<bool> done{false};

    bidbg dbg;
    lcmer_index_t lci{bparam};
    BGZFCoreIndex core_index{.filename = input, .idx = {}};
    overlap_idx_t ovlp_idx;

    gtl::vector<std::thread> workers;
    gtl::vector<int> rc(threads - 1, 0);
    gtl::vector<int> all_read_count(threads - 1, 0);
    for (int i = 0; i < threads - 1; ++i) {
        workers.emplace_back(process_chunks, input, i + 1, lcp_level,
                             std::ref(chunk_queue), std::ref(dbg),
                             std::ref(ovlp_idx), std::ref(core_index),
                             std::ref(lci), std::ref(done), std::ref(rc[i]),
                             std::ref(all_read_count[i]));
    }
    generate_chunks(input, chunk_queue);
    done.store(true);
    for (auto &worker : workers) {
        worker.join();
    }
    int all_reads = 0;
    int used_reads = 0;
    for (int i = 0; i < threads - 1; ++i) {
        used_reads += rc[i];
        all_reads += all_read_count[i];
    }
    fprintf(stderr, "Processed %d reads: Used %d for the dbg and skipped %d\n",
            all_reads, used_reads, all_reads - used_reads);
    return std::make_tuple(lci, core_index, ovlp_idx, dbg);
}

struct sv {
    char *from;
    char *to;
    sv() : from(0), to(0) {}
    sv(char *from, char *to) : from(from), to(to) {}

    void expand_until(char c) {
        while (*to != 0 && *to != c) {
            to++;
        }
    }
    void skip_next(char c) {
        expand_until(c);
        if (*to) ++to;
        from = to;
    }
    template <int N, char DEL>
    void get() {
        for (int i = 0; i < N; ++i) {
            skip_next(DEL);
        }
        expand_until(DEL);
    }
    operator std::string() { return {from, to}; }
};

struct unitig_ar_graph_t {
    gtl::vector<edge> edges;
    gtl::vector<uint32_t> le;
    gtl::vector<gtl::vector<uint32_t>> cores;
    void add_node() {
        le.push_back(-1);
        cores.emplace_back();
    }
    void expand_run(uint32_t core) { cores.back().push_back(core); }
    void add_edge(uint32_t from, uint32_t to) {
        edges.emplace_back(to, le[from], from, 1);
        le[from] = edges.size() - 1;
    }
};

// unitig_ar_graph_t unitigs_gfatools(const reverse_lcmer_idx_t &idx) {
//     FILE *gfatools_pipe = popen("gfatools asm -u temp.gfa", "r");
//     char *line = NULL;
//     size_t lwc = 0;
//     unitig_ar_graph_t g;

//     while (getline(&line, &lwc, gfatools_pipe) != -1) {
//         sv s{line + 2, line + 2};
//         switch (line[0]) {
//             case 'S': {
//                 g.add_node();
//                 break;
//             }
//             case 'A': {
//                 s.get<3, '\t'>();
//                 int lc_idx;
//                 std::from_chars(s.from, s.to, lc_idx);
//                 if (g.cores.back().empty()) {
//                     for (uint32_t core : idx.at(lc_idx)->data) {
//                         g.expand_run(core);
//                     }
//                 } else {
//                     g.expand_run(idx.at(lc_idx)->data[0]);
//                 }
//                 break;
//             }
//             case 'L': {
//                 s.expand_until('\t');
//                 int u1, u2;
//                 auto result = std::from_chars(s.from + 3, s.to - 1, u1);
//                 s.get<2, '\t'>();
//                 result = std::from_chars(s.from + 3, s.to - 1, u2);
//                 g.add_edge(u1 - 1, u2 - 1);
//                 break;
//             }
//             default:
//                 break;
//         }
//     }
//     pclose(gfatools_pipe);
//     return g;
// }

std::tuple<uint8_t *, int, int> find_core_seq(BAMReader &reader,
                                              const BGZFCoreIndex &core_idx,
                                              uint32_t core, bam1_t *rec) {
    auto seq_location = core_idx.idx.find(core);
    if (seq_location == core_idx.idx.end()) {
        return {};
    }
    // int64_t here = bgzf_tell(reader.fp->fp.bgzf);
    sam_hdr_t *hdr = NULL;
    if (seq_location->second.seek > 0) {
        if (bgzf_seek(reader.fp->fp.bgzf, seq_location->second.seek, SEEK_SET) <
            0) {
            std::cerr << "bgzf_seek failed on chunk vstart="
                      << seq_location->second.seek << "\n";
            std::exit(EXIT_FAILURE);
        }
    } else {
        // reader.fp = bam_
    }

    // if(reader.seek(seq_location->second.seek) == -1){
    //     fprintf(stderr, "Seek error\n");
    // }
    int ret;
    unsigned i;
    // for (i = 1; i < seq_location->second.idx; ++i) {
    ret = bam_read1(reader.fp->fp.bgzf, rec);
    if (ret == -1) {
        fprintf(stderr, "End of stream in find_core_seq %d %u %lld\n", ret, i,
                seq_location->second.seek);
        return {};
    }
    // }
    if (ret < -1) {
        fprintf(stderr, "Error in find_core_seq %d %u\n", ret, core);
        throw std::string{"BAD"};
    }
    if (hdr) {
        sam_hdr_destroy(hdr);
    }
    return {bam_get_seq(rec), seq_location->second.from,
            seq_location->second.to};
}

/*    bgzf_block seek;
    uint32_t idx;
    uint32_t from;
    uint32_t to;
*/

void to_basespace(FILE *out, const unitig_ar_graph_t &g,
                  const BGZFCoreIndex &core_idx,
                  const overlap_idx_t &overlaps) {
    BAMReader reader(core_idx.filename);
    std::cout << "Reading from " << core_idx.filename << "\n";
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(),
                                                         bam_destroy1);
    static uint8_t alph[16] = {0,  65, 67, 0, 71, 0, 0, 0,
                               84, 0,  0,  0, 0,  0, 0, 78};

    struct seek_job_t {
        size_t i;
        int64_t seek;
        size_t idx;
        int32_t from;
        int32_t to;
    };
    for (size_t i = 0; i < g.cores.size(); ++i) {
        fprintf(out, "S\tU%zul\t", i);
        int ovlp = 0;
        uint64_t two_cores = g.cores[i].front();
        int len = 0;
        gtl::vector<char *> core_seqs(g.cores.size(), NULL);
        gtl::vector<seek_job_t> seekj;
        for (int j = 0; j < g.cores[i].size(); ++j) {
            auto [_seek, _idx, _f, _t] = core_idx.idx.at(g.cores[i][j]);
            seekj.emplace_back(j, _seek, _idx, _f, _t);
        }
        std::sort(seekj.begin(), seekj.end(),
                  [](const seek_job_t &a, const seek_job_t &b) {
                      if (a.seek == b.seek) {
                          return a.idx < b.idx;
                      }
                      return a.seek < b.seek;
                  });
        for (auto it = g.cores[i].begin() + 1; it != g.cores[i].end(); ++it) {
            auto [s, fro, to] = find_core_seq(reader, core_idx, *it, rec.get());
            // fprintf(stderr, "%d\t%d\n", fro, to);
            for (int j = ovlp + fro; j < to; ++j) {
                fputc(alph[bam_seqi(s, j)], out);
                ++len;
            }
            two_cores = (two_cores << 32) | *it;
            auto oit = overlaps.find(two_cores);
            if (oit != overlaps.end())
                ovlp = oit->second;
            else
                ovlp = 0;
        }
        fprintf(out, "\tLN:i:%d\n", len);
        auto it = g.le[i];
        while (it != -1) {
            fprintf(out, "L\tU%zul\t+\tU%ul\t+\t0M\n", i, g.edges[it].to);
            it = g.edges[it].next;
        }
    }
}

void print_dbg_v4(FILE *out, bidbg &dbg, lcmer_index_t &lci) {
    for (const auto &[lval, edges] : dbg) {
        int sum_cc = edges.cc;
        for(const auto &p : edges.merged){
            sum_cc += p.second;
        }
        fprintf(out, "S\tC%llu\t*\tLN:i:%zu\tRC:i:%d\n", lval,
                LCMER_L * (edges.merged.size() + 1) * 100, sum_cc);
        for (uint64_t nl : edges.prev) {
            fprintf(out, "L\tC%llu\t+\tC%llu\t+\t%dM\n", nl, lval, 100 * (LCMER_L - 1));
        }
        // for (uint64_t nl : edges.next) {
        //     fprintf(out, "L\tC%llu\t-\tC%llu\t-\t0M\n", lval, nl);
        // }
    }
}

void simplify_dbg(bidbg &dbg) {
    for (auto &pair : dbg) {
        if (pair.second.cc == -176) {
            continue;
        }
        while (pair.second.next.size() == 1) {
            auto np = dbg.find(*pair.second.next.begin());
            if (np->second.prev.size() != 1) {
                if (np->second.prev.size() > 1 &&
                    pair.second.merged.size() > 0) {
                    // for(const auto &p : np->second.prev){
                    //     std::cout << p << "\t" <<
                    //     pair.second.merged.back().first << "\n";
                    // }
                    np->second.prev.erase(pair.second.merged.back().first);
                    np->second.prev.emplace(pair.first);
                }
                break;
            }
            pair.second.merged.emplace_back(np->first, np->second.cc);
            pair.second.merged.insert(pair.second.merged.end(),
                                      np->second.merged.begin(),
                                      np->second.merged.end());
            pair.second.next = np->second.next;

            np->second.cc = -176;
            if (np->second.merged.size() > 0) {
                break;
            }
        }
    }
    for (auto it = dbg.begin(), end = dbg.end(); it != end;) {
        if (it->second.cc == -176) {
            dbg._erase(it++);
        } else {
            ++it;
        }
    }
}

struct bgzf_sort_helper {
    int i;
    coreaddress_in_bgzf address;
};

void fill_seqs(gtl::vector<bgzf_sort_helper> &to_sort, int l, int r,
               const BGZFCoreIndex &bgzf, gtl::vector<char *> &ret,
               BAMReader &reader) {
    auto seek_pos = to_sort[l].address.seek;
    auto seek_end = to_sort[r].address.seek;
    if (seek_pos > 0) {
        if (bgzf_seek(reader.fp->fp.bgzf, seek_pos, SEEK_SET) < 0) {
            std::cerr << "bgzf_seek failed on chunk vstart=" << seek_pos
                      << "\n";
            std::exit(EXIT_FAILURE);
        }
    }
    sam_hdr_t *hdr = NULL;
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(),
                                                         bam_destroy1);
    int ri = 0;
    while (true) {
        int64_t here = hts_utell(reader.fp);
        if (here >= seek_end) {
            break;
        }
        int rv = sam_read1(reader.fp, hdr, rec.get());
        if(rv < 0){
            fprintf(stderr, "Something wrong reading %d %lld %lld\n", l, seek_pos, here);
        }
        uint8_t *sd = bam_get_seq(rec.get());
        for(int i = l; i < r; ++i){
            ret[to_sort[i].i] = (char*) malloc(to_sort[i].address.to-to_sort[i].address.from + 1);
            for(int j = to_sort[i].address.from ; j < to_sort[i].address.to; ++j){
                ret[to_sort[i].i][j-to_sort[i].address.from] = seq_nt16_str[bam_seqi(sd, j)];
            }
        }
    }
}
gtl::vector<char *> unitigs2bases(bidbg &dbg, uint64_t uidx,
                                  const BGZFCoreIndex &bgzf,
                                  const lcmer_index_t &lci) {
    const bidbg::value_type &val = *dbg.find(uidx);
    BAMReader reader(bgzf.filename);

    gtl::vector<char *> ret(val.second.merged.size() + LCMER_L, NULL);
    gtl::vector<bgzf_sort_helper> to_sort;
    to_sort.reserve(val.second.merged.size() + LCMER_L);
    const lcmer<LCMER_L> &lbeg = lci.rev.at(val.first);
    int i = 0;
    std::cout << uidx << "\n";
    for (uint32_t v : lbeg.data) {
        // std::cout << "\t" << v << std::endl;
        // std::cout << "\t" << bgzf.idx.at(v).seek << "\n";
        to_sort.emplace_back(i, bgzf.idx.at(v));
        ++i;
    }
    std::sort(to_sort.begin(), to_sort.end(),
              [](const bgzf_sort_helper &a, const bgzf_sort_helper &b) {
                  if (a.address.seek == b.address.seek)
                      return a.address.idx < b.address.idx;
                  return a.address.seek < b.address.seek;
              });
    size_t l = 0;
    size_t r = 0;
    for (; r < to_sort.size(); ++r) {
        if (to_sort[r].address.seek != to_sort[l].address.seek) {
            fill_seqs(to_sort, l, r, bgzf, ret, reader);
            l = r;
        }
    }
    return ret;
}

void convert2basespace(FILE *out, bidbg &dbg, const BGZFCoreIndex &bgzf,
                       const lcmer_index_t &lci) {
    for (const auto &p : dbg) {
        const auto &vec = unitigs2bases(dbg, p.first, bgzf, lci);
        fprintf(out, "S\tU%llu\t", p.first);
        int bc = 0;
        for (char *seg : vec) {
            int printed;
            fprintf(out, "%s%n", seg, &printed);
            bc += printed;
            free(seg);
        }
        fprintf(out, "\tRC:i:%d\tLN:i:%d\n", p.second.cc, bc);
    }
}
int run_assembly(const std::string &input, const std::string &output,
                 int threads, int lcp_level, bloom_param bpd) {
    INIT_TIMERS(HIFI_ASSEMBLY);
    bloom_parameters bparam;
    bparam.false_positive_probability = bpd.false_positive_rate;
    bparam.projected_element_count = bpd.projected_count;
    bparam.random_seed = bpd.seed;

    bparam.compute_optimal_parameters();
    // bloom_filter bfil(bparam);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Bloom Filter setup");

    auto [lci, bgzfci, ovlp, dbg] =
        construct_dbg_v4(input, threads, lcp_level, bparam);
    // // auto [dbg, BGZF_idx, ovlp_idx] =
    // //     construct_dbg_v2(input, threads, lcp_level, bfil);
    // auto [dbg, idx, ridx, BGZF_idx, ovlp_idx] = construct_dbg_v3(input,
    // threads, lcp_level, bfil);qq
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG construction finished");
    // simplify_dbg(dbg);
    simplify_dbg(dbg);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG simplify finished");

    FILE *out = fopen("temp.gfa", "w");
    // auto idx = print_dbg(out, dbg);
    print_dbg_v4(out, dbg, lci);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG printing finished");
    fclose(out);
    out = fopen(output.c_str(), "w");

    convert2basespace(out, dbg, bgzfci, lci);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG to base-space finished");
    fclose(out);
    // auto unitig_graph = unitigs_gfatools(ridx);
    // TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Unitigs finished");

    // to_basespace(out, unitig_graph, BGZF_idx, ovlp_idx);
    // fclose(out);
    // TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Converted to base-space");

    // FILE *topo = fopen("topology.txt", "w");
    // topology(topo, dbg);
    // fclose(topo);
    // TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Topology finished");
    return 0;
}