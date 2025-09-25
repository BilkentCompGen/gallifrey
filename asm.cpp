#include "asm.h"

#include <chrono>
#include <string_view>

// #include <algorithm>
// #include <_strings.h>
#include <concurrentqueue.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include <cstdint>
#include <cstdio>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <thread>
#include <tuple>

#include "bgzf.h"
#include "cascadingbloom.h"
#include "config.h"
#include "ext/bloom_filter.hpp"
#include "ext/kseq.h"
#include "ext/sam_internal.h"
#include "gtl/phmap_fwd_decl.hpp"
#include "lcmer.h"
#include "lcptools/lcptools_ho.h"

using mutex_type = std::mutex;
int hts_useek(htsFile *fp, off_t uoffset, int where) {
    if (fp->is_bgzf)
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
        return (hseek(fp->fp.hfile, uoffset, SEEK_SET) >= 0) ? 0 : -1;
}
off_t hts_utell(htsFile *fp) {
    if (fp->is_bgzf)
        return bgzf_utell(fp->fp.bgzf);
    else
        return htell(fp->fp.hfile);
}
struct bidbg_val {
    gtl::flat_hash_set<uint64_t> prev;
    // gtl::flat_hash_set<uint64_t> next;
    int cc;
    gtl::vector<std::pair<uint64_t, int>> merged;
    char *seq = NULL;
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
    uint64_t, const lcmer<LCMER_L> *, gtl::priv::hash_default_hash<uint64_t>,
    gtl::priv::hash_default_eq<uint64_t>,
    std::allocator<std::pair<const uint64_t, const lcmer<LCMER_L> *>>, 4,
    mutex_type>;
struct lcmer_index_t {
    lcmer_index_map_t map;
    reverse_lcmer_idx_t rev;
    bloom_filter bf;
    lcmer_index_t(const bloom_parameters &bparam) : bf{bparam} {}

    void compute_inverse() {
        for (const auto &[k, v] : map) {
            rev.try_emplace(v, &k);
        }
    }
    uint64_t get(const lcmer<LCMER_L> &l) const {
        auto it = map.find(l);
        if (it != map.end()) {
            return it->second;
        }
        return 0;
    }
    uint64_t emplace(const lcmer<LCMER_L> &l, uint64_t &idx, int tid) {
        auto new_hash = map.hash(l);
        if (!bf.contains(new_hash)) {
            bf.insert(new_hash);
            return 0;
        }
        auto vv = map.try_emplace_with_hash(new_hash, l, idx | tid);
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
                   uint64_t &lcmer_index, bidbg &dbg) {
    lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

    gtl::vector<uint64_t> ll(1, 0);
    ll.reserve(lr.cores.size() - LCMER_L + 3);
    ll.push_back(lci.emplace(l, lcmer_index, thread_id));
    for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
        uint32_t in_core = lr.cores[i];
        l.update_lcmer(in_core);
        ll.push_back(lci.emplace(l, lcmer_index, thread_id));
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
                    ll[i], [&](bidbg::value_type &val) { val.second.cc++; },
                    [&](const bidbg::constructor &ctr) {
                        ctr(ll[i], bidbg_val{{}, 1, {}, NULL});
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
                        ctr(ll[i], bidbg_val{{ll[i - 1]}, 1, {}, NULL});
                    });
            } else {
                dbg.lazy_emplace_l(
                    ll[i],
                    [&](bidbg::value_type &val) {
                        val.second.prev.emplace(ll[i - 1]);
                        val.second.cc++;
                    },
                    [&](const bidbg::constructor &ctr) {
                        ctr(ll[i], bidbg_val{{ll[i - 1]}, 1, {}, NULL});
                    });
            }
        }
    }
}

struct BAMKit {
    /**
     * @brief Generates chunks of SAM/BAM blocks and enqueues them into a
     * concurrent queue.
     *
     * This function instantiates a BAMReader for the given input filename,
     * iterates over its blocks using a BlockIterator, and groups consecutive
     * blocks up to a specified chunk size (BAMREAD_CHUNK_SIZE). It then
     * attempts to enqueue each generated chunk into a
     * moodycamel::ConcurrentQueue of Chunk objects.
     *
     * @param input The path to the input SAM/BAM file.
     * @param chunk_queue A reference to a moodycamel::ConcurrentQueue where
     * chunks will be enqueued.
     */
    static auto generate_chunks(
        std::string input, moodycamel::ConcurrentQueue<Chunk> &chunk_queue) {
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
    static auto process_chunks(const std::string &input, int thread_id,
                               int lcp_level,
                               moodycamel::ConcurrentQueue<Chunk> &chunk_queue,
                               bidbg &dbg, lcmer_index_t &lci,
                               std::atomic<bool> &done,
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
                    std::cerr
                        << "bgzf_seek failed on chunk vstart=" << chunk.vstart
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
                process_cores(lr, lci, thread_id, lcmer_index, dbg);
            }
            all_reads += src;
            // printf("%zu\n", );
        }
    }

    static auto generate_data(
        std::string input, int threads,
        moodycamel::ConcurrentQueue<bam1_t *> &data_queue,
        moodycamel::ConcurrentQueue<bam1_t *> &recycle_queue) {
        samFile *fp = sam_open(input.c_str(), "r");

        if (!fp) {
            std::cerr << "Error: Could not open " << input << std::endl;
            return;
        }
        // Use htslib's thread pool for fast background decompression
        htsThreadPool pool = {NULL, 0};
        if (threads > 1) {
            pool.pool = hts_tpool_init(threads);
            hts_set_thread_pool(fp, &pool);
        }

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
                    std::this_thread::sleep_for(std::chrono::nanoseconds(100));
                }
                rec = nullptr;
            } else {
                // End of file or error. Recycle the unused record and exit.
                bam_destroy1(rec);
                break;
            }
        }
        bam_hdr_destroy(hdr);
        sam_close(fp);
        if (pool.pool) {
            hts_tpool_destroy(pool.pool);
        }
    }
    static auto process_data(
        int thread_id, int lcp_level,
        moodycamel::ConcurrentQueue<bam1_t *> &data_queue,
        moodycamel::ConcurrentQueue<bam1_t *> &recycle_queue, bidbg &dbg,
        lcmer_index_t &lci, std::atomic<bool> &done, int &processed_read_count,
        int &all_reads) {
        uint64_t lcmer_index = 0;
        bam1_t *data_p = NULL;
        while (true) {
            if (!data_queue.try_dequeue(data_p)) {
                if (done.load()) {
                    break;
                }
                std::this_thread::sleep_for(std::chrono::nanoseconds(100));
                continue;
            }
            all_reads++;
            lcp_read lr = parse_cores_from_bam_record(
                bam_get_seq(data_p), data_p->core.l_qseq, lcp_level);
            if (lr.cores.size() < LCMER_L + 2) {
                // printf("Skipping read %lld, cores.size()=%zu\n", here,
                // lr.cores.size());
                goto pd_cleanup;
            }
            ++processed_read_count;
            process_cores(lr, lci, thread_id, lcmer_index, dbg);
        pd_cleanup:
            if (!recycle_queue.try_dequeue(data_p)) {
                bam_destroy1(data_p);
                data_p = NULL;
            }
        }
    }
    static auto reverse_process_cores(uint8_t *seq_nibble, const lcp_read &lr,
                                      const lcmer_index_t &lci, bidbg &dbg) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        uint64_t val = lci.get(l);
        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            int ovlp = lr.pos[i - 1].second - lr.pos[i].first;
            int slen = lr.pos[i].second - lr.pos[i - 1].second;
            dbg.modify_if(val, [&](bidbg::value_type &val) {
                if (val.second.seq == NULL) {
                    if (slen > 0) {
                        val.second.seq = (char *)malloc(slen + 2);
                        nibble2base(seq_nibble + ovlp / 2, val.second.seq,
                                    slen);
                        if (ovlp % 2 == 1) {
                            val.second.seq[0] = 'N';
                        }
                        val.second.seq[slen] = '\0';
                    } else {
                        val.second.seq = (char *)malloc(1);
                        val.second.seq[0] = '\0';
                    }
                }
            });
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            val = lci.get(l);
        }
    }
    static void convert2basespace(std::string input, FILE *out, bidbg &dbg,
                                  const lcmer_index_t &lci, int threads,
                                  int lcp_level) {
        moodycamel::ConcurrentQueue<bam1_t *> data_queue{
            static_cast<size_t>(threads), 1, 1};
        moodycamel::ConcurrentQueue<bam1_t *> recycle_queue{
            static_cast<size_t>(threads), 1, 1};
        std::atomic<bool> done{false};
        gtl::vector<std::thread> workers;
        int data_threads = threads / 3;
        for (int i = 0; i < threads - data_threads; ++i) {
            workers.emplace_back([&]() {
                bam1_t *data_p;
                while (data_queue.size_approx() > 0 || !done.load()) {
                    if (!data_queue.try_dequeue(data_p)) {
                        std::this_thread::sleep_for(
                            std::chrono::nanoseconds(10));
                        continue;
                    }
                    lcp_read lr = parse_cores_from_bam_record(
                        bam_get_seq(data_p), data_p->core.l_qseq, lcp_level);
                    if (lr.cores.size() < LCMER_L + 2) {
                        continue;
                    } else {
                        reverse_process_cores(bam_get_seq(data_p), lr, lci,
                                              dbg);
                    }
                    if (!recycle_queue.try_dequeue(data_p)) {
                        bam_destroy1(data_p);
                        data_p = NULL;
                    }
                }
            });
        }
        generate_data(input, data_threads, data_queue, recycle_queue);
        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }
        bam1_t *item;
        while (recycle_queue.try_dequeue(item)) {
            bam_destroy1(item);
        }
        assert(!data_queue.try_dequeue(item));
    }
};

struct FASTXKit {
    KSEQ_INIT(int, read)

    static auto reverse_process_cores(const char *seq, const lcp_read &lr,
                                      const lcmer_index_t &lci, bidbg &dbg) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        uint64_t val = lci.get(l);
        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            int slen = lr.pos[i].second - lr.pos[i - 1].second;
            dbg.modify_if(val, [&](bidbg::value_type &val) {
                if (val.second.seq == NULL) {
                    if (slen > 0) {
                        val.second.seq = (char *)malloc(slen + 1);
                        strncpy(val.second.seq, seq + lr.pos[i - 1].second,
                                slen);
                        val.second.seq[slen] = '\0';
                    } else {
                        val.second.seq = (char *)malloc(1);
                        val.second.seq[0] = '\0';
                    }
                }
            });
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            val = lci.get(l);
        }
    }
    static auto generate_chunks(
        std::string input, moodycamel::ConcurrentQueue<Chunk> &chunk_queue) {
        FILE *fp = fopen(input.c_str(), "rb");
        fseek(fp, 0L, SEEK_END);
        int64_t sz = ftell(fp);
        fseek(fp, 0L, SEEK_SET);
        // Lets say average read is 20K bases. We want CHUNK_SIZE reads in each
        // chunk. so we want to seek intervals of 20K * CHUNK_SIZE;
        int64_t chunk_size = 20000 * BAMREAD_CHUNK_SIZE;
        int64_t current_seek = 0;
        int64_t prev_seek = 0;
        int64_t new_pos = 0;
        char *buffer = NULL;
        size_t cap = 128000;
        do {
            fseek(fp, chunk_size, SEEK_CUR);
            while (getline(&buffer, &cap, fp) != -1) {
                if (buffer[0] == '@' || buffer[0] == '>') {
                    break;
                }
                new_pos = ftell(fp);
            }

            while (!chunk_queue.try_enqueue(Chunk{prev_seek, new_pos}));
            prev_seek = new_pos;

            current_seek += chunk_size;
        } while (current_seek < sz - chunk_size);
        while (!chunk_queue.try_enqueue(Chunk{prev_seek, sz}));
        free(buffer);
    }
    static auto process_chunks(const std::string &input, int thread_id,
                               int lcp_level,
                               moodycamel::ConcurrentQueue<Chunk> &chunk_queue,
                               bidbg &dbg, lcmer_index_t &lci,
                               std::atomic<bool> &done,
                               int &processed_read_count, int &all_reads) {
        uint64_t lcmer_index = 0;
        Chunk chunk;
        FILE *fp = fopen(input.c_str(), "rb");
        int64_t here = 0;

        kseq_t *seq = kseq_init(fileno(fp));
        while (chunk_queue.size_approx() > 0 || !done.load()) {
            if (!chunk_queue.try_dequeue(chunk)) {
                continue;
            }
            // printf("chunk %lld\t%lld\n", chunk.vstart, chunk.vend);
            fseek(fp, chunk.vstart - here, SEEK_CUR);
            while (kseq_read(seq) > 0) {
                // printf("%s", buffer);

                here = ftell(fp);
                lcp_read lr =
                    parse_cores_from_kseq(seq->seq.s, seq->seq.l, lcp_level);
                ++all_reads;
                if (lr.cores.size() < LCMER_L + 2) {
                    // printf("Skipping read %lld, cores.size()=%zu\n", here,
                    // lr.cores.size());
                    continue;
                }
                ++processed_read_count;
                process_cores(lr, lci, thread_id, lcmer_index, dbg);
                if (here >= chunk.vend) {
                    break;
                }
            }
            // break;
            fflush(stdout);
            // printf("%zu\n", );
        }
        kseq_destroy(seq);
    }

    static void convert2basespace(std::string input, FILE *out, bidbg &dbg,
                                  const lcmer_index_t &lci, int threads,
                                  int lcp_level) {
        moodycamel::ConcurrentQueue<Chunk> chunk_queue{100000, 1, 1};
        std::atomic<bool> done{false};
        gtl::vector<std::thread> workers;
        for (int i = 1; i < threads; ++i) {
            workers.emplace_back([&]() {
                Chunk chunk;
                FILE *fp = fopen(input.c_str(), "rb");
                int64_t here = 0;

                kseq_t *seq = kseq_init(fileno(fp));
                while (chunk_queue.size_approx() > 0 || !done.load()) {
                    if (!chunk_queue.try_dequeue(chunk)) {
                        continue;
                    }
                    // printf("chunk %lld\t%lld\n", chunk.vstart,
                    // chunk.vend);
                    fseek(fp, chunk.vstart - here, SEEK_CUR);
                    while (kseq_read(seq) > 0) {
                        // printf("%s", buffer);

                        here = ftell(fp);
                        lcp_read lr = parse_cores_from_kseq(
                            seq->seq.s, seq->seq.l, lcp_level);
                        if (lr.cores.size() < LCMER_L + 2) {
                            // printf("Skipping read %lld,
                            // cores.size()=%zu\n", here, lr.cores.size());
                            continue;
                        }
                        reverse_process_cores(seq->seq.s, lr, lci, dbg);
                        if (here >= chunk.vend) {
                            break;
                        }
                    }
                    // break;
                    fflush(stdout);
                    // printf("%zu\n", );
                }
                kseq_destroy(seq);
            });
        }
        generate_chunks(input, chunk_queue);
        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }
    }
};
template <class Kit>
auto construct_dbg(std::string input, int threads, int lcp_level,
                   const bloom_parameters &bparam) {
    moodycamel::ConcurrentQueue<Chunk> chunk_queue{100000, 1, 1};
    std::atomic<bool> done{false};

    bidbg dbg;
    lcmer_index_t lci{bparam};
    gtl::vector<std::thread> workers;
    gtl::vector<int> rc(threads - 1, 0);
    gtl::vector<int> all_read_count(threads - 1, 0);
    for (int i = 0; i < threads - 1; ++i) {
        workers.emplace_back(Kit::process_chunks, input, i + 1, lcp_level,
                             std::ref(chunk_queue), std::ref(dbg),
                             std::ref(lci), std::ref(done), std::ref(rc[i]),
                             std::ref(all_read_count[i]));
    }
    Kit::generate_chunks(input, chunk_queue);
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
    return std::make_tuple(std::move(lci), std::move(dbg));
}

template <class Kit>
auto construct_dbgv2(std::string input, int threads, int lcp_level,
                     const bloom_parameters &bparam) {
    moodycamel::ConcurrentQueue<bam1_t *> data_queue{
        static_cast<size_t>(threads), 1, 1};
    moodycamel::ConcurrentQueue<bam1_t *> recycle_queue{
        static_cast<size_t>(threads), 1, 1};
    std::atomic<bool> done{false};

    bidbg dbg;
    lcmer_index_t lci{bparam};
    gtl::vector<std::thread> workers;
    gtl::vector<int> rc(threads - 1, 0);
    gtl::vector<int> all_read_count(threads - 1, 0);
    int data_threads = threads / 3;
    for (int i = 0; i < threads - data_threads; ++i) {
        workers.emplace_back(Kit::process_data, i + 1, lcp_level,
                             std::ref(data_queue), std::ref(recycle_queue),
                             std::ref(dbg), std::ref(lci), std::ref(done),
                             std::ref(rc[i]), std::ref(all_read_count[i]));
    }
    Kit::generate_data(input, data_threads, data_queue, recycle_queue);
    done.store(true);
    for (auto &worker : workers) {
        worker.join();
    }
    Kit::process_data(1, lcp_level, std::ref(data_queue),
                      std::ref(recycle_queue), std::ref(dbg), std::ref(lci),
                      std::ref(done), std::ref(rc[0]),
                      std::ref(all_read_count[0]));
    int all_reads = 0;
    int used_reads = 0;
    for (int i = 0; i < threads - 1; ++i) {
        used_reads += rc[i];
        all_reads += all_read_count[i];
    }
    bam1_t *item;
    while (recycle_queue.try_dequeue(item)) {
        bam_destroy1(item);
    }
    while (data_queue.try_dequeue(item)) {
        bam_destroy1(item);
    }
    fprintf(stderr, "Processed %d reads: Used %d for the dbg and skipped %d\n",
            all_reads, used_reads, all_reads - used_reads);
    return std::make_tuple(std::move(lci), std::move(dbg));
}



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

void print_dbg_seqs(FILE *out, const bidbg &dbg) {
    for (const auto &[lval, edges] : dbg) {
        char *s = edges.seq;
        if (s == NULL) {
            fprintf(out, "S\tC%llu\t%s\tLN:i:%d\tRC:i:%d\n", lval, "", 0,
                    edges.cc);
        } else {
            while (*s == 'N') {
                ++s;
            }
            fprintf(out, "S\tC%llu\t%s\tLN:i:%zu\tRC:i:%d\n", lval, s,
                    strlen(s), edges.cc);
        }

        for (uint64_t nl : edges.prev) {
            if (dbg.contains(nl))
                fprintf(out, "L\tC%llu\t+\tC%llu\t+\t%dM\n", nl, lval, 0);
        }
        free(edges.seq);
    }
}
void print_dbg_v4(FILE *out, bidbg &dbg, lcmer_index_t &lci) {
    for (const auto &[lval, edges] : dbg) {
        int sum_cc = edges.cc;
        for (const auto &p : edges.merged) {
            sum_cc += p.second;
        }
        fprintf(out, "S\tC%llu\t*\tLN:i:%zu\tRC:i:%d\n", lval,
                LCMER_L * (edges.merged.size() + 1) * 100, sum_cc);
        for (uint64_t nl : edges.prev) {
            fprintf(out, "L\tC%llu\t+\tC%llu\t+\t%dM\n", nl, lval,
                    100 * (LCMER_L - 1));
        }
        // for (uint64_t nl : edges.next) {
        //     fprintf(out, "L\tC%llu\t-\tC%llu\t-\t0M\n", lval, nl);
        // }
    }
}

enum class input_type { BAM, FASTX, FASTXGZ, FASTXBZ2, UNKNOWN };

input_type identify_input(const std::string &input) {
    std::string_view sv(input);
    if (sv.ends_with(".bam")) {
        return input_type::BAM;
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

// bool test_bamseq_conversion(const std::string &input) {
//     BAMReader reader(input);
//     static constexpr int tab[] = {4, 0, 1, 4, 2, 4, 4, 4,
//                                     3, 4, 4, 4, 4, 4, 4, 4};
//     std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(),
//                                                          bam_destroy1);

//     int ret = sam_read1(reader.fp, NULL, rec.get());
//     if(ret < 0){
//         fprintf(stderr, "Error reading BAM file: %s\n", input.c_str());
//         return false;
//     }
//     uint8_t *converted = onehotto2bitDNA(bam_get_seq(rec.get()),
//     rec->core.l_qseq); for(int i = 0; i < rec->core.l_qseq; i++){
//         if(converted[i] !=  tab[bam_seqi(bam_get_seq(rec.get()), i)]){
//             fprintf(stderr, "Diff %d: %d %d\n", i, converted[i],
//             tab[bam_seqi(bam_get_seq(rec.get()), i)]);
//             // return false;
//         }
//     }
//     fprintf(stderr, "Good!\n");
//     return true;
// }

int run_assembly(const std::string &input, const std::string &output,
                 int threads, int lcp_level, bloom_param bpd) {
    INIT_TIMERS(HIFI_ASSEMBLY);


    // return test_bamseq_conversion(input);
    bloom_parameters bparam;
    bparam.false_positive_probability = bpd.false_positive_rate;
    bparam.projected_element_count = bpd.projected_count;
    bparam.random_seed = bpd.seed;
    bparam.compute_optimal_parameters();

    // bloom_filter bfil(bparam);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Bloom Filter setup");

    auto inty = identify_input(input);

    if (inty == input_type::UNKNOWN) {
        fprintf(stderr, "Unknown input type %s\n", input.c_str());
        return -1;
    }
    auto [lci, dbg] = [&](const input_type &i) {
        switch (i) {
            case input_type::BAM:
                return construct_dbgv2<BAMKit>(input, threads, lcp_level,
                                               bparam);
            case input_type::FASTX:
                return construct_dbg<FASTXKit>(input, threads, lcp_level,
                                               bparam);
            default:
                return construct_dbg<FASTXKit>(input, threads, lcp_level,
                                               bparam);
                fprintf(stderr, "Unknown input type %s\n", input.c_str());
                // return construct_dbg<BAMKit>(input, threads, lcp_level,
                // bparam);
        }
    }(inty);

    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG construction finished");
    // simplify_dbg(dbg);

    // // TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG simplify finished");
    // {
    //     FILE *out = fopen("temp.gfa", "w");
    //     // auto idx = print_dbg(out, dbg);
    //     print_dbg_v4(out, dbg, lci);
    //     TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG printing finished");
    //     fclose(out);
    // }

    auto result = [&](input_type i) {
        FILE *out = fopen(output.c_str(), "w");
        switch (i) {
            case input_type::BAM:
                BAMKit::convert2basespace(input, out, dbg, lci, threads,
                                          lcp_level);
                break;
            case input_type::FASTX:
                // lci.compute_inverse();
                FASTXKit::convert2basespace(input, out, dbg, lci, threads,
                                            lcp_level);
                break;
            default:
                fclose(out);
                fprintf(stderr, "Unknown input type %d\n", i);
                return -1;
        }
        print_dbg_seqs(out, dbg);
        fclose(out);
        TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG to base-space finished");
        return 0;
    }(inty);
    return result;
}

// // TESTS
// #include "ext/doctest.h"
// #include "gtl/btree.hpp"
// #define CHUNK_SIZE 256

// // Build chunks from raw blocks, preserving exact end offsets
// std::vector<Chunk> indexChunks(const std::string &fn, off_t off) {
//     BlockIteratorN it(fn, off);
//     std::vector<Block> blocks;
//     Block b;
//     while (it.next(b)) blocks.push_back(b);

//     std::vector<Chunk> chunks;
//     size_t n = blocks.size();
//     for (size_t i = 0; i < n; i += CHUNK_SIZE) {
//         size_t end_idx = std::min(i + CHUNK_SIZE, n) - 1;
//         chunks.push_back({blocks[i].vstart(), blocks[end_idx].vend()});
//     }
//     return chunks;
// }

// void processChunk(const std::string &fn, bam_hdr_t *hdr, const Chunk &c,
//                   size_t &cnt, gtl::btree_multiset<int> &lens) {
//     htsFile *fp = hts_open(fn.c_str(), "r");
//     if (!fp) return;
//     bgzf_seek(fp->fp.bgzf, c.vstart, SEEK_SET);
//     auto rec = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>(bam_init1(),
//                                                                 bam_destroy1);
//     while (true) {
//         int64_t start = bgzf_tell(fp->fp.bgzf);
//         // stop if next record would start beyond the chunk
//         if (start >= c.vend) break;
//         int ret = sam_read1(fp, hdr, rec.get());
//         if (ret < 0) break;
//         ++cnt;
//         lens.insert(rec->core.l_qseq);
//     }
//     hts_close(fp);
// }
// TEST_CASE("Full coverage test") {
//     std::string fn = "../tests/test_small.bam";
//     htsFile *fp = hts_open(fn.c_str(), "r");
//     REQUIRE(fp);
//     auto hdr = std::unique_ptr<bam_hdr_t, decltype(&sam_hdr_destroy)>(
//         sam_hdr_read(fp), sam_hdr_destroy);
//     size_t exp = 0;
//     gtl::btree_multiset<int> exl;
//     auto r = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>(bam_init1(),
//                                                               bam_destroy1);
//     while (sam_read1(fp, hdr.get(), r.get()) >= 0) {
//         ++exp;
//         exl.insert(r->core.l_qseq);
//     }
//     hts_close(fp);

//     // chunks and process
//     fp = hts_open(fn.c_str(), "r");
//     hdr.reset(sam_hdr_read(fp));
//     off_t off = bgzf_tell(fp->fp.bgzf) >> 16;
//     hts_close(fp);
//     auto ch = indexChunks(fn, off);
//     size_t pr = 0;
//     gtl::btree_multiset<int> prl;
//     for (auto &c : ch) processChunk(fn, hdr.get(), c, pr, prl);

//     CHECK(exp == pr);
//     CHECK(exl == prl);
// }
