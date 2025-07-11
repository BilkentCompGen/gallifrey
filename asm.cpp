#include "asm.h"
#include "bgzf.h"

#include <concurrentqueue.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>

#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <iostream>

#include "endian.h"  // for le16toh()
#include "lcmer.h"

struct dbg_value {
    size_t count;
    size_t index;
    gtl::flat_hash_set<uint32_t> next;
};

using mutex_type = std::mutex;

using dbg_value_type = dbg_value;
using dbg_map_t = gtl::parallel_flat_hash_map<
    lcmer<LCMER_L>, dbg_value_type, std::hash<lcmer<LCMER_L>>,
    gtl::priv::hash_default_eq<lcmer<LCMER_L>>,
    std::allocator<std::pair<const lcmer<LCMER_L>, dbg_value_type>>, 4,
    mutex_type>;

// move the dbg_map emplace lambda to a separate function.
inline void emplace2dbg(dbg_map_t &dbg, lcmer<LCMER_L> &l,
                        uint64_t &lcmer_index, int thread_id) {
    dbg.lazy_emplace_l(
        l, [](dbg_map_t::value_type &val) { val.second.count++; },
        [&](const dbg_map_t::constructor &ctor) {
            ctor(l, dbg_value_type{1UL, lcmer_index | thread_id, {}});
            lcmer_index += 128;
        });
}
inline void emplace2dbg2(dbg_map_t &dbg, lcmer<LCMER_L> &l,
                         uint64_t &lcmer_index, bool &skipped,
                         uint32_t out_core, int thread_id) {
    dbg.lazy_emplace_l(
        l,
        [&](dbg_map_t::value_type &val) {
            val.second.count++;
            if (!skipped) val.second.next.insert(out_core);
        },
        [&](const dbg_map_t::constructor &ctor) {
            if (skipped)
                ctor(l, dbg_value_type{1UL, lcmer_index | thread_id, {}});
            else
                ctor(l,
                     dbg_value_type{1UL, lcmer_index | thread_id, {out_core}});
            lcmer_index += 128;
        });
}

dbg_map_t construct_dbg(std::string input, int threads, int lcp_level,
                        bloom_filter &bfil) {
    std::atomic<uint64_t> count;
    using worktype = std::tuple<uint64_t, uint8_t *>;

    moodycamel::ConcurrentQueue<worktype> work_queue(32 * threads, 1, 0);

    std::vector<std::thread> consumers;
    //   counter core_counter;
    dbg_map_t dbg;
    std::atomic<size_t> skipped_reads_count{0};
    std::atomic<bool> done{false};
    for (int tid = 0; tid < threads; ++tid) {
        consumers.emplace_back(
            [&](int t) {
                uint64_t lcmer_index = 0;
                worktype work;
                while (!done.load() || work_queue.size_approx() > 0) {
                    if (!work_queue.try_dequeue(work)) {
                        continue;
                    }
                    lcp_read lr = parse_cores_from_bam_record(
                        std::get<1>(work), std::get<0>(work), lcp_level);
                    // cout << bamdata->core.l_qseq << " to " << lr.cores.size()
                    // <<
                    // "\n";
                    // core_counter[lr.cores.size()]++;
                    if (lr.cores.size() < LCMER_L + 2) {
                        skipped_reads_count.fetch_add(1);
                        free(std::get<1>(work));
                        continue;
                    }
                    lcmer<LCMER_L> l{lr.cores.begin()};
                    bool skipped = false;
                    size_t i;
                    for (i = LCMER_L; i < lr.cores.size(); ++i) {
                        if (bfil.contains(l)) {
                            break;
                        }
                        bfil.insert(l);
                        l.update_lcmer(lr.cores[i]);
                    }
                    dbg.lazy_emplace_l(
                        l,
                        [](dbg_map_t::value_type &val) { val.second.count++; },
                        [&](const dbg_map_t::constructor &ctor) {
                            ctor(l, dbg_value_type{1UL, lcmer_index | t, {}});
                            lcmer_index += 128;
                        });
                    for (; i < lr.cores.size(); ++i) {
                        uint32_t out_core = l.update_lcmer(lr.cores[i]);
                        if (!bfil.contains(l)) {
                            bfil.insert(l);
                            skipped = true;
                            continue;
                        }

                        dbg.lazy_emplace_l(
                            l,
                            [&](dbg_map_t::value_type &val) {
                                val.second.count++;
                                if (!skipped) val.second.next.insert(out_core);
                            },
                            [&](const dbg_map_t::constructor &ctor) {
                                if (skipped)
                                    ctor(l, dbg_value_type{
                                                1UL, lcmer_index | t, {}});
                                else
                                    ctor(l,
                                         dbg_value_type{
                                             1UL, lcmer_index | t, {out_core}});
                                lcmer_index += 128;
                            });
                        skipped = false;
                        count.fetch_add(1);
                    }
                    free(std::get<1>(work));
                }
            },
            tid);
    }
    uint8_t *data = NULL;

    BAMReader reader(input);
    bam1_t *bamdata = reader.init_data();

    while (reader.read(bamdata)) {
        data = bam_get_seq(bamdata);
        uint8_t *data_copy = (uint8_t *)malloc(bamdata->core.l_qseq / 2 +
                                               1);  // TODO Fix the extra memory
        memcpy(data_copy, data, bamdata->core.l_qseq / 2 + 1);
        while (!work_queue.try_enqueue(
            std::make_tuple(bamdata->core.l_qseq, data_copy)));
    }
    done.store(true);
    for (std::thread &t : consumers) {
        t.join();
    }
    printf("Processed %llu reads\n", count.load());
    printf("Skipped %lu reads due to insufficient cores\n",
           skipped_reads_count.load());
    //   for (auto &p : core_counter) {
    //     fprintf(stderr, "CC\t%lu\t%lu\n", p.first, p.second);
    //   }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    return dbg;
}

struct Chunk {
    int64_t vstart;
    int64_t vend;
};

constexpr size_t CHUNK_SIZE = 128;

struct Block {
    off_t start;
    uint32_t len;
    constexpr int64_t vstart() const noexcept {
        return static_cast<int64_t>(start) << 16;
    }
    constexpr int64_t vend() const noexcept {
        return (static_cast<int64_t>(start + len)) << 16;
    }
};

class BlockIterator {
   public:
    BlockIterator(htsFile *hts_fp) {
        bgzf_ = hts_fp->fp.bgzf;
        hfp_ = bgzf_->fp;
        if (hseek(hfp_, 0, SEEK_SET) < 0) {
            std::perror("hseek");
            valid_ = false;
        }
    }

    bool next(Block &block) {
        if (!valid_) return false;
        // use htell to get current hFILE offset
        off_t start = htell(hfp_);
        unsigned char hdr[18];
        ssize_t got = hread(hfp_, hdr, sizeof(hdr));
        if (got < static_cast<ssize_t>(sizeof(hdr))) return false;
        if (hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 0x08 ||
            hdr[3] != 0x04)
            return false;
        uint16_t bsize = le16toh(*reinterpret_cast<uint16_t *>(hdr + 16));
        uint32_t block_len = bsize + 1;
        block.start = start;
        block.len = block_len;
        // skip raw file pointer to next block
        if (hseek(hfp_, start + block_len, SEEK_SET) < 0) {
            valid_ = false;
            return false;
        }
        return true;
    }

   private:
    BGZF *bgzf_ = nullptr;
    hFILE *hfp_ = nullptr;
    bool valid_ = true;
};

// void process_chunk(htsFile *fp, const Chunk &chunk) {

//     if (chunk.vstart > 0) {
//         if (bgzf_seek(fp->fp.bgzf, chunk.vstart, SEEK_SET) < 0) {
//             std::cerr << "bgzf_seek failed on chunk vstart=" << chunk.vstart
//                       << "\n";
//             std::exit(EXIT_FAILURE);
//         }
//     }
//     std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(),
//                                                          bam_destroy1);
//     while (true) {
//         int64_t here = bgzf_tell(fp->fp.bgzf);
//         if (here >= chunk.vend) {
//             // std::print(std::cerr, "Bopped {} {}\n", chunk.vstart,
//             // chunk.vend);
//             break;
//         }
//         int ret = sam_read1(fp, hdr, rec.get());
//         if (ret < 0) {
//             std::print(std::cerr, "Popped {} {}\n", chunk.vstart,
//             chunk.vend); break;
//         }
//         std::print(cout, "{}\t{}\n", string{bam_get_qname(rec.get())},
//                    here << '\n');
//     }
// }

dbg_map_t construct_dbg_v2(std::string input, int threads, int lcp_level,
                           bloom_filter &bfil) {
    using worktype = Chunk;
    moodycamel::ConcurrentQueue<worktype> work_queue(32 * threads, 1, 0);
    std::atomic<bool> done{false};
    std::atomic<int64_t> skipped_reads_count{0};
    std::atomic<int64_t> count{0};
    dbg_map_t dbg;
    auto work = [&](int tid) {
        uint64_t lcmer_index = 0;
        worktype chunk;
        BAMReader reader{input};
        uint8_t *seqdata;
        std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(),
                                                             bam_destroy1);
        sam_hdr_t *hdr = NULL;
        while (!done.load() || work_queue.size_approx() > 0) {
            if (!work_queue.try_dequeue(chunk)) {
                continue;
            }
            if (chunk.vstart > 0) {
                if (bgzf_seek(reader.fp->fp.bgzf, chunk.vstart, SEEK_SET) < 0) {
                    std::cerr
                        << "bgzf_seek failed on chunk vstart=" << chunk.vstart
                        << "\n";
                    std::exit(EXIT_FAILURE);
                }
            } else {
                // hdr = sam_hdr_read(reader.fp);
            }

            while (true) {
                int64_t here = bgzf_tell(reader.fp->fp.bgzf);
                if (here >= chunk.vend) {
                    // std::print(std::cerr, "Bopped {} {}\n", chunk.vstart,
                    // chunk.vend);
                    break;
                }
                int ret = sam_read1(reader.fp, hdr, rec.get());
                if (ret < 0) {
                    std::cerr << std::format("Popped {} {}\n", chunk.vstart,
                                             chunk.vend);
                    break;
                }

                lcp_read lr = parse_cores_from_bam_record(
                    bam_get_seq(rec.get()), rec->core.l_qseq, lcp_level);

                if (lr.cores.size() < LCMER_L + 2) {
                    skipped_reads_count.fetch_add(1);
                    continue;
                }
                lcmer<LCMER_L> l{lr.cores.begin()};
                bool skipped = false;
                size_t i;
                for (i = LCMER_L; i < lr.cores.size(); ++i) {
                    if (bfil.contains(l)) {
                        break;
                    }
                    bfil.insert(l);
                    l.update_lcmer(lr.cores[i]);
                }
                std::optional<lcmer<LCMER_L>> waiter = l;
                // uint32_t wait_core = out_core;
                // emplace2dbg(dbg, l, lcmer_index, tid);
                for (; i < lr.cores.size(); ++i) {
                    uint32_t out_core = l.update_lcmer(lr.cores[i]);
                    if (!bfil.contains(l)) {
                        bfil.insert(l);
                        skipped = true;
                        continue;
                    }
                    if (skipped) {
                        waiter = l;
                        skipped = false;
                        continue;
                    }
                    if (waiter.has_value()) {  // TODO Validate this approach
                        emplace2dbg(dbg, waiter.value(), lcmer_index, tid);
                        waiter = std::nullopt;
                    }
                    emplace2dbg2(dbg, l, lcmer_index, skipped, out_core, tid);
                }
                skipped = false;
                count.fetch_add(1);
                // std::print(cout, "{}\t{}\n",
                // string{bam_get_qname(rec.get())},
                //    here << '\n');
            }
        }
    };
    gtl::vector<std::thread> consumers;
    for (int i = 0; i < threads; ++i) {
        consumers.emplace_back(work, i);
    }
    BAMReader reader{input};
    BlockIterator iter(reader.fp);

    Block blk;
    // fetch chunks until no more blocks
    while (iter.next(blk)) {
        // first block in chunk
        int64_t first_vstart = blk.vstart();
        int64_t last_vend = blk.vend();

        size_t count = 1;
        // accumulate up to CHUNK_SIZE blocks
        while (count < CHUNK_SIZE && iter.next(blk)) {
            last_vend = blk.vend();
            ++count;
        }
        while (!work_queue.try_enqueue({first_vstart, last_vend}));
        // chunks.push_back({first_vstart, last_vend});
    }
    done.store(true);
    for (auto &t : consumers) {
        t.join();
    }
    printf("Processed %lld reads\n", count.load());
    printf("Skipped %lld reads due to insufficient cores\n",
           skipped_reads_count.load());
    return dbg;
}

void print_dbg(FILE *out, const dbg_map_t &dbg) {
    for (const auto &p : dbg) {
        fprintf(out, "S\t%lu\t*\tLN:i:%d\tRC:i:%lu\n", p.second.index, LCMER_L,
                p.second.count);
    }
    for (const auto &p : dbg) {
        lcmer ln = lcmer<LCMER_L>::right_shift(p.first);
        for (uint32_t next : p.second.next) {
            ln.setN(next);
            if (dbg.find(ln) == dbg.end()) {
                fprintf(stderr, "Cannot find the lcmer!\n");
                continue;
            }
            fprintf(out, "L\t%lu\t+\t%lu\t+\t%dM\n", dbg.at(ln).index,
                    p.second.index, LCMER_L - 1);
        }
    }
}

int run_assembly(const std::string &input, const std::string &output, int threads, int lcp_level, bloom_parameters bparam) {
    INIT_TIMERS(HIFIASM);
    bparam.compute_optimal_parameters();
    bloom_filter bfil(bparam);
    TIME_CHECKPOINT(HIFIASM, stderr, "Bloom Filter setup");

    dbg_map_t dbg = construct_dbg_v2(input, threads, lcp_level, bfil);
    TIME_CHECKPOINT(HIFIASM, stderr, "DBG construction finished");
    FILE *out = fopen(output.c_str(), "w");
    print_dbg(out, dbg);
    TIME_CHECKPOINT(HIFIASM, stderr, "DBG printing finished");

    fclose(out);
    return 0;
}