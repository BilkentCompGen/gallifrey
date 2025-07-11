#include "asm.h"

#include <concurrentqueue.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/sam.h>

#include <cstdio>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <iostream>

#include "bgzf.h"
#include "endian.h"  // for le16toh()
#include "ext/bloom_filter.hpp"
#include "gtl/phmap_fwd_decl.hpp"
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

struct bgzfseqaddress {
    int64_t seek;
    uint32_t idx;
    uint32_t from;
    uint32_t to;
};
using core2bgzf_index_t = gtl::parallel_flat_hash_map<
    uint32_t, bgzfseqaddress, gtl::priv::hash_default_hash<uint32_t>,
    gtl::priv::hash_default_eq<uint32_t>,
    std::allocator<std::pair<const uint32_t, bgzfseqaddress>>, 4, mutex_type>;

struct BGZFCoreIndex {
    core2bgzf_index_t idx;
    void emplace(uint32_t core, bgzfseqaddress a) { idx.try_emplace(core, a); }
};

inline void emplace2dbg_first(dbg_map_t &dbg, lcmer<LCMER_L> &l,
                              uint64_t &lcmer_index, int thread_id) {
    dbg.lazy_emplace_l(
        l, [](dbg_map_t::value_type &val) { val.second.count++; },
        [&](const dbg_map_t::constructor &ctor) {
            ctor(l, dbg_value_type{1UL, lcmer_index | thread_id, {}});
            lcmer_index += 128;
        });
}

inline void emplace2dbg_linked(dbg_map_t &dbg, lcmer<LCMER_L> &l,
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

struct Chunk {
    int64_t vstart;
    int64_t vend;
};

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

dbg_map_t construct_dbg_v2(std::string input, int threads, int lcp_level,
                           bloom_filter &bfil) {
    using worktype = Chunk;
    moodycamel::ConcurrentQueue<worktype> work_queue(32 * threads, 1, 0);
    std::atomic<bool> done{false};
    gtl::vector<int64_t> skipped_reads_count(threads, 0);
    gtl::vector<int64_t> count(threads, 0);
    dbg_map_t dbg;

    BGZFCoreIndex core_idx;

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
            unsigned sr1_cnt = 0;
            while (true) {
                int64_t here = bgzf_tell(reader.fp->fp.bgzf);
                if (here >= chunk.vend) {
                    break;
                }
                int ret = sam_read1(reader.fp, hdr, rec.get());
                ++sr1_cnt;
                if (ret == -1) {
                    fprintf(stderr, "End of file %s at %lld %lld\n",
                            reader.filename.c_str(), chunk.vstart, chunk.vend);
                    break;
                } else if (ret < -1) {
                    fprintf(
                        stderr,
                        "sam_read1 failed on %s %lld %lld possibly due to a "
                        "read "
                        "being longer than block size in the previous block\n",
                        reader.filename.c_str(), chunk.vstart, chunk.vend);
                    break;
                }
                lcp_read lr = parse_cores_from_bam_record(
                    bam_get_seq(rec.get()), rec->core.l_qseq, lcp_level);

                if (lr.cores.size() < LCMER_L + 2) {
                    skipped_reads_count[tid] += 1;
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
                for (int q = i - LCMER_L; q < i; ++q) {
                    core_idx.emplace(lr.cores[q],
                                     {chunk.vstart, sr1_cnt, lr.pos[q].first,
                                      lr.pos[q].second});
                }
                std::optional<lcmer<LCMER_L>> waiter = l;

                for (; i < lr.cores.size(); ++i) {
                    uint32_t out_core = l.update_lcmer(lr.cores[i]);
                    if (!bfil.contains(l)) {
                        bfil.insert(l);
                        skipped = true;
                        continue;
                    }
                    core_idx.emplace(lr.cores[i],
                                     {chunk.vstart, sr1_cnt, lr.pos[i].first,
                                      lr.pos[i].second});
                    if (skipped) {
                        waiter = l;
                        skipped = false;
                        continue;
                    }
                    if (waiter.has_value()) {  // TODO Validate this approach
                        emplace2dbg_first(dbg, waiter.value(), lcmer_index,
                                          tid);
                        waiter = std::nullopt;
                    }
                    emplace2dbg_linked(dbg, l, lcmer_index, skipped, out_core,
                                       tid);
                }
                skipped = false;
                count[tid] += 1;
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
        while (count < BAMREAD_CHUNK_SIZE && iter.next(blk)) {
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
    int64_t count_accu = 0;
    for(int64_t c : count){
        count_accu += c;
    }
    int64_t skip_accu = 0;
    for(int64_t s : skipped_reads_count){
        skip_accu += s;
    }
    printf("Processed %lld reads\n", count_accu);
    printf("Skipped %lld reads due to insufficient cores\n",
           skip_accu);
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

int run_assembly(const std::string &input, const std::string &output,
                 int threads, int lcp_level, bloom_param bpd) {
    INIT_TIMERS(HIFI_ASSEMBLY);
    bloom_parameters bparam;
    bparam.false_positive_probability = bpd.false_positive_rate;
    bparam.projected_element_count = bpd.projected_count;
    bparam.random_seed = bpd.seed;

    bparam.compute_optimal_parameters();
    bloom_filter bfil(bparam);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "Bloom Filter setup");

    dbg_map_t dbg = construct_dbg_v2(input, threads, lcp_level, bfil);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG construction finished");
    FILE *out = fopen(output.c_str(), "w");
    print_dbg(out, dbg);
    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG printing finished");

    fclose(out);
    return 0;
}