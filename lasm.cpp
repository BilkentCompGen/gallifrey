#include <CLI/CLI.hpp>
#include <cstdint>

#include "CLI/CLI.hpp"
#define LCPTOOLS_IMPL
#include <concurrentqueue.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <string>

#include "config.h"
#include "endian.h"  // for le16toh()
#include "ext/bloom_filter.hpp"
#include "ext/lcptools/lcptools_ho.h"
#include "histogram.h"
#include "lcmer.h"
#include "syncmer.h"

using std::string;

struct core_view {
    core *from;
    core *to;
};

// https://github.com/ruanjue/wtdbg2/blob/b77c5657c8095412317e4a20fe3668f5bde6b1ac/dna.h#L108
template <class T>
T reverse_complement_kmer(T kmer, int K) {
    uint32_t shift_count = (sizeof(T) * 8) - 2 * K;
    kmer = ~kmer;
    kmer = ((kmer & 0x3333333333333333LLU) << 2) |
           ((kmer & 0xCCCCCCCCCCCCCCCCLLU) >> 2);
    kmer = ((kmer & 0x0F0F0F0F0F0F0F0FLLU) << 4) |
           ((kmer & 0xF0F0F0F0F0F0F0F0LLU) >> 4);
    if constexpr (sizeof(T) == 8) {
        kmer = __builtin_bswap64(kmer);
    } else if constexpr (sizeof(T) == 4) {
        kmer = __builtin_bswap32(kmer);
    } else if constexpr (sizeof(T) == 2) {
        kmer = __builtin_bswap16(kmer);
    }
    kmer = kmer >> shift_count;
    return kmer;
}

template <int _K, typename INT_TYPE = uint32_t>
struct KmerViewOnNibbles {
    using int_type = INT_TYPE;
    static constexpr int k = _K;
    int_type value;
    static constexpr int_type KMASK = ((int_type)1 << k) - 1;
    uint8_t nibble_mask = 15;
    KmerViewOnNibbles() : value{0} {}

    int fill(const uint8_t *data_address) {
        int i;
        for (i = 0; i < k; i++) {
            move_get(data_address, i);
        }
        return i;
    }
    int_type move_get(const uint8_t *data_address, int index) {
        value <<= 2;
        if constexpr (KMASK != -1) value &= KMASK;
        value |= bam_seqi(data_address, index);
        return value;
    }
};

void init_core1bam(struct core *cr, const uint8_t *data, uint64_t distance,
                   uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= ((distance - 2) << 6);
    cr->label |= (bam_seqi(data, si) << 4);
    cr->label |= (bam_seqi(data, si + distance - 2) << 2);
    cr->label |= (bam_seqi(data, si + distance - 1));
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

template <typename _T>
void init_core_kmer(struct core *cr, _T Kval, uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= Kval;
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * SYNCMER_K;
}

template <typename _T>
void init_core_kmer_xor(struct core *cr, _T Kval, uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= (Kval ^ reverse_complement_kmer(Kval, SYNCMER_K));
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * SYNCMER_K;
}
syncmer::syncmer_precomputed<
    syncmer::SYNCMER_RULE<SYNCMER_K, SYNCMER_S, SYNCMER_T>>
    table;
uint64_t parse1_syncmers(const uint8_t *data, int64_t len, struct core *cores) {
    int core_index = 0;
    KmerViewOnNibbles<SYNCMER_K, uint32_t> kmer_view;
    int64_t it1 = kmer_view.fill(data);
    // find lcp cores
    for (; it1 + SYNCMER_K < len; it1++) {
        if (table[kmer_view.value]) {
            init_core_kmer<uint32_t>(&(cores[core_index]), kmer_view.value, it1,
                                     it1 + SYNCMER_K);
            core_index++;
        }
        kmer_view.move_get(data, it1);
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}

uint64_t parse1_bam(const uint8_t *data, int64_t len, struct core *cores) {
    int64_t it1 = 0;
    int64_t it2 = len;

    int core_index = 0;
    // int64_t last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < len; it1++) {
        if (bam_seqi(data, it1) == bam_seqi(data, it1 + 1)) {
            continue;
        }

        // // check for RINT core
        // if (bam_seqi(data, (it1 + 1)) == bam_seqi(data, (it1 + 2))) {

        //   // count middle characters
        //   uint32_t middle_count = 1;
        //   int64_t temp = it1 + 2;
        //   while (temp < len && bam_seqi(data, (temp - 1)) == bam_seqi(data,
        //   temp)) {
        //     temp++;
        //     middle_count++;
        //   }
        //   if (temp != len) {
        //     // check if there is any SSEQ cores left behind
        //     if (it2 < it1 && last_invalid_char_index < it2 - 1) {
        //       init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 -
        //       1,
        //                     it1 + 1);
        //       core_index++;
        //     }

        //     // create RINT core
        //     it2 = it1 + 2 + middle_count;
        //     init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
        //     core_index++;

        //     continue;
        //   }
        // }

        // if (bam_seqi(data, it1) > bam_seqi(data, (it1 + 1)) &&
        //     bam_seqi(data, (it1 + 1)) < bam_seqi(data, (it1 + 2))) {

        //   // check if there is any SSEQ cores left behind
        //   if (it2 < it1 && last_invalid_char_index < it2 - 1) {
        //     init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
        //                   it1 + 1);
        //     core_index++;
        //   }

        //   // create LMIN core
        //   it2 = it1 + 3;
        //   init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
        //   core_index++;

        //   continue;
        // }

        if (0 == it1) {
            continue;
        }

        // check for LMAX
        if (it1 + 3 < len && bam_seqi(data, it1) < bam_seqi(data, (it1 + 1)) &&
            bam_seqi(data, (it1 + 1)) > bam_seqi(data, (it1 + 2)) &&
            bam_seqi(data, (it1 - 1)) <= bam_seqi(data, (it1)) &&
            bam_seqi(data, (it1 + 2)) >= bam_seqi(data, (it1 + 3))) {
            // // check if there is any SSEQ cores left behind
            // if (it2 < it1 && last_invalid_char_index < it2 - 1) {
            //   init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 -
            //   1,
            //                 it1 + 1);
            //   core_index++;
            // }

            // create LMAX core
            it2 = it1 + 3;
            init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
            core_index++;

            continue;
        }
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}

/// @brief Initializes an lps struct with a given sequence in the bam compressed
/// format
void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len) {
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores =
        (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
    lps_ptr->size = parse1_syncmers(data, len, lps_ptr->cores);
}

struct lcp_read {
    gtl::vector<uint32_t> cores;
    gtl::vector<std::pair<uint32_t, uint32_t>> pos;
    void reserve(size_t size) {
        cores.reserve(size);
        pos.reserve(size);
    }
};
struct dbg_value {
    size_t count;
    size_t index;
    gtl::flat_hash_set<uint32_t> next;
};
lcp_read parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                     int lcp_level) {
    // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

    lps str;
    init_lps_bam(&str, data, len);
    lps_deepen(&str, lcp_level);
    lcp_read lr;
    lr.reserve(str.size);
    for (int i = 0; i < str.size; ++i) {
        // fprintf(stderr, "%u\n", str.cores[i].bit_size);
        lr.cores.emplace_back(str.cores[i].label);
        lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
    }
    free_lps(&str);
    return lr;
}
struct BAMReader {
    string filename;
    samFile *fp;
    sam_hdr_t *hdr;

    BAMReader(const string &filename, int threads = 2)
        : fp{sam_open(filename.c_str(), "rb")} {
        this->filename = filename;
        // bgzf_mt(fp->fp.bgzf, 2, 0);
        // hts_set_threads(fp,threads);
        if (!(hdr = sam_hdr_read(fp))) {
            throw "Cannot read header from file!";
        }
    }
    ~BAMReader() {
        if (hdr) {
            sam_hdr_destroy(hdr);
        }
        if (fp) {
            sam_close(fp);
        }
    }

    bool read(bam1_t *bamdata) { return sam_read1(fp, hdr, bamdata) >= 0; }
    off_t tell() {
        if (fp->is_bgzf) {
            return bgzf_tell(fp->fp.bgzf);
        }
        return -199;
    }
    int seek(off_t offset) {
        if (fp->is_bgzf) {
            // return hseek(fp->fp.bgzf->fp, offset, SEEK_SET);
            return bgzf_seek(fp->fp.bgzf, offset, SEEK_SET);
        } else
            return -199;
    }
    bam1_t *init_data() {
        bam1_t *data = bam_init1();
        if (data == NULL) {
            throw "Not enough memory";
        }
        return data;
    }
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

dbg_map_t construct_dbg(string input, int threads, int lcp_level,
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

dbg_map_t construct_dbg_v2(string input, int threads, int lcp_level,
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

                emplace2dbg(dbg, l, lcmer_index, tid);
                for (; i < lr.cores.size(); ++i) {
                    uint32_t out_core = l.update_lcmer(lr.cores[i]);
                    if (!bfil.contains(l)) {
                        bfil.insert(l);
                        skipped = true;
                        continue;
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

int asm_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("hifi", "run hifi assembly");
    auto &app = *app_;
    std::shared_ptr<string> input{new string};
    std::shared_ptr<string> output{new string};
    std::shared_ptr<int> lcp_level{new int};
    std::shared_ptr<int> threads{new int{8}};
    app.add_option("input", *input, "Input bam file")->required();
    app.add_option("output", *output, "Output gfa file")->required();
    app.add_option("-l,--level", *lcp_level, "LCP level to compute")
        ->required();
    app.add_option("-t,--threads", *threads, "Number of threads")
        ->capture_default_str();

    std::shared_ptr<bloom_parameters> bparam{new bloom_parameters};
    bparam->projected_element_count = 1'000'000;
    bparam->false_positive_probability = 0.0001;
    bparam->random_seed = 0xA5A5A5A5;
    auto bfg =
        app.add_option_group("bloomfilter", "Parameters to setup bloom filter");
    bfg->add_option("--projected-element-count",
                    bparam->projected_element_count, "Projected element count")
        ->capture_default_str();
    bfg->add_option("--false-positive-probability",
                    bparam->false_positive_probability,
                    "False positive probability")
        ->capture_default_str();
    bfg->add_option("--random-seed", bparam->random_seed, "Random seed")
        ->capture_default_str();

    app.callback([=]() {
        // cout << "Input: " << *input << endl;
        // cout << "Output: " << *output << endl;

        bparam->compute_optimal_parameters();
        bloom_filter bfil(*bparam);

        dbg_map_t dbg = construct_dbg_v2(*input, *threads, *lcp_level, bfil);

        FILE *out = fopen(output->c_str(), "w");
        print_dbg(out, dbg);
        fclose(out);
    });

    return 0;
}

template <class T>
struct stat_accu {
    T count;
    T mean;
    T dS;
    T min;
    T max;

    stat_accu() : count(0), mean(0), dS(0), min(0), max(0) {}

    void update(T val) {
        count++;
        auto mean_dif = (val - mean) / count;
        auto new_mean = mean + mean_dif;

        auto dS_dif = (val - new_mean) * (val - mean);
        auto new_dS = dS + dS_dif;

        mean = new_mean;
        dS = new_dS;
    }

    constexpr T variance() const { return dS / count; }

    constexpr T stddev() const { return sqrt(variance()); }
};

void stats_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("bam-stats", "Print stats on input bam");
    auto &app = *app_;
    std::shared_ptr<string> input{new string};
    std::shared_ptr<uint64_t> update_frequency{new uint64_t{100}};
    std::shared_ptr<int> lcmer_level{new int{5}};
    app.add_option("input", *input, "Input bam file")->required();
    app.add_option("--update-frequency", *update_frequency, "Update frequency")
        ->capture_default_str();
    app.add_option("-l,--level", *lcmer_level, "Depth of lcmers to be computed")
        ->capture_default_str();

    app.callback([=]() {
        BAMReader reader(*input);
        bam1_t *bamdata = reader.init_data();
        // uint8_t *data;
        // uint64_t read_count;
        uint64_t base_count = 0;

        gtl::flat_hash_map<string, stat_accu<double>> stats;
        uint64_t read_count = 0;
        printf("\n");
        histo hist_len{"Read lengths", 100, 50, 20};
        histo hist_cores{"Core counts", 10, 10, 5};
        histo hist_core_lens{"Core lengths", 200, 100, 5};
        while (reader.read(bamdata)) {
            lcp_read lcps = parse_cores_from_bam_record(
                bam_get_seq(bamdata), bamdata->core.l_qseq, *lcmer_level);
            stats["length"].update(bamdata->core.l_qseq);
            stats["core-count"].update(lcps.cores.size());
            for (const auto &len : lcps.pos) {
                hist_core_lens.push(len.second - len.first);
                stats["core-length"].update(len.second - len.first);
            }
            base_count += bamdata->core.l_qseq;
            hist_len.push(bamdata->core.l_qseq);

            hist_cores.push(lcps.cores.size());

            if ((read_count % *update_frequency) == 0) {
                printf(
                    "\rcount: %llu\tbase-count: %llu\t length(µ,σ): %.2lf "
                    "%.2lf\tcores(µ,σ): %.2lf %.2lf\t core-lengths(µ,σ): %.2lf "
                    "%.2lf    "
                    "                       \n",
                    read_count, base_count, stats["length"].mean,
                    stats["length"].stddev(), stats["core-count"].mean,
                    stats["core-count"].stddev(), stats["core-length"].mean,
                    stats["core-length"].stddev());
                fflush(stdout);
                int reroll_count = 1;
                reroll_count += hist_len.print(stdout);
                reroll_count += hist_cores.print(stdout);
                reroll_count += hist_core_lens.print(stdout);
                hist_len.reroll(stdout, reroll_count);
            }
            ++read_count;
        }
        printf(
            "\rcount: %llu\tbase-count: %llu\t length(µ,σ): %.2lf "
            "%.2lf\tcores(µ,σ): %.2lf %.2lf\t core-lengths(µ,σ): %.2lf "
            "%.2lf    "
            "                       \n",
            read_count, base_count, stats["length"].mean,
            stats["length"].stddev(), stats["core-count"].mean,
            stats["core-count"].stddev(), stats["core-length"].mean,
            stats["core-length"].stddev());
        fflush(stdout);
        int reroll_count = 1;
        hist_len.print(stdout);
        hist_cores.print(stdout);
        hist_core_lens.print(stdout);
        return 0;
    });
}

void test_htslib_file_movements(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("test-hts",
                                       "Testing HTSLib Multi File movements");
    auto &app = *app_;
    struct parameters {
        string input_file;
        int threads;
        int block_size;
    };
    // Allocate parameters as heap and set default values
    std::shared_ptr<parameters> params{new parameters{"", 8, 4}};
    app.add_option("input", params->input_file, "Input file")->required();
    app.add_option("-t,--threads", params->threads, "Number of threads")
        ->capture_default_str();
    app.add_option("-b,--block-size", params->block_size, "Block size")
        ->capture_default_str();

    app.callback([=]() {

    });
}
int main(int argc, char **argv) {
    CLI::App app("lasm");
    app.require_subcommand(1);
    app.set_help_all_flag("--help-all", "Expand all help");
    asm_program(app);
    stats_program(app);
    test_htslib_file_movements(app);
    CLI11_PARSE(app, argc, argv);
    return 0;
}