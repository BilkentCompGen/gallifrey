
#include <concurrentqueue.h>
#include <stdio.h>

#include <chrono>
#include <gtl/phmap.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <string>
#include <iostream>
#include <mutex>



#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include "asm.h"
#include "bgzf.h"
#include "config.h"
#include "bloomfilter.h"
#include "lcmer.h"



using u64 = uint64_t;

using std::string;
enum class input_type { BAM, FASTX, FASTXGZ, FASTXBZ2, UNKNOWN };

#define NPH_BLOCKS 12
struct dbg_index {
    using mutex_type = std::mutex;
    struct bidbg_val {
        gtl::flat_hash_set<u64> prev;
        int cc;
        char *seq = NULL;
    };
    using dbg_map_t = gtl::parallel_flat_hash_map<
        u64, bidbg_val, std::hash<u64>, gtl::priv::hash_default_eq<u64>,
        std::allocator<std::pair<const u64, bidbg_val>>, NPH_BLOCKS,
        mutex_type>;
    using lcmer_index_t = gtl::parallel_flat_hash_map<
        lcmer<LCMER_L>, u64, std::hash<lcmer<LCMER_L>>,
        gtl::priv::hash_default_eq<lcmer<LCMER_L>>,
        std::allocator<std::pair<const lcmer<LCMER_L>, u64>>, NPH_BLOCKS,
        mutex_type>;

    dbg_map_t dbg;
    lcmer_index_t lci;
    bloom_interface<bloom_ext::RegisterBlockedBloomFilter<1>> bf;

    dbg_index(bloom_param bparam) : bf{bparam.false_positive_rate, bparam.projected_count, bparam.seed} {}

    u64 get(const lcmer<LCMER_L> &l) const {
        auto it = lci.find(l);
        if (it != lci.end()) {
            return it->second;
        }
        return 0;
    }

    u64 emplace(const lcmer<LCMER_L> &l, u64 &idx, int tid) {
        auto new_hash = lci.hash(l);
        if (!bf.contains(new_hash)) {
            bf.insert(new_hash);
            return 0;
        }
        auto vv = lci.try_emplace_with_hash(new_hash, l, idx | tid);
        if (vv.second) {  // Means new item on the table
            idx += 128;
        }
        return vv.first->second;
    }
};

struct assembler {
    string input;
    string output;
    int threads;
    int lcp_level;
    bloom_param bparam;
    assembler(const string &input, const string &output, int threads,
              int lcp_level, bloom_param bparam)
        : input{input},
          output{output},
          threads{threads},
          lcp_level{lcp_level},
          bparam{bparam} {}
  static auto process_cores(const lcp_read &lr, dbg_index &dbgi, int thread_id,
                   uint64_t &lcmer_index) {
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
        if (ll[i] == 0) {
            continue;
        }
        if (ll[i - 1] == 0) {
            if (ll[i + 1] == 0) {
                continue;
            } else {
                dbgi.dbg.lazy_emplace_l(
                    ll[i], [&](dbg_index::dbg_map_t::value_type &val) { val.second.cc++; },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{{}, 1, NULL});
                    });
            }
        } else {
        if (ll[i + 1] == 0) {
                dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](dbg_index::dbg_map_t::value_type &val) {
                        val.second.prev.emplace(ll[i - 1]);
                        val.second.cc++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{{ll[i - 1]}, 1, NULL});
                    });
            } else {
                dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&]( dbg_index::dbg_map_t::value_type &val) {
                        val.second.prev.emplace(ll[i - 1]);
                        val.second.cc++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{{ll[i - 1]}, 1, NULL});
                    });
            }
        }
    }
}



    static std::pair<int, int> get_thread_dist(int threads) {
        if (threads == 0) {
            return {0, 0};
        }
        if (threads < 4) {
            return {1, threads - 1};
        }
        return {threads / 2, threads - threads / 2};
    }
    using cqt = moodycamel::ConcurrentQueue<Chunk>;
    using dqt = moodycamel::ConcurrentQueue<bam1_t *>;

    static auto process_data_from_chunks(
        int thread_id, int lcp_level,
        string filename, cqt &chunk_queue,
        dbg_index &dbgi, std::atomic<bool> &done, int &processed_read_count,
        int &all_reads) {
        uint64_t lcmer_index = 0;
        

    }
   static auto process_data(
        int thread_id, int lcp_level,
        dqt &data_queue,
        dqt &recycle_queue, dbg_index &dbgi, std::atomic<bool> &done, int &processed_read_count,
        int &all_reads) {
        uint64_t lcmer_index = 0;
        bam1_t *data_p = NULL;
        lps str = {0};
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
                bam_get_seq(data_p), data_p->core.l_qseq, lcp_level, &str);
            if (lr.cores.size() < LCMER_L + 2) {
                // printf("Skipping read %lld, cores.size()=%zu\n", here,
                // lr.cores.size());
                goto pd_cleanup;
            }
            ++processed_read_count;
            process_cores(lr, dbgi, thread_id, lcmer_index);
        pd_cleanup:
            if (!recycle_queue.try_dequeue(data_p)) {
                bam_destroy1(data_p);
                data_p = NULL;
            }
        }
        free_lps(&str);
    }

    auto generate_chunks(
        std::string input, cqt &chunk_queue) {
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
        // Use htslib's thread pool for fast background decompression
        htsThreadPool pool = {NULL, 0};
        if (threads > 1) {
            pool.pool = hts_tpool_init(threads);
            // hts_set_thread_pool(fp, &pool);
            hts_set_cache_size(fp, 1024 * 1024 * 1024);
            hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool);
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
    auto construct_dbg() {
        dqt data_queue{static_cast<size_t>(threads) * 1024, 1, 1};
        dqt recycle_queue{static_cast<size_t>(threads) * 1024, 1, 1};
        std::atomic<bool> done{false};

        dbg_index dbgi{bparam};

        gtl::vector<std::thread> workers;
        auto [data_thread_count, worker_thread_count] =
            get_thread_dist(threads);
        gtl::vector<int> valid_read_counts(worker_thread_count, 0);
        gtl::vector<int> all_read_counts(worker_thread_count, 0);

        for (int tid = 1; tid <= worker_thread_count; ++tid) {
           workers.emplace_back(process_data, tid, lcp_level,
                             std::ref(data_queue), std::ref(recycle_queue),
                             std::ref(dbgi), std::ref(done),
                             std::ref(valid_read_counts[tid-1]), std::ref(all_read_counts[tid-1]));
        }

        generate_data(input, std::ref(data_queue), std::ref(recycle_queue),
                      data_thread_count);

        done.store(true);
        for (auto &worker : workers) {
            worker.join();
        }

        return  dbgi;
    }

    auto convert2basespace() -> int { return 0; }
};

//TODO do not use extensions to identify input type.
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

void print_raw_dbg(FILE *out, const dbg_index &dbgi) {
    for (const auto &[lval, edges] : dbgi.dbg) {
        fprintf(out, "S\tC%llu\t%s\tLN:i:%d\tRC:i:%d\n", lval, "*", 100, edges.cc);
    }

    for (const auto &[lval, edges] : dbgi.dbg) {
        for (const auto &p : edges.prev) {
            fprintf(out, "L\tC%llu\t+\tC%llu\t+\t%dM\n", p, lval, 1);
        }
    }

}

int run_assembly(const std::string &input, const std::string &output,
                 int threads, int lcp_level, bloom_param bpd) {
    INIT_TIMERS(HIFI_ASSEMBLY);

    auto inty = identify_input(input);

    if (inty == input_type::UNKNOWN) {
        fprintf(stderr, "Unknown input type %s\n", input.c_str());
        return -1;
    }
    assembler as(input, output, threads, lcp_level, bpd);
    auto dbg = as.construct_dbg();

    TIME_CHECKPOINT(HIFI_ASSEMBLY, stderr, "DBG construction finished");

    FILE *out = fopen(output.c_str(), "w");
    print_raw_dbg(out, dbg);

    fclose(out);

    return 0;
}
