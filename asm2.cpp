#include <blockingconcurrentqueue.h>
#include <fcntl.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <stdio.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gtl/phmap.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include <gtl/vector.hpp>
#include <iostream>
#include <string>
#include <utility>

#include "asm.h"
#include "bgzf.h"
#include "config.h"
#include "ext/doctest.h"
#include "graph.h"
#include "lcmer.h"
#include "lcptools/lcptools_ho.h"
#include "syncmer.h"

using std::string;
enum class input_type { SAM, BAM, FASTX, FASTXGZ, FASTXBZ2, UNKNOWN };

template <class T>
struct mptr {
    static constexpr size_t vmask = 3;
    static constexpr size_t pmask = ~vmask;
    T *ptr;
    ~mptr() {
        delete unwrap();
        ptr = 0;
    }
    mptr(T val) : ptr{new T{val}} {}
    mptr(T val, int message) : ptr{new T{val}} { write_message(message); }
    mptr(T *ptr) : ptr{ptr} {}
    mptr(T *ptr, int message) : ptr{ptr} { write_message(message); }
    mptr(mptr &&other) {
        ptr = other.ptr;
        other.ptr = 0;
    }
    mptr(const mptr &other) = delete;
    mptr &operator=(const mptr &other) = delete;
    mptr &&operator=(mptr &&other) {
        this->ptr = other.ptr;
        other.ptr = 0;
    }
    constexpr const T *unwrap() const { return (T *)((size_t)(ptr)&pmask); }
    constexpr T *unwrap() { return (T *)((size_t)(ptr)&pmask); }
    constexpr T operator*() const {
        if (ptr) return *unwrap();
        throw 42;
    }
    constexpr int get_message() const { return (size_t)(ptr)&vmask; }
    void write_message(int message) {
        ptr = (T *)((size_t)(ptr) | (message & vmask));
    }
};

struct assembler {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSyncmerRule<uint32_t, SYNCMER_HASHER>;
    string input;
    string output;
    decltype(asm_config::threads) threads;
    int lcp_level;
    double density;
    decltype(asm_config::bloom) bloom;
    uint64_t seed;
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
                        val.second.core_count++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{1});
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
                        ctr(ll[i], dbg_index::bidbg_val{0U, ll[i - 1]});
                    });
            } else {
                return dbgi.dbg.lazy_emplace_l(
                    ll[i],
                    [&](dbg_index::dbg_map_t::value_type &val) {
                        val.second.emplace_prev(ll[i - 1]);
                        val.second.core_count++;
                    },
                    [&](const dbg_index::dbg_map_t::constructor &ctr) {
                        ctr(ll[i], dbg_index::bidbg_val{1U, ll[i - 1]});
                    });
            }
        }
    }
    static gtl::vector<uint64_t> find_lcmers(const lcp_read &lr,
                                             const dbg_index &dbgi) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        gtl::vector<uint64_t> ll(1, 0);
        ll.reserve(lr.cores.size() - LCMER_L + 3);
        dbg_index::lci_type lval = dbgi.get(l);
        ll.push_back(lval.id);
        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            lval = dbgi.get(l);
            ll.push_back(lval.id);
        }
        ll.push_back(0);
        return ll;
    }
    static gtl::vector<std::tuple<u64, u64, u64>> find_lcmer_ranges(
        const lcp_read &lr, const dbg_index &dbgi) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64>> llr;

        dbg_index::lci_type lval = dbgi.get(l);
        llr.emplace_back(lval.id, 0, 0);
        for (size_t i = 2; i < lr.cores.size() - LCMER_L - 1; ++i) {
            uint32_t in_core = lr.cores[i + LCMER_L];
            l.update_lcmer(in_core);
            lval = dbgi.get(l);
            if (lval.id == std::get<0>(llr.back())) {  // Extend the previous
                std::get<2>(llr.back()) = i;
            } else if (lval.id != 0) {
                llr.emplace_back(lval.id, i, i);
            }
        }
        return llr;
    }
    static void process_cores(const lcp_read &lr, dbg_index &dbgi,
                              int thread_id, uint64_t &lcmer_index) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;
        gtl::vector<std::tuple<u64, u64, u64>> llr;

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
                // Lets try atomic
                dbgi.dbg.modify_if(
                    ll[i], [&lr, i](dbg_index::dbg_map_t::value_type &val) {
                        val.second.seq.emplace_back(
                            lr.data, lr.pos[i + LCMER_L - 1].second,
                            lr.pos[i + LCMER_L].second, val.first);
                    });
                // dbgi.dbg[ll[i]].seq.emplace_back(lr.data,
                //                                  lr.pos[i + LCMER_L -
                //                                  1].second, lr.pos[i +
                //                                  LCMER_L].second);
            }
        }
    }
    static void process_cores_for_basespace(
        const lcp_read &lr, dbg_index &dbgi, [[maybe_unused]] int thread_id,
        [[maybe_unused]] uint64_t &lcmer_index) {
        lcmer<LCMER_L> l{lr.cores.begin() + 1};  // Skip first core;

        uint64_t val = dbgi.get(l).id;

        for (size_t i = LCMER_L + 1; i < lr.cores.size() - 1; ++i) {
            uint64_t from = lr.pos[i - 1].second;
            uint64_t to = lr.pos[i].second;

            dbgi.dbg.modify_if(val, [&](dbg_index::dbg_map_t::value_type &val) {
                if (val.second.seq.size() == 0) {
                    if (to - from > 0) {
                        val.second.seq.emplace_back(lr.data, from, to,
                                                    val.first);
                        // assert((to & 3) + 4 < lr.data_length);
                    }
                }
            });
            uint32_t in_core = lr.cores[i];
            l.update_lcmer(in_core);
            val = dbgi.get(l).id;
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
    static auto process_data(int thread_id, int lcp_level, dqt &data_queue,
                             dqt &recycle_queue, dbg_index &dbgi,
                             std::atomic<bool> &done, int &processed_read_count,
                             int &all_reads,
                             decltype(assembler::process_cores) Process_Func,
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
            parse_cores_from_bam_record(bam_get_seq(data_p),
                                        data_p->core.l_qseq, lcp_level, &str,
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

    [[nodiscard]] auto generate_data(const std::string &input, dqt &data_queue,
                                     dqt &recycle_queue, int threads) {
        samFile *fp = sam_open(input.c_str(), "r");

        if (!fp) {
            std::cerr << "Error: Could not open " << input << std::endl;
            return false;
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
        return true;
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

        if (!generate_data(input, std::ref(data_queue), std::ref(recycle_queue),
                           threads.io)) {
            exit(-1);
        }

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
    template <class map, class F>
    static void for_each_map(const map &mp, const F &fCallback, int threads) {
        std::vector<std::thread> workers;
        for (int tid = 0; tid < threads; ++tid) {
            workers.emplace_back(
                [&](const map &mp, int tid) {
                    for (size_t i = tid; i < mp.subcnt(); i += threads) {
                        mp.with_submap(i, [&](const map::EmbeddedSet &set) {
                            for (const auto &p : set) {
                                fCallback(p);
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

    size_t count_seqv_len(const gtl::vector<SequenceContainerType> &seqv) {
        size_t sum = 0;
        for (const SequenceContainerType &s : seqv) {
            sum += s.len;
        }
        return sum;
    }
    // I need to recompute the LCP, it is likely more expensive to find relavant
    // nodes.
    void iterative_re_dbg(dbg_index &dbgi) {
        gtl::vector<std::pair<u64, size_t>> node_id_len;
        for (const auto &[k, v] : dbgi.dbg) {
            node_id_len.emplace_back(k, v.seq.size());
        }
        std::sort(
            node_id_len.begin(), node_id_len.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });
        // Check this again TODO.
        dbg_index new_index(dbgi.bf.implementation.m, dbgi.bf.implementation.n,
                            dbgi.bf.random_seed);

        auto it = node_id_len.begin();
        assert(it != node_id_len.end());
        lps str;
        auto calculate_reverse_utig = [](auto node) { return lps{}; };
        auto move_forward_unitig = [&](u64 id) {
            auto node_it = dbgi.dbg.find(id);
            if (node_it != dbgi.dbg.end()) {
                int len = count_seqv_len(node_it->second.seq);
                prepare_lps(&str, len);
                init_lps_crumbled(&str, node_it->second.seq, syncmer_functor);
                lps_deepen(&str, lcp_level);
            }
        };
        auto rev_contig_present_in_the_index =
            [&](const dbg_index::dbg_map_t::value_type &node, auto ctg) {
                return true;
            };
        move_forward_unitig(node_id_len.front().first);
        ++it;
        while (it != node_id_len.end()) {
            auto rev_ctg = calculate_reverse_utig(*it);
            if (rev_contig_present_in_the_index(*it, rev_ctg)) {
            } else {
                move_forward_unitig(it->first);
            }
            ++it;
        }
    }
    struct rev_match {
        u64 node;
        u64 from;
        u64 to;
        u64 size;
    };
    void fold_dbg_graph(dbg_index &dbgi) {
        fprintf(stderr, "dbg_size: %zu\n", dbgi.dbg.size());
        fprintf(stderr, "lci_size: %zu\n", dbgi.lci.size());

        std::atomic<int> all_n{0};
        std::atomic<int> hits{0};
        std::atomic<int> bad{0};
        gtl::parallel_flat_hash_map_d<u64, gtl::vector<rev_match>, HT_SHARD_BC,
                                      std::mutex>
            matched;
        for_each_map(
            dbgi.dbg,
            [&](const auto &node) {
                all_n.fetch_add(1);
                if (dbgi.active_set.find(node.first) == dbgi.active_set.end()) {
                    bad.fetch_add(1);
                    return;
                }

                lps str = {0, 0, 0, 0};
                lcp_read lr;
                // dbgi.dbg.for_each([&](auto &node) {
                if (node.second.redirected()) {
                    bad.fetch_add(1);
                    return;
                }
                // good_n.fetch_add(1);
                lr.clear();

                gtl::vector<SequenceContainerType> rev_seqs;
                for (auto it = node.second.seq.rbegin();
                     it != node.second.seq.rend(); ++it) {
                    rev_seqs.emplace_back(it->reverse_complement());
                }
                init_lps_crumbled(&str, rev_seqs, syncmer_functor);
                lps_deepen(&str, lcp_level);

                for (int i = 0; i < str.size; ++i) {
                    lr.cores.emplace_back(str.cores[i].label);
                    lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
                }
                if (lr.cores.size() < LCMER_L + 1) {
                    bad.fetch_add(1);
                    // good_n.fetch_sub(1);
                    free_lps(&str);
                    return;
                }
                auto found_ranges = find_lcmer_ranges(lr, dbgi);
                for (auto [id, from, to] : found_ranges) {
                    if (id != 0) {
                        hits.fetch_add(1);
                        matched[node.first].emplace_back(id, from, to,
                                                         rev_seqs.size());
                        // printf("%llu to %llu\t%llu-%llu seqs: %zu %% %lf\n",
                        //        node.first, id, from, to, lr.cores.size(),
                        //        (double)(to - from + 1) /
                        //            (lr.cores.size() - LCMER_L));
                    }
                }
                free_lps(&str);
            },
            threads.work);

        // for (const auto &[k, v] : dbgi.dbg) {
        //     if (!v.redirected()) {
        //         printf("+\t%llu\t%zu\t%zu\t%u\n", k, v.seq.size(),
        //                v.prev.size(), v.forward_edge_count);
        //     }
        // }
        gtl::flat_hash_map<u64, u64> reciprocal;
        for (const auto &[k, v] : matched) {
            for (const auto &match : v) {
                // if (k < match.node) {  // smaller node is the key
                auto it = matched.find(
                    match.node);  // and matched node also has overlaps
                if (it != matched.end()) {
                    auto ff = std::find_if(
                        it->second.begin(), it->second.end(),
                        [k](const auto &val) { return val.node == k; });
                    if (ff != it->second.end()) {  // and has the reverse match
                        reciprocal.emplace(k, match.node);

                        // printf("%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n", k,
                        //        match.node, match.from, match.to,
                        //        ff->from, ff->to);
                        printf(
                            "C%llu\t%llu\t%llu\t%llu\t+\tC%llu\t%llu\t%"
                            "llu\t%llu\t%llu\t%llu\t%d\n",
                            k, match.size + LCMER_L, match.from, match.to,
                            match.node, ff->size + LCMER_L, ff->from, ff->to,
                            ff->to - ff->from, ff->to - ff->from, 250);
                    } else { // Handle after reciprocals are good
                        // printf(
                        //     "C%llu\t%llu\t%llu\t%llu\t+\tC%llu\t%d\t%d\t%d\t%"
                        //     "d\t%d\t%d\n",
                        //     k, match.size + LCMER_L, match.from, match.to,
                        //     match.node, -1, -1, -1, -1, -1, 250);
                    }
                    // }
                }

                // printf("%llu\t%llu\t%llu\t%llu\n", k, match.node, match.from,
                //        match.to);
            }
        }
        fprintf(stderr, "all:%d\ngood: %d nodes\nhits: %d\n", all_n.load(),
                all_n.load() - bad.load(), hits.load());
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

int run_assembly(const asm_config &cfg) {
    if (!cfg.silent) cfg.print_cfg(stderr);
    // free(malloc(1024ULL * 1024 * 1024 * 20)); INTIMIDATE the OS for the
    // ram

    INIT_TIMERS(HIFI_ASSEMBLY);

    // auto inty = identify_input(cfg.input);

    // if (inty == input_type::UNKNOWN) {
    //     fprintf(stderr, "Unknown input type %s\n", cfg.input.c_str());
    //     return -1;
    // }

    assembler as(cfg);
    dbg_index dbg = as.construct_dbg();

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG construction");

    dbg.count_forward_edges();
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr,
                       "Count forward edges");

    if (cfg.print_precleanup_graph) {
        FILE *out = fopen((cfg.output + ".raw").c_str(), "w");
        print_dbg_seqs(out, dbg);
        fclose(out);

        IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG printing");
    }

    GraphCleaning::run(dbg, cfg);

    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Graph cleaning");

    as.fold_dbg_graph(dbg);
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "Fold graph");

    FILE *out = fopen(cfg.output.c_str(), "w");

    print_dbg_seqs(out, dbg);
    fclose(out);
    IF_TIME_CHECKPOINT(!cfg.silent, HIFI_ASSEMBLY, stderr, "DBG printing");

    return 0;
}

void tests_mock_fill_seq(dbg_index &dbgi) {
    for (auto &node : dbgi.dbg) {
        node.second.seq.emplace_back();
    }
}

// TEST_CASE("Graph simplification simple") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     int MOCK_CORE_COUNT = 40;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     lcmer<LCMER_L> l(mock_cores.begin());
//     u64 idx = 0;

//     gtl::vector<u64> ll(1, 0);
//     ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L + 1; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);

//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     SUBCASE("Check ll vec") {
//         for (size_t i = 1; i < ll.size() - 1; ++i) {
//             CHECK(ll[i] == (i - 1) * 128);
//         }
//     }
//     SUBCASE("Pre simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//             CHECK(node.second.core_count < 2);
//             CHECK(node.second.prev.size() < 2);
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.lci.size() == MOCK_CORE_COUNT - 2 - LCMER_L + 1);
//     }

//     tests_mock_fill_seq(dbgi);
//     for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
//         n.second.core_count *= 10;
//     }
//     dbgi.count_forward_edges();

//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});
//     SUBCASE("Post simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//         CHECK(dbgi.dbg.begin()->second.seq.size() == 30);
//     }

//     // // Write Dot file
//     // FILE *dot_file = fopen("graph.dot", "w");

//     // if (dot_file) {
//     //     // print_raw_dbg(dot_file, dbgi);
//     //     fprintf(dot_file, "digraph G {\n");
//     //     for (auto &node : dbgi.dbg) {
//     //         fprintf(dot_file, "\t%" PRIu64 " [label=\"%" PRIu64
//     "\"];\n",
//     //         node.first,
//     //                 node.first);
//     //         for (auto &prev : node.second.prev) {
//     //             fprintf(dot_file, "\t%" PRIu64 " -> %" PRIu64 ";\n",
//     prev,
//     //             node.first);
//     //         }
//     //     }
//     //     fprintf(dot_file, "}\n");
//     //     fclose(dot_file);
//     // }
// }

// TEST_CASE("Graph simplification Bubble") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     int MOCK_CORE_COUNT = 40;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     lcmer<LCMER_L> l(mock_cores.begin());
//     u64 idx = 0;

//     gtl::vector<u64> ll(1, 0);
//     ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     ll.clear();
//     mock_cores[MOCK_CORE_COUNT / 2] = 123123123;
//     l = lcmer<LCMER_L>(mock_cores.begin());

//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     dbgi.count_forward_edges();
//     for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
//         n.second.core_count *= 10;
//     }
//     tests_mock_fill_seq(dbgi);
//     print_dbg_dot(dbgi, 0);

//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});
//     print_dbg_dot(dbgi, 1);

//     SUBCASE("Post simplify dbg iter 2") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//     }

//     print_dbg_dot(dbgi, 2);
// }

// TEST_CASE("Graph simplification Tips (Outgoing)") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     int MOCK_CORE_COUNT = 30;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     lcmer<LCMER_L> l(mock_cores.begin());
//     u64 idx = 0;

//     gtl::vector<u64> ll(1, 0);
//     ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }
//     for (auto &n : dbgi.dbg) {  // Main Path has 10X depth
//         n.second.core_count *= 10;
//     }
//     ll.clear();
//     mock_cores[12] = 123123123;
//     l = lcmer<LCMER_L>(mock_cores.begin());

//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() / 2 - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     dbgi.count_forward_edges();

//     tests_mock_fill_seq(dbgi);
//     print_dbg_dot(dbgi, 3);

//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});
//     SUBCASE("Post simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//     }
//     print_dbg_dot(dbgi, 4);
// }
// TEST_CASE("Graph simplification Tips (Incoming)") {
//     dbg_index dbgi{0.95, 10, 0};
//     gtl::vector<u32> mock_cores;
//     lcmer<LCMER_L> l;
//     gtl::vector<u64> ll;

//     u64 idx = 0;
//     int MOCK_CORE_COUNT = 30;
//     for (int i = 0; i < MOCK_CORE_COUNT; i++) {
//         mock_cores.push_back(MurmurHash3_32(&i, sizeof(int), 156));
//     }
//     for (int DEPTH = 0; DEPTH < 10; ++DEPTH) {
//         l = lcmer<LCMER_L>(mock_cores.begin());

//         ll.push_back(0);

//         ll.reserve(MOCK_CORE_COUNT - LCMER_L + 3);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//         for (int i = LCMER_L; i < MOCK_CORE_COUNT - 1; ++i) {
//             uint32_t in_core = mock_cores[i];
//             l.update_lcmer(in_core);
//             ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//         }
//         ll.push_back(0);
//         for (size_t i = 1; i <= ll.size() - 1; i++) {
//             assembler::insert_lcmer_helper(ll, i, dbgi);
//         }

//         ll.clear();
//     }
//     ll.push_back(0);
//     mock_cores[8] = 123123123;
//     l = lcmer<LCMER_L>(mock_cores.begin() + 7);

//     ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     for (int i = LCMER_L + 7; i < MOCK_CORE_COUNT - 1; ++i) {
//         uint32_t in_core = mock_cores[i];
//         l.update_lcmer(in_core);
//         ll.push_back(dbgi.emplace_nobloomfilter(l, idx, 0));
//     }
//     ll.push_back(0);
//     for (size_t i = 1; i <= ll.size() / 2 - 1; i++) {
//         assembler::insert_lcmer_helper(ll, i, dbgi);
//     }

//     dbgi.count_forward_edges();

//     tests_mock_fill_seq(dbgi);
//     print_dbg_dot(dbgi, 6);
//     GraphCleaning::run(dbgi, {.graph_cleaning = {1, 10, 10, 2}});

//     SUBCASE("Post simplify dbg") {
//         for (auto &node : dbgi.dbg) {
//             INFO(node.first, "\t.\t.\t", node.second.seq.size());
//         }
//         CAPTURE(dbgi.dbg.size());
//         CHECK(dbgi.dbg.size() == 1);
//     }
//     print_dbg_dot(dbgi, 7);
// }
