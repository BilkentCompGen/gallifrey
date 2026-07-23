#ifndef _BAMDOCTOR_H_
#define _BAMDOCTOR_H_

#include <string>
#include <vector>
#include <functional>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <stdio.h>
#include "config.h"
#include <thread>
#include <atomic>
#include <concurrentqueue.h>
#include <gtl/vector.hpp>
#include "lcmer.h"

#include "syncmer.h"

enum class SamReadTestingMode {
    BASIC,
    POOL,
};

inline int run_sam_read_testing(const std::string &input, int threads, int lcp_level, SamReadTestingMode mode, const std::string &job_type, double density){
    INIT_TIMERS(BAM_READING);
    //Explain the program nicely formatted brief
    printf("Input: %s\n", input.c_str());
    printf("Threads: %d\n", threads);
    printf("LCP level: %d\n", lcp_level);
    printf("Mode: %s\n", mode == SamReadTestingMode::BASIC ? "BASIC" : "POOL");


    // Read bam file with htslib
    samFile *fp = sam_open(input.c_str(), "r");
    if(!fp){
        fprintf(stderr, "Error: Could not open %s\n", input.c_str());
        return 1;
    }
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if(!hdr){
        fprintf(stderr, "Error: Could not read header from file %s\n", input.c_str());
        sam_close(fp);
        return 1;
    }
    bam1_t *bamdata = bam_init1();
    if(!bamdata){
        fprintf(stderr, "Error: Could not initialize bamdata\n");
        sam_hdr_destroy(hdr);
        sam_close(fp);
        return 1;
    }
    htsThreadPool pool = {NULL, 0};
    if(mode == SamReadTestingMode::BASIC){
        hts_set_threads(fp, threads);
    }else if(mode == SamReadTestingMode::POOL){

        if(threads > 1){
            pool.pool = hts_tpool_init(threads);
            if(!pool.pool){
                fprintf(stderr, "Error: Could not initialize thread pool\n");
                bam_destroy1(bamdata);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                return 1;
            }
            if(hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool) < 0){
                fprintf(stderr, "Error: Could not set thread pool\n");
                hts_tpool_destroy(pool.pool);
                bam_destroy1(bamdata);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                return 1;
            }
        }
    }
    moodycamel::ConcurrentQueue<bam1_t *> data_queue {1024, 1, static_cast<size_t>(threads)};
    moodycamel::ConcurrentQueue<bam1_t *> recycle_queue {1024, 1, static_cast<size_t>(threads)};
    gtl::vector<std::thread> workers;
    std::atomic<bool> done{false};
    // One slot per worker; sized by the (CLI-clamped) thread count instead of a
    // fixed 1024 so a large -t can never index out of bounds.
    std::vector<int> counts(threads, 0);
    // Use the same rule the assembly pipeline uses (GenericSyncmerRule) so the
    // diagnostic counts measure what production computes (H8).
    syncmer_params sr_params;
    sr_params.k = 8;
    sr_params.density = density;
    syncmer::GenericSyncmerRule<uint32_t, syncmer::MurMur> sr{sr_params, 423};
    // Per-worker scratch (charsi/pos_scratch) is passed in as arguments rather
    // than captured, so the single shared `job` closure stays race-free.
    auto job = [&]() -> std::function<int(bam1_t *, lps *, charsi &, pos_scratch &)> {
        if (job_type == "read_count"){
            return [](bam1_t *, lps *, charsi &, pos_scratch &) -> int {
                return 1;
            };
        }
        else if (job_type == "core_count"){
            return [&, lcp_level](bam1_t *bamdata, lps *str, charsi &chs, pos_scratch &pos_buf) -> int {
                lcp_read lr;
                lr.clear();
                parse_cores_from_bam_record(bam_get_seq(bamdata), bamdata->core.l_qseq, lcp_level, str, sr, chs, lr, pos_buf);
                // if (lr.cores.size() < LCMER_L + 2) {
                //     return 0;
                // }

                return (int)lr.cores.size();
            };
        }
        return [](bam1_t *, lps *, charsi &, pos_scratch &) -> int {
            return 0;
        };
    }();

    for(int i = 0; i < threads; ++i){
        workers.emplace_back([&, i](){
            bam1_t *bamdata = NULL;
            lps str = {0, 0, 0, nullptr};
            charsi chs;
            pos_scratch pos_buf;
            while(!done.load() || data_queue.size_approx() > 0){
                if(!data_queue.try_dequeue(bamdata)){
                    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
                    continue;
                }
                counts[i] += job(bamdata, &str, chs, pos_buf);
                // Hand the buffer back for the producer to refill. Only the
                // producer reads into it, and only after it has been recycled,
                // so no worker is still touching it.
                if(!recycle_queue.try_enqueue(bamdata)){
                    bam_destroy1(bamdata);
                }
            }
            if(str.cores != nullptr){
                free_lps(&str);
            }
        });
    }

    // Fill `bamdata`, hand it off to a worker, then obtain the *next* buffer
    // from the recycle queue (or allocate one). This is the correct direction:
    // the producer never reuses a record a worker may still be reading (H3).
    while(sam_read1(fp, hdr, bamdata) >= 0){
        while(!data_queue.try_enqueue(bamdata)){
            std::this_thread::sleep_for(std::chrono::nanoseconds(100));
        }
        if(!recycle_queue.try_dequeue(bamdata)){
            bamdata = bam_init1();
        }
    }
    done.store(true);
    for(auto &worker : workers){
        worker.join();
    }
    size_t count = 0;
    for(int i = 0; i < threads; ++i){
        count += counts[i];
    }
    // Reclaim every buffer. After the joins no worker is running: data_queue is
    // drained by the workers (loop condition keeps them going until empty), but
    // destroy any straggler defensively; recycle_queue holds the buffers handed
    // back; `bamdata` is the spare read buffer left over at EOF.
    bam1_t *item;
    while(data_queue.try_dequeue(item)){
        bam_destroy1(item);
    }
    while(recycle_queue.try_dequeue(item)){
        bam_destroy1(item);
    }
    bam_destroy1(bamdata);
    TIME_CHECKPOINT(BAM_READING, stderr, "Read %zu reads", count);

    sam_hdr_destroy(hdr);
    sam_close(fp);
    if(pool.pool){
        hts_tpool_destroy(pool.pool);
    }

    return 0;
}


#endif
