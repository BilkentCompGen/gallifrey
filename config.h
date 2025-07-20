#ifndef _LASM_CONFIG_H_
#define _LASM_CONFIG_H_

#ifndef LCMER_L
#define LCMER_L 8
#endif

#ifndef SYNCMER_K
#define SYNCMER_K 8
#endif

#ifndef SYNCMER_S
#define SYNCMER_S 1
#endif

#ifndef SYNCMER_T
#define SYNCMER_T 20
#endif

#ifndef SYNCMER_RULE
#define SYNCMER_RULE DensitySyncmerRule
#endif

#ifndef KMER_SPACE_MASK
#define KMER_SPACE_MASK 0xDDDDDDDDDDDDDDDD
#endif 

#ifndef BAMREAD_CHUNK_SIZE
#define BAMREAD_CHUNK_SIZE 256ULL
#endif

#define INIT_TIMERS(name) \
    std::chrono::high_resolution_clock::time_point name##_chrono_start = std::chrono::high_resolution_clock::now(); \
    std::chrono::high_resolution_clock::time_point name##_last_checkpoint = name##_chrono_start; \
    clock_t name##_clock_start = std::clock(); \
    clock_t name##_last_clock = name##_clock_start;


#define TIME_CHECKPOINT(name, stream, fmt, ...) \
    do { \
        auto name##_now = std::chrono::high_resolution_clock::now(); \
        clock_t name##_clock_now = std::clock(); \
        \
        long long name##_total_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(name##_now - name##_chrono_start).count(); \
        long long name##_interval_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(name##_now - name##_last_checkpoint).count(); \
        \
        double name##_total_cpu = 1000.0 * (name##_clock_now - name##_clock_start) / CLOCKS_PER_SEC; \
        double name##_interval_cpu = 1000.0 * (name##_clock_now - name##_last_clock) / CLOCKS_PER_SEC; \
        \
        double name##_total_util = (name##_total_elapsed > 0) ? (name##_total_cpu / name##_total_elapsed) * 100.0 : 0.0; \
        double name##_interval_util = (name##_interval_elapsed > 0) ? (name##_interval_cpu / name##_interval_elapsed) * 100.0 : 0.0; \
        \
        fprintf(stream, "[Checkpoint: %s] ", #name); \
        fprintf(stream, fmt, ##__VA_ARGS__); \
        fprintf(stream, "\n"); \
        fprintf(stream, "  Total:   Wall = %lld ms, CPU = %.2f ms, Utilization = %.2f%%\n", name##_total_elapsed, name##_total_cpu, name##_total_util); \
        fprintf(stream, "  Current: Wall = %lld ms, CPU = %.2f ms, Utilization = %.2f%%\n\n", name##_interval_elapsed, name##_interval_cpu, name##_interval_util); \
        \
        name##_last_checkpoint = name##_now; \
        name##_last_clock = name##_clock_now; \
    } while (0)


#endif