#pragma once

#ifndef _ASM_H_
#define _ASM_H_

#include <inttypes.h>
#include <stdint.h>

#include <cstdio>
#include <expected>
#include <optional>
#include <string>
#include <vector>
// Syncmer selection parameters (Edgar 2021). A k-mer is selected by the
// position µ of its smallest s-mer: a closed syncmer has µ at the first or
// last s-mer, an open syncmer has µ at offset t.
struct syncmer_params {
    // Selection packs k-mers into 32-bit codes.
    static constexpr unsigned max_k = 16;
    // Rule used when neither s nor t is given: open syncmer, s = 4, t = 5.
    static constexpr unsigned long default_s = 4;
    static constexpr unsigned long default_t = 5;

    enum class typeenum { OpenSyncmer, ClosedSyncmer, Invalid };
    constexpr typeenum type() const {
        if (!smer.s.has_value()) {
            return typeenum::Invalid;
        }
        return smer.t.has_value() ? typeenum::OpenSyncmer
                                  : typeenum::ClosedSyncmer;
    }
    unsigned k = 11;
    // s-mer length s and (open syncmer) offset t.
    struct {
        std::optional<unsigned long> s;
        std::optional<unsigned long> t;
    } smer;

    // With no s and no t, select the default open syncmer. An s without a t
    // stays a closed syncmer.
    void apply_default_rule() {
        if (!smer.s.has_value() && !smer.t.has_value()) {
            smer.s = default_s;
            smer.t = default_t;
        }
    }

    // Validate the whole syncmer configuration once, before any syncmer
    // functor is constructed, so callers can surface a clean error.
    std::expected<void, std::string> validate() const {
        if (type() == typeenum::Invalid) {
            return std::unexpected(
                std::string("Provide --submer-s (optionally with --submer-t)"));
        }
        if (k > max_k) {
            return std::unexpected("--syncmer-k (" + std::to_string(k) +
                                   ") must be <= " + std::to_string(max_k));
        }
        if (*smer.s == 0 || *smer.s >= k) {
            return std::unexpected("--submer-s (" + std::to_string(*smer.s) +
                                   ") must be in [1, --syncmer-k " +
                                   std::to_string(k) + ")");
        }
        if (smer.t.has_value() && *smer.t >= k - *smer.s + 1) {
            return std::unexpected("--submer-t (" + std::to_string(*smer.t) +
                                   ") must be < k - s + 1 (" +
                                   std::to_string(k - *smer.s + 1) + ")");
        }
        return {};
    }
    void print_cfg(FILE *f) const {
        fprintf(f, "k: %u\n", k);
        if (!smer.s.has_value()) {
            fprintf(f, "submer-rule: \"none\"\n");
        } else if (smer.t.has_value()) {
            fprintf(f, "submer-rule: \"open-syncmer\"\n");
            fprintf(f, "syncmer-s: %lu\nsyncmer-t: %lu\n", *smer.s, *smer.t);
        } else {
            fprintf(f, "submer-rule: \"closed-syncmer\"\n");
            fprintf(f, "syncmer-s: %lu\n", *smer.s);
        }
    }
};
// Default --random-seed. The seed only drives the syncmer s-mer hash, which
// is 32-bit, so the seed is stored as uint32_t.
inline constexpr uint32_t kDefaultRandomSeed = 1453;
struct asm_config {
    std::vector<std::string> inputs;
    std::string output;
    int lcp_level = 3;
    // LCMER window length (cores) of the DBG keys; overlap detection queries
    // the same index, so it uses the same length.
    int lcmer_l = 6;
    struct {
        int io = 4;
        int work = 8;
    } threads;
    uint32_t seed = kDefaultRandomSeed;
    bool rle = true;
    bool parallel_file_reading = true;
    syncmer_params syncmer;
    bool print_precleanup_graph = false;
    bool print_prelayout_graph = false;
    // Drop crumb nodes (fewer than lcmer_l+1 cores) that continue into a
    // longer node instead of emitting them as C segments.
    bool drop_crumbs = true;
    // Also parse every read's reverse complement and insert its lcmer keys
    // into the Bloom filter only, so the other strand's reads are admitted at
    // their first sighting.
    bool prime_reverse = true;
    bool silent = false;
    struct {
        uint64_t projected_count = 100'000'000;
        double false_positive_rate = 0.001;
    } bloom;
    struct {
        int iterations = 6;
        int max_tip_length = 8 * 2;
        double min_bubble_length_difference = 2;
        // Inner cores of the largest bubble to pop; 0 disables popping.
        size_t max_bubble_length = 1024;
        // Bubbles within this many copies lie at one locus and are popped.
        int ploidy = 1;
    } graph_cleaning;

    void print_cfg(FILE *f) const {
        if (!inputs.empty()) {
            fprintf(f, "Inputs (%zu):\n", inputs.size());
            for (const auto &in : inputs) {
                fprintf(f, "\t%s\n", in.c_str());
            }
        } else {
            fprintf(f, "Inputs (0):\n");
        }
        fprintf(f, "Output: %s\n", output.c_str());
        fprintf(f, "LCP Level: %d\n", lcp_level);
        fprintf(f, "LCMER Length: %d\n", lcmer_l);
        fprintf(f, "Threads:\n\tIO: %d\n\tWork: %d\n", threads.io, threads.work);
        fprintf(f, "Seed: %" PRIu32 "\n", seed);
        fprintf(f, "RLE: %s\n", rle ? "true" : "false");
        fprintf(f, "File Reading Mode: %s\n",
                parallel_file_reading ? "parallel" : "sequential");
        syncmer.print_cfg(f);
        fprintf(f,
                "Bloom Filter Parameters:\n\tProjected Count: %" PRIu64
                "\n\tFalse Positive "
                "Rate: %.3f\n",
                bloom.projected_count, bloom.false_positive_rate);
        fprintf(f, "Drop Crumbs: %s\n", drop_crumbs ? "true" : "false");
        fprintf(f, "Prime Reverse: %s\n", prime_reverse ? "true" : "false");
        fprintf(f,
                "Graph Cleaning Parameters: \n\tIterations: %d\n\tMax Tip "
                "Length: %d\n\tMin Bubble Length Difference: %.2f\n\tMax "
                "Bubble Length: %zu\n\tPloidy: %d\n",
                graph_cleaning.iterations, graph_cleaning.max_tip_length,
                graph_cleaning.min_bubble_length_difference,
                graph_cleaning.max_bubble_length, graph_cleaning.ploidy);
    }
};

int run_assembly(const asm_config &param);
#endif
