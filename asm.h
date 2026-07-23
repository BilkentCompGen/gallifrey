#pragma once

#ifndef _ASM_H_
#define _ASM_H_

#include <inttypes.h>
#include <stdint.h>

#include <cassert>
#include <expected>
#include <optional>
#include <string>
#include <vector>
struct syncmer_params {
    enum class typeenum { Density, MinCodeOpen, MinCodeClosed, Invalid };
    constexpr typeenum type() const {
        if (density.has_value()) {
            return typeenum::Density;
        }
        if (mincode.s.has_value()) {
            if (mincode.t.has_value()) {
                return typeenum::MinCodeOpen;
            } else {
                return typeenum::MinCodeClosed;
            }
        }
        return typeenum::Invalid;
    }
    unsigned k = 8;
    std::optional<double> density;
    struct {
        std::optional<unsigned long> s;
        std::optional<unsigned long> t;
    } mincode;
    std::string spaced_seed_mask = "11111111111111111111111111111111";

    // Validating mask builder: returns the packed mask or an error message.
    // This is the single source of truth for mask parsing; the total
    // compute_spaced_mask() below delegates here.
    std::expected<uint64_t, std::string> compute_spaced_mask_checked() const {
        if (spaced_seed_mask.size() > 32) {
            return std::unexpected(
                "mask string is too long " +
                std::to_string(spaced_seed_mask.size()) + " > 32");
        }
        uint64_t acc = -1;
        for (char c : spaced_seed_mask) {
            acc <<= 2;
            switch (c) {
                case '0':
                    break;
                case '1':
                    acc |= 3;
                    break;
                case '2':
                    acc |= 2;
                    break;
                case '3':
                    acc |= 1;
                    break;
                default:
                    return std::unexpected(
                        std::string("invalid character '") + c +
                        "' in mask string");
            }
        }
        return acc;
    }

    // Total accessor for the hot construction paths (member-init lists). The
    // config is validated once up front via validate(); reaching the error
    // branch here means that contract was broken.
    uint64_t compute_spaced_mask() const {
        auto r = compute_spaced_mask_checked();
        if (!r) {
            assert(false && "mask must be validated before use");
            __builtin_unreachable();
        }
        if ((int64_t)*r != -1) {
            fprintf(stderr,
                    "Using spaced seeds with mask %s built into 0x%" PRIX64
                    "!\n",
                    spaced_seed_mask.c_str(), *r);
        }
        return *r;
    }

    // Validate the whole syncmer configuration once, before any syncmer
    // functor is constructed. Lets callers surface a clean error instead of
    // exit()-ing from deep inside a constructor.
    std::expected<void, std::string> validate() const {
        if (type() == typeenum::Invalid) {
            return std::unexpected(
                std::string("Provide a k-mer selection: either --density or "
                            "--submer-s (optionally with --submer-t)"));
        }
        if (auto r = compute_spaced_mask_checked(); !r) {
            return std::unexpected(r.error());
        }
        // Cross-field numeric checks: reject configurations that select no
        // k-mers or index an s-window out of the k-mer before any functor is
        // constructed (a per-option CLI::Range can't express these relations).
        if (density.has_value() && *density == 0.0) {
            return std::unexpected(std::string(
                "--density must be > 0 (0 selects no k-mers -> empty "
                "assembly)"));
        }
        if (mincode.s.has_value()) {
            if (*mincode.s >= k) {
                return std::unexpected(
                    "--submer-s (" + std::to_string(*mincode.s) +
                    ") must be < --syncmer-k (" + std::to_string(k) + ")");
            }
            if (mincode.t.has_value() && *mincode.t >= k - *mincode.s + 1) {
                return std::unexpected(
                    "--submer-t (" + std::to_string(*mincode.t) +
                    ") must be < k - s + 1 (" +
                    std::to_string(k - *mincode.s + 1) + ")");
            }
        }
        return {};
    }
    void print_cfg(FILE *f) const {
        fprintf(f, "spaced-seed-mask: \"%s\"\n", spaced_seed_mask.c_str());
        fprintf(f, "k: %d\n", k);
        if (density.has_value()) {
            fprintf(f, "kmer-selection: \"density\"\n");
            fprintf(f, "Density: %f\n", *density);
        } else if (mincode.s.has_value()) {
            if (mincode.t.has_value()) {
                fprintf(f, "kmer-selection: \"open-syncmers\"\n");
                fprintf(f, "syncmer-s:%lu\nsyncmer-t:%lu\n", *mincode.s,
                        *mincode.t);
            } else {
                fprintf(f, "kmer-selection: \"closed-syncmers\"\n");
                fprintf(f, "syncmer-s: %lu\n", *mincode.s);
            }
        } else {
            fprintf(f, "kmer-selection: \"none\"\n");
        }
    }
};
struct asm_config {
    std::vector<std::string> inputs;
    std::string output;
    int lcp_level = 3;
    int dbg_lcmer_l = 8;
    int olc_lcmer_l = 8;
    struct {
        int io = 4;
        int work = 8;
    } threads;
    size_t seed = 1453;
    bool rle = true;
    bool parallel_file_reading = true;
    syncmer_params syncmer;
    bool print_precleanup_graph = false;
    bool print_prelayout_graph = false;
    bool silent = false;
    struct {
        uint64_t projected_count = 100'000'000;
        double false_positive_rate = 0.001;
    } bloom;
    struct {
        int iterations = 2;
        int max_tip_length = 8 * 2;
        int min_bubble_length_difference = 2;
        double min_bubble_depth_difference_percentage = 2;
        bool enable_loop_opening = false;
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
        fprintf(f, "DBG LCMER Length: %d\n", dbg_lcmer_l);
        fprintf(f, "OLC LCMER Length: %d\n", olc_lcmer_l);
        fprintf(f, "Threads:\n\tIO: %d\n\tWork: %d\n", threads.io, threads.work);
        fprintf(f, "Seed: %zu\n", seed);
        fprintf(f, "RLE: %s\n", rle ? "true" : "false");
        fprintf(f, "File Reading Mode: %s\n",
                parallel_file_reading ? "parallel" : "sequential");
        syncmer.print_cfg(f);
        fprintf(f,
                "Bloom Filter Parameters:\n\tProjected Count: %" PRIu64
                "\n\tFalse Positive "
                "Rate: %.3f\n",
                bloom.projected_count, bloom.false_positive_rate);
        fprintf(f,
                "Graph Cleaning Parameters: \n\tIterations: %d\n\tMax Tip "
                "Length: %d\n\tMin Bubble Length Difference: %d\n\tMin Bubble "
                "Depth Difference Rate: %.2f\n\tOpen Simple Loops: %s\n",
                graph_cleaning.iterations, graph_cleaning.max_tip_length,
                graph_cleaning.min_bubble_length_difference,
                graph_cleaning.min_bubble_depth_difference_percentage,
                graph_cleaning.enable_loop_opening ? "true" : "false");
    }
};


int run_assembly(const asm_config &param);
#endif
