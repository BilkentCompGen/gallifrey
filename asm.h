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
// Submer selection parameters. A submer is a k-mer selected by a rule:
// a mincode submer (select if the hashed code < H/c), or a syncmer
// (select by the position of the smallest s-mer: closed = first or last
// s-mer, open = the s-mer at offset t).
struct syncmer_params {
    enum class typeenum { MincodeSubmer, OpenSyncmer, ClosedSyncmer, Invalid };
    constexpr typeenum type() const {
        if (density.has_value()) {
            return typeenum::MincodeSubmer;
        }
        if (smer.s.has_value()) {
            if (smer.t.has_value()) {
                return typeenum::OpenSyncmer;
            } else {
                return typeenum::ClosedSyncmer;
            }
        }
        return typeenum::Invalid;
    }
    unsigned k = 8;
    // Mincode-submer parameter: target density d (compression factor c = 1/d).
    std::optional<double> density;
    // Syncmer parameters: s-mer length s and (open syncmer) offset t.
    struct {
        std::optional<unsigned long> s;
        std::optional<unsigned long> t;
    } smer;
    // Apply the rule to the strand-canonical k-mer (min of the k-mer and its
    // reverse complement) so a k-mer and its reverse complement are co-selected.
    bool canonical = false;
    std::string spaced_seed_mask = "11111111111111111111111111111111";

    // Validating mask builder: returns the packed mask or an error message.
    // Single source of truth for mask parsing; compute_spaced_mask() delegates
    // here.
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

    // Total accessor for member-init lists. The config is validated up front
    // via validate(); the error branch here means that contract was broken.
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

    // Effective k after spaced-seed compaction: each '0' mask char removes
    // its position, so selection runs on the dense (k - #zeros)-mer. Only
    // positions within the k-mer window count (chars beyond k are shifted
    // above the window and never apply). '2'/'3' keep their position.
    std::expected<unsigned, std::string> effective_k_checked() const {
        auto r = compute_spaced_mask_checked();
        if (!r) {
            return std::unexpected(r.error());
        }
        unsigned kept = 0;
        for (unsigned i = 0; i < k; ++i) {
            kept += ((*r >> (2 * i)) & 3) != 0;
        }
        return kept;
    }

    // Validate the whole syncmer configuration once, before any syncmer
    // functor is constructed, so callers can surface a clean error.
    std::expected<void, std::string> validate() const {
        if (type() == typeenum::Invalid) {
            return std::unexpected(
                std::string("Provide a k-mer selection: either --density or "
                            "--submer-s (optionally with --submer-t)"));
        }
        auto ek = effective_k_checked();
        if (!ek) {
            return std::unexpected(ek.error());
        }
        const unsigned k_eff = *ek;
        if (k_eff == 0) {
            return std::unexpected(
                "--spaced-seed-mask removes every k-mer position (effective "
                "k = 0)");
        }
        // Cross-field numeric checks: reject configs that select no k-mers or
        // index an s-window out of the compacted k-mer. s/t bounds are against
        // the effective k (masked positions do not exist in the selection).
        if (density.has_value() && *density == 0.0) {
            return std::unexpected(std::string(
                "--density must be > 0 (0 selects no k-mers -> empty "
                "assembly)"));
        }
        if (smer.s.has_value()) {
            if (*smer.s >= k_eff) {
                return std::unexpected(
                    "--submer-s (" + std::to_string(*smer.s) +
                    ") must be < effective k (" + std::to_string(k_eff) +
                    " = --syncmer-k " + std::to_string(k) + " minus " +
                    std::to_string(k - k_eff) + " masked positions)");
            }
            if (smer.t.has_value() && *smer.t >= k_eff - *smer.s + 1) {
                return std::unexpected(
                    "--submer-t (" + std::to_string(*smer.t) +
                    ") must be < effective-k - s + 1 (" +
                    std::to_string(k_eff - *smer.s + 1) + ")");
            }
        }
        return {};
    }
    void print_cfg(FILE *f) const {
        fprintf(f, "spaced-seed-mask: \"%s\"\n", spaced_seed_mask.c_str());
        fprintf(f, "k: %d\n", k);
        if (auto ek = effective_k_checked(); ek && *ek != k) {
            fprintf(f, "effective-k (after spaced-seed compaction): %u\n",
                    *ek);
        }
        if (density.has_value()) {
            fprintf(f, "submer-rule: \"mincode-submer\"\n");
            fprintf(f, "density: %f\n", *density);
        } else if (smer.s.has_value()) {
            if (smer.t.has_value()) {
                fprintf(f, "submer-rule: \"open-syncmer\"\n");
                fprintf(f, "syncmer-s: %lu\nsyncmer-t: %lu\n", *smer.s,
                        *smer.t);
            } else {
                fprintf(f, "submer-rule: \"closed-syncmer\"\n");
                fprintf(f, "syncmer-s: %lu\n", *smer.s);
            }
        } else {
            fprintf(f, "submer-rule: \"none\"\n");
        }
        fprintf(f, "canonical: %s\n", canonical ? "true" : "false");
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
        int iterations = 6;
        int max_tip_length = 8 * 2;
        double min_bubble_length_difference = 2;
        size_t max_bubble_arm_length = 128;
        double min_bubble_depth_difference_percentage = 2;
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
                "Length: %d\n\tMin Bubble Length Difference: %.2f\n\tMax "
                "Bubble Arm Length: %zu\n\tMin Bubble Depth Difference "
                "Rate: %.2f\n",
                graph_cleaning.iterations, graph_cleaning.max_tip_length,
                graph_cleaning.min_bubble_length_difference,
                graph_cleaning.max_bubble_arm_length,
                graph_cleaning.min_bubble_depth_difference_percentage);
    }
};


int run_assembly(const asm_config &param);
#endif
