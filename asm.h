#pragma once

#ifndef _ASM_H_
#define _ASM_H_
#include <inttypes.h>
#include <stdint.h>

#include <string>
struct syncmer_params {
    static constexpr unsigned long invalid_param_value = 127127127;
    enum class typeenum { Density, MinCodeOpen, MinCodeClosed, Invalid };
    constexpr typeenum type() const {
        if (invalid_density != invalid_param_value) {
            return typeenum::Density;
        }
        if (mincode.s != invalid_param_value) {
            if (mincode.t != invalid_param_value) {
                return typeenum::MinCodeOpen;
            } else {
                return typeenum::MinCodeClosed;
            }
        }
        return typeenum::Invalid;
    }
    unsigned k = 8;
    union {
        double density;
        unsigned long invalid_density = invalid_param_value;
    };
    struct {
        unsigned long s = invalid_param_value;
        unsigned long t = invalid_param_value;
    } mincode;
    std::string spaced_seed_mask = "11111111111111111111111111111111";
    uint64_t compute_spaced_mask() const {
        if (spaced_seed_mask.size() > 32) {
            fprintf(stderr, "Error: mask string is too long %zu > 32\n",
                    spaced_seed_mask.size());
            exit(EXIT_FAILURE);
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
                    fprintf(stderr,
                            "Error: invalid character '%c' in mask string\n",
                            c);
                    exit(EXIT_FAILURE);
            }
        }
        if ((int64_t)acc != -1) {
            fprintf(stderr,
                    "Using spaced seeds with mask %s built into 0x%" PRIX64
                    "!\n",
                    spaced_seed_mask.c_str(), acc);
        }
        return acc;
    }
};
struct asm_config {
    std::string input;
    std::string output;
    int lcp_level = 3;
    struct {
        int io = 4;
        int work = 8;
    } threads;
    size_t seed = 1453;
    syncmer_params syncmer;
    struct {
        uint64_t projected_count = 100'000'000;
        double false_positive_rate = 0.001;
    } bloom;
    struct {
        int iterations = 2;
        int max_tip_length = LCMER_L * 2;
        int min_bubble_length_difference = 2;
        double min_bubble_depth_difference_percentage = 2;
    } graph_cleaning;
};

int run_assembly(const asm_config &param);
#endif
