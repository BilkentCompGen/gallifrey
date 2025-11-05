#pragma once

#ifndef _ASM_H_
#define _ASM_H_
#include <stdint.h>

#include <string>

struct asm_parameters {
    std::string input;
    std::string output;
    int lcp_level = 3;
    int threads = 8;
    size_t random_seed = 1453;
    struct {
        int k;
        union {
            double density = 0.03;
            int s;
            int t;
        };
    } syncmer;
    struct {
        uint64_t projected_count = 100'000'000;
        double false_positive_rate = 0.001;
    } bloom;
    struct {
        int iterations = 2;
        int max_tip_length = LCMER_L * 2;
        int min_bubble_length_difference = 2;
    } graph_cleaning;
};

int run_assembly(const asm_parameters &param);
#endif
