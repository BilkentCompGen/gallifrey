#pragma once
#include <sys/_types/_u_int64_t.h>
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
        double density = 0.03;
        struct syncmer{
            int k;
            int s;
            int t;
        };
        struct {
            uint64_t projected_count = 100'000'000;
            double false_positive_rate = 0.001;
        } bloom;
};

int run_assembly(const asm_parameters &param);
#endif
