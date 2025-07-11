#pragma once
#ifndef _ASM_H_
#define _ASM_H_
#include <stdint.h>
#include <string>
struct bloom_param{
    uint64_t projected_count;
    double false_positive_rate;
    uint64_t seed;
}; 
int run_assembly(const std::string &input, const std::string &output, int threads, int lcp_level, bloom_param bparam);
#endif