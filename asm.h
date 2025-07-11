#pragma once
#ifndef _ASM_H_
#define _ASM_H_
#include "ext/bloom_filter.hpp"
int run_assembly(const std::string &input, const std::string &output, int threads, int lcp_level, bloom_parameters bparam);
#endif