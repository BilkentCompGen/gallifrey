#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H
#define KOMIHASH_NS_CUSTOM komihash
#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <random>
#include <immintrin.h>
#include "ext/komihash.h"
namespace bloom_ext{
#include "ext/bloom_filter.hpp"
#include "ext/savebuffer_bh.h"
}

template<class BF>
concept BF_Insert_Query_Concept = 
    requires  (BF bf) {
        { bf.Query(uint32_t{}, uint32_t{}) } -> std::same_as<bool>;
        { bf.Insert(uint32_t{}, uint32_t{}) } -> std::same_as<void>;
        // { BF(double{}, uint64_t{}) } -> BF;
    };


template<class BloomImplementation>
struct bloom_interface{
    double false_positive_rate;
    double projected_element_count;
    unsigned int random_seed;
    BloomImplementation implementation;
    bloom_interface(double false_positive_rate, double projected_element_count, unsigned int random_seed) :
     false_positive_rate(false_positive_rate), projected_element_count(projected_element_count), random_seed(random_seed), implementation(false_positive_rate, projected_element_count, random_seed) {}
    bool contains(const auto &element) const{
        return implementation.contains(element);
    }
    bool insert(const auto &element){
        return implementation.insert(element);
    }
};



template<>
struct bloom_interface<bloom_ext::bloom_filter> {
    bloom_ext::bloom_parameters parameters;
    bloom_ext::bloom_filter implementation;

    bloom_interface(double false_positive_rate, size_t projected_element_count, uint64_t random_seed) 
    {
        parameters.false_positive_probability = false_positive_rate;
        parameters.projected_element_count = projected_element_count;
        parameters.random_seed = random_seed;
        parameters.compute_optimal_parameters();
        implementation = bloom_ext::bloom_filter(parameters);
    }
    bool contains(const auto &element) const{
        return implementation.contains(element);
    }
    void insert(const auto &element){
        implementation.insert(element);
    }
};

template<BF_Insert_Query_Concept BF>
struct bloom_interface<BF>{
    BF implementation;
    uint64_t random_seed;
    bloom_interface(double false_positive_rate, size_t projected_element_count, uint64_t random_seed) : 
        implementation(projected_element_count, false_positive_rate) {}
    bool contains(const auto &element) const{
        uint64_t h = komihash::komihash(&element, sizeof(element), random_seed);

        return implementation.Query(h >> 32, h & 0xFFFFFFFF);
    }
    void insert(const auto &element){
        uint64_t h = komihash::komihash(&element, sizeof(element), random_seed);
        implementation.Insert(h >> 32, h & 0xFFFFFFFF);
    }
};

#endif