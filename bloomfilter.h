#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H

#include <inttypes.h>

#include <atomic>
#include <bit>
#include <cmath>
#include <concepts>
#include <memory>

namespace bloom_ext {
// Copyright (c) 2023 Sasha Krassovsky

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// https://save-buffer.github.io/bloom_filter.html

#ifdef __cpp_lib_hardware_interference_size
using std::hardware_constructive_interference_size;
#else
// 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │
// ...
constexpr std::size_t hardware_constructive_interference_size = 64;
#endif

// Lock-free blocked Bloom filter (the blocked layout of the post above): the
// bit array is a flat array of atomic words, Insert uses fetch_or and Query
// uses load (both relaxed; nothing is published through the filter).
template <int BLOCK_BYTES = hardware_constructive_interference_size>
struct BlockedAtomicBloomFilter {
    static constexpr int block_bytes = BLOCK_BYTES;
    static constexpr int block_bits = block_bytes * 8;
    using block_type = unsigned long long;
    static constexpr int word_bits = sizeof(block_type) * 8;
    static constexpr int words_per_block = block_bytes / sizeof(block_type);
    static_assert(std::atomic<block_type>::is_always_lock_free,
                  "atomic block word must be lock-free for this filter to be "
                  "useful");

    BlockedAtomicBloomFilter(int64_t n, float eps)
        : n(n),
          epsilon(eps),
          m{ComputeNumBits()},
          k{ComputeNumHashFns()},
          num_blocks{(m + block_bits - 1) / block_bits},
          num_words{num_blocks * words_per_block},
          bv{std::make_unique<std::atomic<block_type>[]>(num_words)} {
        fprintf(stderr,
                "BloomFilter(atomic): n=%" PRId64 ", eps=%.2e, m=%" PRIu64
                ", k=%d, num_blocks=%" PRIu64 ", num_words=%" PRIu64 "\n",
                n, epsilon, m, k, num_blocks, num_words);
    }

    uint64_t ComputeNumBits() {
        // Keep the log math in double, then clamp: a very large n must not
        // overflow the double->uint64 cast (UB) nor request an absurd
        // allocation. Ceiling ~1e13 bits.
        double bits = -1.44 * static_cast<double>(n) * std::log2(epsilon) + 0.5;
        if (bits < 1.0) bits = 1.0;
        constexpr double kMaxBits = 1e13;
        if (bits > kMaxBits) bits = kMaxBits;
        return static_cast<uint64_t>(bits);
    }

    int ComputeNumHashFns() {
        return static_cast<int>(-std::log2(epsilon) + 0.5);
    }

    uint64_t BlockWordBase(uint32_t h1) const {
        return static_cast<uint64_t>(h1 % num_blocks) * words_per_block;
    }

    void Insert(uint32_t h1, uint32_t h2) {
        const uint64_t base = BlockWordBase(h1);
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % block_bits;
            bv[base + bit_pos / word_bits].fetch_or(
                block_type(1) << (bit_pos % word_bits),
                std::memory_order_relaxed);
        }
    }

    bool Query(uint32_t h1, uint32_t h2) const {
        bool result = true;
        const uint64_t base = BlockWordBase(h1);
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % block_bits;
            result &= (bv[base + bit_pos / word_bits].load(
                           std::memory_order_relaxed) >>
                       (bit_pos % word_bits)) &
                      1;
        }
        return result;
    }

    void Reset() {
        for (uint64_t i = 0; i < num_words; i++)
            bv[i].store(0, std::memory_order_relaxed);
    }

    // Fraction of filter bits set; Query false-positive rate is about
    // fill^(k - 1) (Insert and Query probe k - 1 bits).
    double fill_fraction() const {
        uint64_t set = 0;
        for (uint64_t i = 0; i < num_words; i++)
            set += std::popcount(bv[i].load(std::memory_order_relaxed));
        return static_cast<double>(set) /
               (static_cast<double>(num_words) * word_bits);
    }

    // Distinct keys inserted, estimated from fill_fraction(): each insert
    // sets k - 1 uniform bits, so fill = 1 - exp(-(k - 1) * keys / bits).
    // The block index and the probe offsets both derive from h1, so with an
    // even num_blocks a block's keys share an h1 parity and reach fewer of
    // its bits: the fill grows slower and this reads several percent low.
    double distinct_keys_for_fill(double fill) const {
        const double bits = static_cast<double>(num_words) * word_bits;
        return -bits * std::log1p(-fill) / (k - 1);
    }

    int64_t n;
    float epsilon;

    uint64_t m;
    int k;
    uint64_t num_blocks;
    uint64_t num_words;
    std::unique_ptr<std::atomic<block_type>[]> bv;
};
}  // namespace bloom_ext

template <class BF>
concept BF_Insert_Query_Concept = requires(BF bf) {
    { bf.Query(uint32_t{}, uint32_t{}) } -> std::same_as<bool>;
    { bf.Insert(uint32_t{}, uint32_t{}) } -> std::same_as<void>;
    // { BF(double{}, uint64_t{}) } -> BF;
};

struct MockBF {
    double false_positive_rate;
    double projected_element_count;
    unsigned int random_seed;

    MockBF(double false_positive_rate,
           unsigned long long projected_element_count,
           unsigned long long random_seed)
        : false_positive_rate(false_positive_rate),
          projected_element_count(projected_element_count),
          random_seed(random_seed) {}
    bool contains(const auto &) const { return true; }
    bool insert(const auto &) { return false; }
};

template <class BloomImplementation>
struct bloom_interface {
    double false_positive_rate;
    double projected_element_count;
    unsigned int random_seed;

    BloomImplementation implementation;
    bloom_interface(double false_positive_rate, double projected_element_count,
                    unsigned int random_seed)
        : false_positive_rate(false_positive_rate),
          projected_element_count(projected_element_count),
          random_seed(random_seed),
          implementation(false_positive_rate, projected_element_count,
                         random_seed) {}
    bool contains(const auto &element) const {
        return implementation.contains(element);
    }
    bool insert(const auto &element) { return implementation.insert(element); }
};

// The filter is seed-independent: the probe positions come from the 64-bit
// lcmer key alone (h1 = high word, h2 = low word), so for a given key set
// --random-seed does not change which keys collide in the filter. The seed
// argument only keeps the constructor uniform with bloom_interface<MockBF>;
// threading it into the probes would change the assembly.
template <BF_Insert_Query_Concept BF>
struct bloom_interface<BF> {
    BF implementation;

    bloom_interface(double false_positive_rate, size_t projected_element_count,
                    uint64_t /*random_seed*/)
        : implementation(projected_element_count, false_positive_rate) {}
    bool contains(const uint64_t &h) const {
        return implementation.Query(h >> 32, h & 0xFFFFFFFF);
    }
    void insert(const uint64_t &h) {
        implementation.Insert(h >> 32, h & 0xFFFFFFFF);
    }
};

#endif
