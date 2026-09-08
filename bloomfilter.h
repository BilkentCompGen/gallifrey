#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H
#include <mutex>
#define KOMIHASH_NS_CUSTOM komihash

#include <algorithm>
#include <atomic>
#include <inttypes.h>
#include <memory>
#include <vector>
#include <cmath>

// #include "ext/komihash.h"
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
using std::hardware_destructive_interference_size;
#else
// 64 bytes on x86-64 │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ __cacheline_aligned │
// ...
constexpr std::size_t hardware_constructive_interference_size = 64;
constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

struct BasicBloomFilter {
    BasicBloomFilter(int n, float eps) : n(n), epsilon(eps) {
        m = ComputeNumBits();
        k = ComputeNumHashFns();
        bv.resize((m + 7) / 8);
    }

    int ComputeNumBits() {
        return static_cast<int>(-1.44 * n * std::log2(epsilon) + 0.5);
    }

    int ComputeNumHashFns() {
        return static_cast<int>(-std::log2(epsilon) + 0.5);
    }

    void Insert(uint32_t h1, uint32_t h2) {
        for (int i = 0; i < k; i++) {
            uint32_t hash = (h1 + i * h2) % m;
            uint64_t bit_idx = hash % 8;
            uint64_t byte_idx = hash / 8;
            bv[byte_idx] |= (1 << bit_idx);
        }
    }

    bool Query(uint32_t h1, uint32_t h2) const {
        bool result = true;
        for (int i = 0; i < k; i++) {
            uint32_t hash = (h1 + i * h2) % m;
            uint64_t bit_idx = hash % 8;
            uint64_t byte_idx = hash / 8;
            result &= (bv.at(byte_idx) >> bit_idx) & 1;
        }
        return result;
    }

    void Reset() { std::fill(bv.begin(), bv.end(), 0); }

    int n;
    float epsilon;

    int m;
    int k;
    std::vector<uint8_t> bv;
};

template <int BLOCK_BYTES = hardware_constructive_interference_size>
struct BlockedBloomFilter {
    static constexpr int block_bytes = BLOCK_BYTES;
    static constexpr int block_bits = block_bytes * 8;
    // Storage word only: every access goes through the byte pointer from
    // GetBlock(), so this type just sizes and aligns the allocation.
    using block_type = unsigned long long;
    static_assert(block_bytes % sizeof(block_type) == 0,
                  "block must be a whole number of storage words");
    BlockedBloomFilter(int n, float eps) : n(n), epsilon(eps) {
        m = ComputeNumBits();
        k = ComputeNumHashFns();
        num_blocks = (m + block_bits - 1) / block_bits;
        bv.resize(num_blocks * block_bytes / sizeof(block_type));
    }

    int ComputeNumBits() {
        return static_cast<int>(-1.44 * n * std::log2(epsilon) + 0.5);
    }

    int ComputeNumHashFns() {
        return static_cast<int>(-std::log2(epsilon) + 0.5);
    }

    const uint8_t *GetBlock(uint32_t h1, uint32_t) const {
        uint32_t block_idx = h1 % num_blocks;
        uint32_t byte_idx = block_idx * block_bytes;
        return ((const uint8_t *)bv.data()) + byte_idx;
    }
    uint8_t *GetBlock(uint32_t h1, uint32_t) {
        uint32_t block_idx = h1 % num_blocks;
        uint32_t byte_idx = block_idx * block_bytes;
        return ((uint8_t *)bv.data()) + byte_idx;
    }

    void Insert(uint32_t h1, uint32_t h2) {
        uint8_t *block = GetBlock(h1, h2);
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % block_bits;
            uint64_t bit_idx = bit_pos % 8;
            uint64_t byte_idx = bit_pos / 8;
            block[byte_idx] |= (1 << bit_idx);
        }
    }

    bool Query(uint32_t h1, uint32_t h2) const {
        bool result = true;
        const uint8_t *block = GetBlock(h1, h2);
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % block_bits;
            uint64_t bit_idx = bit_pos % 8;
            uint64_t byte_idx = bit_pos / 8;
            result &= (block[byte_idx] >> bit_idx) & 1;
        }
        return result;
    }

    void Reset() { std::fill(bv.begin(), bv.end(), 0); }

    int n;
    float epsilon;

    int num_blocks;
    int m;
    int k;
    std::vector<block_type> bv;
};
template <int LOCKS_PER_2PN_BLOCKS = 0, class mutex_type = std::mutex,
          int BLOCK_BYTES = hardware_constructive_interference_size>
struct BlockedLockedBloomFilter {
    static constexpr int block_bytes = BLOCK_BYTES;
    static constexpr int block_bits = block_bytes * 8;
    static_assert(LOCKS_PER_2PN_BLOCKS < 32);
    static constexpr int locks_per_n_blocks = 1 << LOCKS_PER_2PN_BLOCKS;
    using block_type = unsigned long long;
    BlockedLockedBloomFilter(int n, float eps)
        : n(n),
          epsilon(eps),
          m{ComputeNumBits()},
          k{ComputeNumHashFns()},
          num_blocks{(m + block_bits - 1) / block_bits},
          num_locks{1 + num_blocks / locks_per_n_blocks},
          bv(num_blocks * block_bytes / sizeof(block_type)),
          locks(new mutex_type[num_locks]) {
              fprintf(stderr,
                      "BloomFilter: n=%d, eps=%.2e, m=%" PRIu64
                      ", k=%d, num_blocks=%" PRIu64 ", num_locks=%" PRIu64
                      "\n",
                      n, epsilon, m, k, num_blocks, num_locks);
          }

    uint64_t ComputeNumBits() {
        return static_cast<uint64_t>(-1.44 * n * std::log2(epsilon) + 0.5);
    }

    int ComputeNumHashFns() {
        return static_cast<int>(-std::log2(epsilon) + 0.5);
    }

    const uint8_t *GetBlock(uint32_t h1, uint32_t) const {
        uint32_t block_idx = h1 % num_blocks;
        uint32_t byte_idx = block_idx * block_bytes;
        return ((const uint8_t *)bv.data()) + byte_idx;
    }
    std::pair<uint8_t *, uint32_t> GetBlock(uint32_t h1, uint32_t) {
        uint32_t block_idx = h1 % num_blocks;
        uint32_t byte_idx = block_idx * block_bytes;
        return {(uint8_t *)bv.data() + byte_idx,
                block_idx >> LOCKS_PER_2PN_BLOCKS};
    }

    void Insert(uint32_t h1, uint32_t h2) {
        auto [block, lock_idx] = GetBlock(h1, h2);
        std::lock_guard guard(locks[lock_idx]);
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % block_bits;
            uint64_t bit_idx = bit_pos % 8;
            uint64_t byte_idx = bit_pos / 8;
            block[byte_idx] |= (1 << bit_idx);
        }
    }

    bool Query(uint32_t h1, uint32_t h2) const {
        bool result = true;
        const uint8_t *block = GetBlock(h1, h2);
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % block_bits;
            uint64_t bit_idx = bit_pos % 8;
            uint64_t byte_idx = bit_pos / 8;
            result &= (block[byte_idx] >> bit_idx) & 1;
        }
        return result;
    }

    void Reset() { std::fill(bv.begin(), bv.end(), 0); }

    int n;
    float epsilon;

    uint64_t m;
    int k;
    uint64_t num_blocks;
    uint64_t num_locks;
    std::vector<block_type> bv;
    std::unique_ptr<mutex_type[]> locks;
};
// Lock-free blocked Bloom filter. Same block geometry and sizing as
// BlockedLockedBloomFilter, but the bit array is a flat array of atomic
// words: Insert uses fetch_or, Query uses load (both relaxed; nothing is
// published through the filter). No per-block mutex array is needed;
// atomic<block_type> is the same 8 bytes as the plain word.
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

    int64_t n;
    float epsilon;

    uint64_t m;
    int k;
    uint64_t num_blocks;
    uint64_t num_words;
    std::unique_ptr<std::atomic<block_type>[]> bv;
};

template <int Compensation>
struct RegisterBlockedBloomFilter {
    RegisterBlockedBloomFilter(int n, float eps) : n(n), epsilon(eps) {
        m = ComputeNumBits();
        k = ComputeNumHashFns();
        num_blocks = (m + 64 - 1) / 64;
        bv.resize(num_blocks);
    }

    int ComputeNumBits() {
        auto bits_per_val = -1.44 * std::log2(epsilon) + Compensation;
        return static_cast<int>(bits_per_val * n + 0.5);
    }

    int ComputeNumHashFns() {
        return static_cast<int>(-std::log2(epsilon) + 0.5);
    }

    uint64_t *GetBlock(uint32_t h1, uint32_t) {
        uint32_t block_idx = h1 % num_blocks;
        return &bv[block_idx];
    }
    const uint64_t *GetBlock(uint32_t h1, uint32_t) const {
        uint32_t block_idx = h1 % num_blocks;
        return &bv.at(block_idx);
    }

    uint64_t ConstructMask(uint32_t h1, uint32_t h2) const {
        uint64_t mask = 0;
        for (int i = 1; i < k; i++) {
            uint32_t bit_pos = (h1 + i * h2) % 64;
            mask |= (1ull << bit_pos);
        }
        return mask;
    }

    void Insert(uint32_t h1, uint32_t h2) {
        uint64_t *block = GetBlock(h1, h2);
        *block |= ConstructMask(h1, h2);
    }

    bool Query(uint32_t h1, uint32_t h2) const {
        const uint64_t *block = GetBlock(h1, h2);
        uint64_t mask = ConstructMask(h1, h2);
        return (*block & mask) == mask;
    }

    void Reset() { std::fill(bv.begin(), bv.end(), 0); }

    int n;
    float epsilon;

    int num_blocks;
    int m;
    int k;
    std::vector<uint64_t> bv;
};
}  // namespace bloom_ext

template <class BF>
concept BF_Insert_Query_Concept = requires(BF bf) {
    { bf.Query(uint32_t{}, uint32_t{}) } -> std::same_as<bool>;
    { bf.Insert(uint32_t{}, uint32_t{}) } -> std::same_as<void>;
    // { BF(double{}, uint64_t{}) } -> BF;
};
struct null_mutex {
    void lock() {}
    void try_lock() {}
    void unlock() {}

    nullptr_t native_handle() { return nullptr; }
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

template <BF_Insert_Query_Concept BF>
struct bloom_interface<BF> {
    BF implementation;
    uint64_t random_seed;

    bloom_interface(double false_positive_rate, size_t projected_element_count,
                    uint64_t random_seed)
        : implementation(projected_element_count, false_positive_rate),
          random_seed(random_seed) {}
    bool contains(const uint64_t &h) const {
        return implementation.Query(h >> 32, h & 0xFFFFFFFF);
    }
    void insert(const uint64_t &h) {
        implementation.Insert(h >> 32, h & 0xFFFFFFFF);
    }
};

#endif
