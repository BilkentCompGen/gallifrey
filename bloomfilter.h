#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H

#include <inttypes.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <concepts>
#include <memory>
#include <numbers>

// MurmurHash3 64-bit finalizer.
constexpr uint64_t fmix64(uint64_t k) {
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33;
    return k;
}

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

    // Bits for n keys at false-positive rate epsilon: the smallest block
    // count whose modelled rate (ModelFpr) at n / num_blocks keys per block
    // is at most epsilon. Keep the math in double, then clamp: a very large n
    // must not overflow the double->uint64 cast (UB) nor request an absurd
    // allocation. Without probes (epsilon near 1, the filter is off) the
    // unblocked optimum n log2(1/epsilon) / ln 2 is kept.
    uint64_t ComputeNumBits() {
        // Ceiling of about 1.2 TB of filter.
        constexpr double kMaxBits = 1e13;
        const int probes = ComputeNumHashFns() - 1;
        double bits =
            -static_cast<double>(n) * std::log2(epsilon) / std::numbers::ln2 +
            0.5;
        if (probes >= 1) {
            const double blocks = std::ceil(static_cast<double>(n) /
                                            MaxLoadForFpr(probes, epsilon));
            bits = blocks * block_bits;
        }
        bits = std::clamp(bits, 1.0, kMaxBits);
        return static_cast<uint64_t>(bits);
    }

    int ComputeNumHashFns() {
        return static_cast<int>(-std::log2(epsilon) + 0.5);
    }

    // Expected false-positive rate of a filter holding lambda keys per block
    // on average, each setting `probes` independent uniform positions of its
    // block. Block loads are Poisson(lambda); a block holding j keys has s
    // bits set with the distribution of j * probes uniform throws into
    // block_bits cells, and a fresh key is a false positive with probability
    // (s / block_bits)^probes. Exact for this design, so it also prices the
    // block-load variance that the unblocked fill^probes ignores.
    static double ModelFpr(int probes, double lambda) {
        if (probes < 1) return 1.0;
        if (lambda <= 0) return 0.0;
        constexpr int B = block_bits;
        std::array<double, B + 1> dist{};
        std::array<double, B + 1> hit{};
        for (int s = 0; s <= B; ++s)
            hit[s] = std::pow(static_cast<double>(s) / B, probes);
        dist[0] = 1.0;
        // Poisson mass beyond lambda + 12 sd + 30 is far below double
        // precision for every lambda.
        constexpr double kTailSds = 12.0;
        constexpr double kTailPad = 30.0;
        const int jmax =
            static_cast<int>(lambda + kTailSds * std::sqrt(lambda) + kTailPad);
        const double log_lambda = std::log(lambda);
        double fpr = 0.0;
        for (int j = 0; j <= jmax; ++j) {
            const double w =
                std::exp(-lambda + j * log_lambda - std::lgamma(j + 1.0));
            double h = 0.0;
            for (int s = 0; s <= B; ++s) h += dist[s] * hit[s];
            fpr += w * h;
            // One more key: `probes` throws, each lands on an unset bit with
            // probability (B - s) / B. Descending s keeps the update in place.
            for (int t = 0; t < probes; ++t) {
                for (int s = B - 1; s >= 0; --s) {
                    dist[s + 1] += dist[s] * (B - s) / B;
                    dist[s] *= static_cast<double>(s) / B;
                }
            }
        }
        return fpr;
    }

    // Largest mean keys per block whose ModelFpr is at most eps, by bisection
    // on a log scale.
    static double MaxLoadForFpr(int probes, double eps) {
        // Below one key per million blocks every rate of interest is met.
        constexpr double kMinLoad = 1e-6;
        // At 4 block_bits keys per block practically every bit is set.
        constexpr double kMaxLoad = 4.0 * block_bits;
        // 60 halvings of the log range reach double precision.
        constexpr int kSteps = 60;
        double lo = kMinLoad, hi = kMaxLoad;
        if (ModelFpr(probes, hi) <= eps) return hi;
        for (int i = 0; i < kSteps; ++i) {
            const double mid = std::sqrt(lo * hi);
            if (ModelFpr(probes, mid) <= eps)
                lo = mid;
            else
                hi = mid;
        }
        return lo;
    }

    // Probe positions of a key: its block's first word and the block's bits
    // to set or test as one mask per word. The block comes from the top bits
    // of fmix64(key) by multiply-shift; the k - 1 positions are independent
    // uniform draws of log2(block_bits) bits each, the first two from the
    // mixed key's low bits (the block barely depends on them), the rest from
    // fmix64 of the mixed key plus multiples of a Weyl constant. Two
    // positions may coincide. Positions drawn from an arithmetic progression
    // instead allow only block_bits * block_bits / 2 patterns per block, and
    // the keys of a block that share a pattern are certain false positives
    // of each other, which floors the false-positive rate however large the
    // filter.
    static constexpr int pos_bits = std::countr_zero(unsigned{block_bits});
    static_assert(std::has_single_bit(unsigned{block_bits}));
    // The block comes from the mixed key's top bits, so its low bits are free
    // for two positions (18 of 64 bits; the block depends on them only
    // through a carry).
    static constexpr int kPositionsFromKey = 2;
    // Positions one further fmix64 provides.
    static constexpr int kPositionsPerMix = 64 / pos_bits;
    // 2^64 / golden ratio: successive salts are far apart.
    static constexpr uint64_t kWeyl = 0x9e3779b97f4a7c15ULL;
    // A GCC/Clang extension; __extension__ keeps -Wpedantic quiet about it.
    __extension__ typedef unsigned __int128 u128;
    struct Probes {
        uint64_t base;
        std::array<block_type, words_per_block> mask;
    };
    Probes ProbesOf(uint64_t key) const {
        const uint64_t x = fmix64(key);
        const auto block =
            static_cast<uint64_t>((static_cast<u128>(x) * num_blocks) >> 64);
        Probes p{block * words_per_block, {}};
        uint64_t bits = x;
        uint64_t salt = x;
        int avail = kPositionsFromKey;
        for (int i = 1; i < k; i++) {
            if (avail == 0) {
                salt += kWeyl;
                bits = fmix64(salt);
                avail = kPositionsPerMix;
            }
            const auto pos = static_cast<uint32_t>(bits & (block_bits - 1));
            bits >>= pos_bits;
            --avail;
            p.mask[pos / word_bits] |= block_type(1) << (pos % word_bits);
        }
        return p;
    }

    void Insert(uint64_t key) { InsertProbes(ProbesOf(key)); }

    bool Query(uint64_t key) const { return QueryProbes(ProbesOf(key)); }

    // Inserts key; returns whether it was absent before. Mixes once, reads
    // the whole block without branching (a present key, the common case,
    // stops there) and fetch_ors only the words that miss a bit.
    bool InsertIfAbsent(uint64_t key) {
        const Probes p = ProbesOf(key);
        std::array<block_type, words_per_block> missing;
        block_type any = 0;
        for (int w = 0; w < words_per_block; w++) {
            missing[w] =
                ~bv[p.base + w].load(std::memory_order_relaxed) & p.mask[w];
            any |= missing[w];
        }
        if (!any) return false;
        for (int w = 0; w < words_per_block; w++) {
            if (missing[w])
                bv[p.base + w].fetch_or(missing[w], std::memory_order_relaxed);
        }
        return true;
    }

    void InsertProbes(const Probes &p) {
        for (int w = 0; w < words_per_block; w++) {
            if (p.mask[w])
                bv[p.base + w].fetch_or(p.mask[w], std::memory_order_relaxed);
        }
    }

    bool QueryProbes(const Probes &p) const {
        bool result = true;
        for (int w = 0; w < words_per_block; w++) {
            result &= (bv[p.base + w].load(std::memory_order_relaxed) &
                       p.mask[w]) == p.mask[w];
        }
        return result;
    }

    // Fraction of filter bits set.
    double fill_fraction() const {
        uint64_t set = 0;
        for (uint64_t i = 0; i < num_words; i++)
            set += std::popcount(bv[i].load(std::memory_order_relaxed));
        return static_cast<double>(set) /
               (static_cast<double>(num_words) * word_bits);
    }

    // Mean keys per block for a fill: a key leaves a given bit unset with
    // probability (1 - 1 / block_bits)^(k - 1) and block loads are Poisson,
    // so fill = 1 - exp(-lambda * (1 - (1 - 1 / block_bits)^(k - 1))).
    double load_for_fill(double fill) const {
        const double per_key = 1.0 - std::pow(1.0 - 1.0 / block_bits, k - 1);
        return -std::log1p(-fill) / per_key;
    }

    // Distinct keys inserted, from the fill (exact in expectation).
    double distinct_keys_for_fill(double fill) const {
        return load_for_fill(fill) * static_cast<double>(num_blocks);
    }

    // Modelled false-positive rate of a query at this fill.
    double fpr_for_fill(double fill) const {
        return ModelFpr(k - 1, load_for_fill(fill));
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
    { bf.Query(uint64_t{}) } -> std::same_as<bool>;
    { bf.Insert(uint64_t{}) } -> std::same_as<void>;
    { bf.InsertIfAbsent(uint64_t{}) } -> std::same_as<bool>;
};

// The filter is seed-independent: the probe positions come from the 64-bit
// lcmer key alone (via fmix64), so for a given key set
// --random-seed does not change which keys collide in the filter. The seed
// argument is accepted and unused; threading it into the probes would change
// the assembly.
template <BF_Insert_Query_Concept BF>
struct bloom_interface {
    BF implementation;

    bloom_interface(double false_positive_rate, size_t projected_element_count,
                    uint64_t /*random_seed*/)
        : implementation(projected_element_count, false_positive_rate) {}
    bool contains(const uint64_t &h) const { return implementation.Query(h); }
    void insert(const uint64_t &h) { implementation.Insert(h); }
    // Inserts h; returns whether it was absent before.
    bool insert_if_absent(const uint64_t &h) {
        return implementation.InsertIfAbsent(h);
    }
};

#endif
