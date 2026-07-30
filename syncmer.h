#ifndef _SYNCMER_H_
#define _SYNCMER_H_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <limits>

#include "asm.h"
#include "gtl/vector.hpp"
#include "kmer_bits.h"
#define KOMIHASH_NS_CUSTOM komihash
#include <komihash.h>
namespace syncmer {
template <class Rule>
concept SyncmerRule = requires(Rule k) {
    { k(unsigned{}) } -> std::same_as<bool>;
};
template <class HF, class INT_TYPE>
concept HashFunction = requires(HF h) {
    { HF{INT_TYPE{}} } -> std::same_as<HF>;
    { h(INT_TYPE{}) } -> std::same_as<INT_TYPE>;
};

struct Komi {
    uint64_t seed;
    Komi(uint64_t seed) : seed{seed} {}
    uint64_t operator()(uint64_t val) const {
        return komihash::komihash(&val, sizeof(uint64_t), seed);
    }
};
// uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed);
struct MurMur {
    static uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed) {
        const uint8_t *data = (const uint8_t *)key;
        const int nblocks = len / 4;

        uint32_t h1 = seed;

        const uint32_t c1 = 0xcc9e2d51;
        const uint32_t c2 = 0x1b873593;

        // Body: Process blocks of 4 bytes at a time
        const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

        for (int i = -nblocks; i; i++) {
            uint32_t k1 = blocks[i];

            k1 *= c1;
            k1 = (k1 << 15) | (k1 >> (32 - 15));
            k1 *= c2;

            h1 ^= k1;
            h1 = (h1 << 15) | (h1 >> (32 - 15));
            h1 = h1 * 5 + 0xe6546b64;
        }

        // Tail: Process remaining bytes
        const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

        uint32_t k1 = 0;

        switch (len & 3) {
            case 3:
                k1 ^= tail[2] << 16;
                [[fallthrough]];
            case 2:
                k1 ^= tail[1] << 8;
                [[fallthrough]];
            case 1:
                k1 ^= tail[0];
                k1 *= c1;
                k1 = (k1 << 15) | (k1 >> (32 - 15));
                k1 *= c2;
                h1 ^= k1;
        }

        // Finalization: Mix the hash to ensure the last few bits are fully
        // mixed
        h1 ^= len;

        /* fmix32 */
        h1 ^= h1 >> 16;
        h1 *= 0x85ebca6b;
        h1 ^= h1 >> 13;
        h1 *= 0xc2b2ae35;
        h1 ^= h1 >> 16;
        return h1;
    }
    unsigned seed;
    MurMur(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const {
        return MurmurHash3_32((const char *)&value, sizeof(value), seed);
    }
};
struct Identity {
    unsigned seed;
    Identity(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const { return value; }
};
// Builds a LuT for the given submer rule; works fine for K <= 11, 12.
// Dispatches all three submer rules (Edgar 2021 terminology): a mincode submer
// (code < H/c) and closed/open syncmers (selected by the position of the
// smallest s-mer µ(κ,s)).
template <class INT_TYPE, HashFunction<INT_TYPE> HASH_FUNCTION>
struct GenericSubmerRule {
    static constexpr double pure_lcp_density_threshold = 0.5;
    using int_type = INT_TYPE;
    syncmer_params cfg;
    int_type seed;
    HASH_FUNCTION hf;
    int_type s_mask;
    int_type k_mask;
    // Spaced-seed compaction: mask char '0' removes its position entirely, so
    // selection runs on the dense (k - m)-mer of kept positions ('1'/'2'/'3'
    // all keep their position; '2'/'3' still bit-mask via k_mask). keep_mask
    // bit i (field at bits 2i) is set when the position is kept; k_eff = k-m.
    uint64_t keep_mask;
    int_type k_eff;
    int64_t table_size;
    bool pure_lcp;
    gtl::vector<bool> table;

    bool canonical;

    constexpr int_type k() const { return cfg.k; }
    constexpr int_type s() const { return cfg.smer.s.value(); }
    constexpr int_type t() const { return cfg.smer.t.value(); }
    constexpr double density() const { return cfg.density.value(); }
    constexpr int_type threshold() const {
        return density() * std::numeric_limits<int_type>::max();
    }
    constexpr auto type() const -> syncmer_params::typeenum {
        return cfg.type();
    }
    static constexpr uint64_t build_keep_mask(int_type km, unsigned k) {
        uint64_t keep = 0;
        for (unsigned i = 0; i < k; ++i) {
            if (((km >> (2 * i)) & 3) != 0) keep |= (1ULL << i);
        }
        return keep;
    }
    static constexpr int_type count_kept(uint64_t keep, unsigned k) {
        int_type c = 0;
        for (unsigned i = 0; i < k; ++i) c += (keep >> i) & 1;
        return c;
    }
    // Gather the kept 2-bit fields of v (ascending field order) into a dense
    // k_eff-mer. Identity when no position is removed, so masks without '0'
    // reproduce the previous AND-mask behavior bit-for-bit.
    constexpr int_type compact(int_type v) const {
        if (keep_mask == ((1ULL << k()) - 1ULL)) return v;
        int_type out = 0;
        unsigned o = 0;
        for (unsigned i = 0; i < k(); ++i) {
            if ((keep_mask >> i) & 1) {
                out |= (int_type)(((v >> (2 * i)) & 3) << (2 * o));
                ++o;
            }
        }
        return out;
    }
    GenericSubmerRule(const syncmer_params &cfg, int_type seed)
        : cfg{cfg},
          seed{seed},
          hf{seed},
          k_mask{(INT_TYPE)cfg.compute_spaced_mask()},
          keep_mask{build_keep_mask(k_mask, cfg.k)},
          k_eff{count_kept(keep_mask, cfg.k)},
          table_size{2 << (2 * k_eff - 1)},
          pure_lcp(false),
          table(table_size, 0),
          canonical(cfg.canonical) {
        if (cfg.smer.s.has_value()) {
            // s bases packed 2 bits each ⇒ the s-mer occupies 2*s low bits.
            s_mask = (int_type)((1ULL << (2 * s())) - 1ULL);
        } else {
            s_mask = -1;
        }
        size_t tc = 0;
        for (int i = 0; i < table_size; ++i) {
            // i IS the compacted k_eff-mer (the LuT lives in compacted space,
            // so it only needs 4^k_eff entries); queries compact into it.
            table[i] = is_submer(i);
            tc += table[i];
        }
        fprintf(stderr, "Real submer density: %lf\n", (double)tc / table_size);
        if ((double)tc / table_size > pure_lcp_density_threshold) {
            fprintf(stderr, "Density of %lf > %lf, switching to pure LCP!\n",
                    (double)tc / table_size, pure_lcp_density_threshold);
            pure_lcp = true;
        }
    }
    // Strand-canonical selection (opt-in): apply the rule to min(κ, rc(κ)) so a
    // k-mer and its reverse complement are always co-selected. Default off ⇒
    // the raw forward value is used (unchanged behavior). The query value is
    // bit-masked and compacted into LuT space; with no '0' mask chars this is
    // `table[val & k_mask]`, bit-identical to the pre-compaction behavior.
    constexpr bool operator()(int_type val) const {
        if (canonical) {
            val = std::min(val, reverse_complement_kmer<int_type>(val, k()));
        }
        return table[compact(val & k_mask)];
    }
    // µ(κ,s) = smallest s-mer (first on ties). Closed syncmer: µ at the first or
    // last s-mer; open syncmer: µ at offset t. Mincode submer: hashed code < H/c.
    // `val` is the COMPACTED k_eff-mer (masked positions already removed by
    // compact()), so window indices and t live in compacted space: there are
    // k_eff-s+1 s-mer windows; window i starts at compacted base i (bit 2*i).
    constexpr bool is_submer(int_type val) const {
        switch (type()) {
            case syncmer_params::typeenum::ClosedSyncmer:  // Intentional
                                                           // fallthrough
            case syncmer_params::typeenum::OpenSyncmer: {
                int_type min_index = 0;
                int_type min_val = std::numeric_limits<int_type>::max();
                for (int_type i = 0; i < k_eff - s() + 1; ++i) {
                    int_type masked_s = (val >> (2 * i)) & s_mask;
                    int_type hashed_s = hf(masked_s);
                    if (hashed_s < min_val) {
                        min_index = i;
                        min_val = hashed_s;
                    }
                }
                if (type() == syncmer_params::typeenum::ClosedSyncmer) {
                    // Closed syncmer: µ at either end (window 0 or the last
                    // window, index k_eff-s).
                    return min_index == 0 || min_index == k_eff - s();
                } else {
                    return t() == min_index;
                }
            }
            case syncmer_params::typeenum::MincodeSubmer: {
                return hf(val) < threshold();
            }
            case syncmer_params::typeenum::Invalid:
            default:
                // Unreachable: syncmer_params::validate() rejects Invalid
                // before any functor is constructed.
                assert(false && "submer type validated before construction");
                __builtin_unreachable();
        }
    }
};

}  // namespace syncmer

#endif
