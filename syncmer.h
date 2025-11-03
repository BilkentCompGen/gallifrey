#ifndef _SYNCMER_H_
#define _SYNCMER_H_

#include <concepts>
#include <limits>

namespace syncmer {
template <class Rule>
concept SyncmerRule = requires(Rule k) {
    { k(unsigned{}) } -> std::same_as<bool>;
};
template <class HF>
concept HashFunction = requires(HF h) {
    { HF{unsigned{}} } -> std::same_as<HF>;
    { h(unsigned{}) } -> std::same_as<unsigned>;
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
                break;
            case 2:
                k1 ^= tail[1] << 8;
                break;
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
struct Identity{
    unsigned seed;
    Identity(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const {
        return value;
    }
};
template <class INT_TYPE, HashFunction HASH_FUNCTION>
struct DensitySyncmerRule {
    double density;
    double threshold;
    INT_TYPE seed;
    HASH_FUNCTION hf;
    DensitySyncmerRule(double density, INT_TYPE seed)
        : density{density},
          threshold{std::numeric_limits<INT_TYPE>::max() * density},
          seed{seed},
          hf{seed} {}
    constexpr bool operator()(INT_TYPE val) const {
        return hf(val) < threshold;
    }
};
}  // namespace syncmer

#endif
