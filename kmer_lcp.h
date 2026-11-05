// LCP deepening of cores whose level 1 is k-mers, on per-field arrays.
// lcptools parses characters generically (struct lps, lps_deepen). Gallifrey's
// level-1 cores are the selected k-mers, so the k-mer compression and the
// array layout live here, and lcptools supplies only the core label hash.
// core_buf::deepen gives exactly lps_deepen's cores (label, start, end,
// bit_rep, bit_size) for the same level-1 k-mer cores.
#pragma once
#include <sys/mman.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gtl/vector.hpp>
#include <type_traits>
#include <utility>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSE4_1__)
#include <smmintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

#include "ext/lcptools/lcptools_ho.h"

// The cores of a sequence as labels and (start, end) base positions.
struct lcp_read {
    gtl::vector<uint32_t> cores;
    gtl::vector<std::pair<uint32_t, uint32_t>> pos;
    const uint8_t *data = nullptr;
    void reserve(size_t size) {
        cores.reserve(size);
        pos.reserve(size);
    }
    void clear() {
        cores.clear();
        pos.clear();
        data = NULL;
    }
};

// One allocation holding a parse's scratch arrays, 64-byte aligned. Large
// slabs are mapped directly, so releasing one returns its memory to the OS at
// once; from malloc they would come from the arenas, which keep the freed
// space.
struct parse_slab {
    // Below this a slab is read-sized, and a mapping's syscalls and page
    // faults would cost more than the arena space it gives back.
    static constexpr size_t kMapBytes = size_t{1} << 20;
    // A cache line: arrays carved from a slab share no line, and the SIMD
    // kernels get aligned loads.
    static constexpr size_t kAlign = 64;
    static constexpr size_t aligned(size_t b) {
        return (b + kAlign - 1) & ~(kAlign - 1);
    }
    std::byte *p = nullptr;
    size_t bytes = 0;
    parse_slab() = default;
    parse_slab(const parse_slab &) = delete;
    parse_slab &operator=(const parse_slab &) = delete;
    ~parse_slab() { release(); }
    // At least n bytes; the old contents are not kept.
    std::byte *ensure(size_t n) {
        if (n <= bytes) return p;
        release();
        n = aligned(n);
        void *q;
        if (n >= kMapBytes) {
            q = mmap(nullptr, n, PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
            if (q == MAP_FAILED) q = nullptr;
        } else {
            q = aligned_alloc(kAlign, n);
        }
        if (q == nullptr) {
            fprintf(stderr, "parse_slab: out of memory\n");
            abort();
        }
        p = static_cast<std::byte *>(q);
        bytes = n;
        return p;
    }
    void release() {
        if (bytes >= kMapBytes) {
            munmap(p, bytes);
        } else {
            free(p);
        }
        p = nullptr;
        bytes = 0;
    }
};

namespace kmer_lcp {

// bit_size of a compressed core with code r (r < 128).
inline uint32_t code_bits(uint32_t r) {
    if (r < 4) return 2;
    return 32 - static_cast<uint32_t>(__builtin_clz(r));
}

// Compressed codes of k-mer cores: c[j] from cores j (left) and j + 1.
// With x = label_j ^ label_j+1 != 0 the first differing bit is ctz(x) <
// kmer_bits and the code is 2 * ctz(x) + that bit of the right core; equal
// cores give 2 * kmer_bits (bit kmer_bits of 0x8000000000000000 | code is
// 0). In the SIMD paths ctz(x) is the exponent of (float)(x & -x).
inline void compress_kmers(const uint32_t *label, int64_t m, uint32_t kmer_bits,
                           uint8_t *c) {
    int64_t j = 0;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    const uint32x4_t eqcode = vdupq_n_u32(2 * kmer_bits);
    for (; j + 8 <= m; j += 8) {
        uint16x4_t half[2];
        for (int h = 0; h < 2; h++) {
            const uint32x4_t l = vld1q_u32(label + j + 4 * h);
            const uint32x4_t r = vld1q_u32(label + j + 4 * h + 1);
            const uint32x4_t x = veorq_u32(l, r);
            const uint32x4_t low = vandq_u32(
                x, vreinterpretq_u32_s32(vnegq_s32(vreinterpretq_s32_u32(x))));
            const uint32x4_t f = vsubq_u32(
                vshrq_n_u32(vreinterpretq_u32_f32(vcvtq_f32_u32(low)), 23),
                vdupq_n_u32(127));
            const uint32x4_t bit = vshrq_n_u32(vtstq_u32(r, low), 31);
            const uint32x4_t v = vaddq_u32(vaddq_u32(f, f), bit);
            half[h] = vmovn_u32(vbslq_u32(vceqzq_u32(x), eqcode, v));
        }
        vst1_u8(c + j, vmovn_u16(vcombine_u16(half[0], half[1])));
    }
#elif defined(__AVX2__)
    const __m256i c127 = _mm256_set1_epi32(127), cff = _mm256_set1_epi32(0xFF);
    const __m256i eqcode = _mm256_set1_epi32(static_cast<int>(2 * kmer_bits));
    const __m256i one = _mm256_set1_epi32(1), zero = _mm256_setzero_si256();
    for (; j + 8 <= m; j += 8) {
        const __m256i l = _mm256_loadu_si256((const __m256i *)(label + j));
        const __m256i r = _mm256_loadu_si256((const __m256i *)(label + j + 1));
        const __m256i x = _mm256_xor_si256(l, r);
        const __m256i low = _mm256_and_si256(x, _mm256_sub_epi32(zero, x));
        // The float is negative for bit 31: drop the sign with the mask.
        const __m256i f = _mm256_sub_epi32(
            _mm256_and_si256(
                _mm256_srli_epi32(_mm256_castps_si256(_mm256_cvtepi32_ps(low)),
                                  23),
                cff),
            c127);
        const __m256i bit = _mm256_andnot_si256(
            _mm256_cmpeq_epi32(_mm256_and_si256(r, low), zero), one);
        const __m256i v = _mm256_add_epi32(_mm256_add_epi32(f, f), bit);
        const __m256i w =
            _mm256_blendv_epi8(v, eqcode, _mm256_cmpeq_epi32(x, zero));
        const __m128i p16 = _mm_packs_epi32(_mm256_castsi256_si128(w),
                                            _mm256_extracti128_si256(w, 1));
        _mm_storel_epi64((__m128i *)(c + j), _mm_packus_epi16(p16, p16));
    }
#elif defined(__SSE2__)
    const __m128i c127 = _mm_set1_epi32(127), cff = _mm_set1_epi32(0xFF);
    const __m128i eqcode = _mm_set1_epi32(static_cast<int>(2 * kmer_bits));
    const __m128i one = _mm_set1_epi32(1), zero = _mm_setzero_si128();
    for (; j + 8 <= m; j += 8) {
        __m128i w[2];
        for (int h = 0; h < 2; h++) {
            const __m128i l =
                _mm_loadu_si128((const __m128i *)(label + j + 4 * h));
            const __m128i r =
                _mm_loadu_si128((const __m128i *)(label + j + 4 * h + 1));
            const __m128i x = _mm_xor_si128(l, r);
            const __m128i low = _mm_and_si128(x, _mm_sub_epi32(zero, x));
            const __m128i f = _mm_sub_epi32(
                _mm_and_si128(
                    _mm_srli_epi32(_mm_castps_si128(_mm_cvtepi32_ps(low)), 23),
                    cff),
                c127);
            const __m128i bit = _mm_andnot_si128(
                _mm_cmpeq_epi32(_mm_and_si128(r, low), zero), one);
            const __m128i v = _mm_add_epi32(_mm_add_epi32(f, f), bit);
            const __m128i e = _mm_cmpeq_epi32(x, zero);
            w[h] =
                _mm_or_si128(_mm_and_si128(e, eqcode), _mm_andnot_si128(e, v));
        }
        const __m128i p16 = _mm_packs_epi32(w[0], w[1]);
        _mm_storel_epi64((__m128i *)(c + j), _mm_packus_epi16(p16, p16));
    }
#endif
    for (; j < m; j++) {
        const uint32_t x = label[j] ^ label[j + 1];
        if (x == 0) {
            c[j] = static_cast<uint8_t>(2 * kmer_bits);
        } else {
            const int f = __builtin_ctz(x);
            c[j] = static_cast<uint8_t>(2 * f + ((label[j + 1] >> f) & 1));
        }
    }
}

// Compressed codes of general cores (levels >= 2), as lcptools'
// core_compress().
inline void compress(const uint64_t *br, const uint8_t *bs, int64_t m,
                     uint8_t *c) {
    for (int64_t j = 0; j < m; j++) {
        const uint64_t l = br[j], r = br[j + 1];
        const uint32_t lb = bs[j], rb = bs[j + 1];
        const uint64_t x = l ^ r;
        uint32_t f = rb;
        if (x != 0) f = static_cast<uint32_t>(__builtin_ctzll(x));
        f = std::min({f, lb, rb});
        c[j] = static_cast<uint8_t>(2 * f + ((r >> f) & 1));
    }
}

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
// One bit per byte of v[0..3] (each byte holds its bit at its lane shift).
inline uint64_t pack64(const uint8x16_t v[4]) {
    return vgetq_lane_u64(
        vreinterpretq_u64_u8(
            vpaddq_u8(vpaddq_u8(vpaddq_u8(v[0], v[1]), vpaddq_u8(v[2], v[3])),
                      vdupq_n_u8(0))),
        0);
}
#endif

// Pair masks of c (readable up to c[64 * words]): bit i of eq/gt/lt is
// c[i] ==/>/< c[i + 1]. Codes are < 128, so signed byte compares work.
inline void pair_masks(const uint8_t *c, int64_t words, uint64_t *eq,
                       uint64_t *gt, uint64_t *lt) {
    for (int64_t w = 0; w < words; w++) {
        const uint8_t *p = c + 64 * w;
        uint64_t e = 0, g = 0, l = 0;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        static const int8_t shift[16] = {0, 1, 2, 3, 4, 5, 6, 7,
                                         0, 1, 2, 3, 4, 5, 6, 7};
        const int8x16_t sh = vld1q_s8(shift);
        uint8x16_t me[4], mg[4], ml[4];
        for (int q = 0; q < 4; q++) {
            const uint8x16_t a = vld1q_u8(p + 16 * q),
                             b = vld1q_u8(p + 16 * q + 1);
            me[q] = vshlq_u8(vshrq_n_u8(vceqq_u8(a, b), 7), sh);
            mg[q] = vshlq_u8(vshrq_n_u8(vcgtq_u8(a, b), 7), sh);
            ml[q] = vshlq_u8(vshrq_n_u8(vcltq_u8(a, b), 7), sh);
        }
        e = pack64(me);
        g = pack64(mg);
        l = pack64(ml);
#elif defined(__AVX2__)
        for (int h = 0; h < 2; h++) {
            const __m256i a = _mm256_loadu_si256((const __m256i *)(p + 32 * h));
            const __m256i b =
                _mm256_loadu_si256((const __m256i *)(p + 32 * h + 1));
            e |= uint64_t{static_cast<uint32_t>(
                     _mm256_movemask_epi8(_mm256_cmpeq_epi8(a, b)))}
                 << (32 * h);
            g |= uint64_t{static_cast<uint32_t>(
                     _mm256_movemask_epi8(_mm256_cmpgt_epi8(a, b)))}
                 << (32 * h);
            l |= uint64_t{static_cast<uint32_t>(
                     _mm256_movemask_epi8(_mm256_cmpgt_epi8(b, a)))}
                 << (32 * h);
        }
#elif defined(__SSE2__)
        for (int q = 0; q < 4; q++) {
            const __m128i a = _mm_loadu_si128((const __m128i *)(p + 16 * q));
            const __m128i b =
                _mm_loadu_si128((const __m128i *)(p + 16 * q + 1));
            e |= uint64_t{static_cast<uint32_t>(
                     _mm_movemask_epi8(_mm_cmpeq_epi8(a, b)))}
                 << (16 * q);
            g |= uint64_t{static_cast<uint32_t>(
                     _mm_movemask_epi8(_mm_cmpgt_epi8(a, b)))}
                 << (16 * q);
            l |= uint64_t{static_cast<uint32_t>(
                     _mm_movemask_epi8(_mm_cmpgt_epi8(b, a)))}
                 << (16 * q);
        }
#else
        for (int i = 0; i < 64; i++) {
            e |= uint64_t{p[i] == p[i + 1]} << i;
            g |= uint64_t{p[i] > p[i + 1]} << i;
            l |= uint64_t{p[i] < p[i + 1]} << i;
        }
#endif
        eq[w] = e;
        gt[w] = g;
        lt[w] = l;
    }
}

// out[i] = label of the core over compressed [b, e), key[i] = b | e << 32:
// lcptools' MurmurHash3_32(seed 42) of (lab[b + 1], lab[e - 1], lab[e],
// e - b - 2), as init_core3 hashes a core's constituents. Lanes are
// independent, so SIMD hides the multiply latency of one hash.
inline void labels(const uint64_t *key, const uint32_t *lab, uint32_t *out,
                   int64_t n) {
    int64_t i = 0;
#if defined(__AVX2__) && !(defined(__ARM_NEON) || defined(__ARM_NEON__))
    const __m256i c1 = _mm256_set1_epi32(static_cast<int>(0xcc9e2d51u));
    const __m256i c2 = _mm256_set1_epi32(0x1b873593);
    const __m256i c5 = _mm256_set1_epi32(5);
    const __m256i ce = _mm256_set1_epi32(static_cast<int>(0xe6546b64u));
    const __m256i f1 = _mm256_set1_epi32(static_cast<int>(0x85ebca6bu));
    const __m256i f2 = _mm256_set1_epi32(static_cast<int>(0xc2b2ae35u));
    const __m256i one = _mm256_set1_epi32(1), two = _mm256_set1_epi32(2);
    const __m256i split = _mm256_setr_epi32(0, 2, 4, 6, 1, 3, 5, 7);
    const int *L = reinterpret_cast<const int *>(lab);
    for (; i + 8 <= n; i += 8) {
        const __m256i k0 = _mm256_permutevar8x32_epi32(
            _mm256_loadu_si256((const __m256i *)(key + i)), split);
        const __m256i k1 = _mm256_permutevar8x32_epi32(
            _mm256_loadu_si256((const __m256i *)(key + i + 4)), split);
        const __m256i B = _mm256_permute2x128_si256(k0, k1, 0x20);
        const __m256i E = _mm256_permute2x128_si256(k0, k1, 0x31);
        __m256i in[4];
        in[0] = _mm256_i32gather_epi32(L, _mm256_add_epi32(B, one), 4);
        in[1] = _mm256_i32gather_epi32(L, _mm256_sub_epi32(E, one), 4);
        in[2] = _mm256_i32gather_epi32(L, E, 4);
        in[3] = _mm256_sub_epi32(_mm256_sub_epi32(E, B), two);
        __m256i h = _mm256_set1_epi32(42);
        for (int q = 0; q < 4; q++) {
            __m256i k = _mm256_mullo_epi32(in[q], c1);
            k = _mm256_or_si256(_mm256_slli_epi32(k, 15),
                                _mm256_srli_epi32(k, 17));
            k = _mm256_mullo_epi32(k, c2);
            h = _mm256_xor_si256(h, k);
            h = _mm256_or_si256(_mm256_slli_epi32(h, 15),
                                _mm256_srli_epi32(h, 17));
            h = _mm256_add_epi32(_mm256_mullo_epi32(h, c5), ce);
        }
        h = _mm256_xor_si256(h, _mm256_set1_epi32(16));
        h = _mm256_xor_si256(h, _mm256_srli_epi32(h, 16));
        h = _mm256_mullo_epi32(h, f1);
        h = _mm256_xor_si256(h, _mm256_srli_epi32(h, 13));
        h = _mm256_mullo_epi32(h, f2);
        h = _mm256_xor_si256(h, _mm256_srli_epi32(h, 16));
        _mm256_storeu_si256((__m256i *)(out + i), h);
    }
#elif defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__SSE4_1__)
    // 4 lanes; the inputs are loaded one by one (no gather).
    for (; i + 4 <= n; i += 4) {
        uint32_t w[4][4];
        for (int t = 0; t < 4; t++) {
            const uint32_t b = static_cast<uint32_t>(key[i + t]);
            const uint32_t e = static_cast<uint32_t>(key[i + t] >> 32);
            w[0][t] = lab[b + 1];
            w[1][t] = lab[e - 1];
            w[2][t] = lab[e];
            w[3][t] = e - b - 2;
        }
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        uint32x4_t h = vdupq_n_u32(42);
        for (int q = 0; q < 4; q++) {
            uint32x4_t k = vmulq_n_u32(vld1q_u32(w[q]), 0xcc9e2d51u);
            k = vsriq_n_u32(vshlq_n_u32(k, 15), k, 17);
            k = vmulq_n_u32(k, 0x1b873593u);
            h = veorq_u32(h, k);
            h = vsriq_n_u32(vshlq_n_u32(h, 15), h, 17);
            h = vaddq_u32(vmulq_n_u32(h, 5), vdupq_n_u32(0xe6546b64u));
        }
        h = veorq_u32(h, vdupq_n_u32(16));
        h = veorq_u32(h, vshrq_n_u32(h, 16));
        h = vmulq_n_u32(h, 0x85ebca6bu);
        h = veorq_u32(h, vshrq_n_u32(h, 13));
        h = vmulq_n_u32(h, 0xc2b2ae35u);
        h = veorq_u32(h, vshrq_n_u32(h, 16));
        vst1q_u32(out + i, h);
#else
        __m128i h = _mm_set1_epi32(42);
        for (int q = 0; q < 4; q++) {
            __m128i k =
                _mm_mullo_epi32(_mm_loadu_si128((const __m128i *)w[q]),
                                _mm_set1_epi32(static_cast<int>(0xcc9e2d51u)));
            k = _mm_or_si128(_mm_slli_epi32(k, 15), _mm_srli_epi32(k, 17));
            k = _mm_mullo_epi32(k, _mm_set1_epi32(0x1b873593));
            h = _mm_xor_si128(h, k);
            h = _mm_or_si128(_mm_slli_epi32(h, 15), _mm_srli_epi32(h, 17));
            h = _mm_add_epi32(_mm_mullo_epi32(h, _mm_set1_epi32(5)),
                              _mm_set1_epi32(static_cast<int>(0xe6546b64u)));
        }
        h = _mm_xor_si128(h, _mm_set1_epi32(16));
        h = _mm_xor_si128(h, _mm_srli_epi32(h, 16));
        h = _mm_mullo_epi32(h, _mm_set1_epi32(static_cast<int>(0x85ebca6bu)));
        h = _mm_xor_si128(h, _mm_srli_epi32(h, 13));
        h = _mm_mullo_epi32(h, _mm_set1_epi32(static_cast<int>(0xc2b2ae35u)));
        h = _mm_xor_si128(h, _mm_srli_epi32(h, 16));
        _mm_storeu_si128((__m128i *)(out + i), h);
#endif
    }
#endif
    for (; i < n; i++) {
        const uint32_t b = static_cast<uint32_t>(key[i]);
        const uint32_t e = static_cast<uint32_t>(key[i] >> 32);
        const uint32_t in[4] = {lab[b + 1], lab[e - 1], lab[e], e - b - 2};
        out[i] = MurmurHash3_32(in, sizeof in, 42);
    }
}

// Bits i + s of mask v at word w (s in {-1, 1, 2}); v[words] is 0.
inline uint64_t at(const uint64_t *v, int64_t w, int s) {
    if (s > 0) return (v[w] >> s) | (v[w + 1] << (64 - s));
    uint64_t carry = 0;
    if (w > 0) carry = v[w - 1] >> 63;
    return (v[w] << 1) | carry;
}

// Bits of word w (positions 64w..64w+63) below position lim.
inline uint64_t below(int64_t w, int64_t lim) {
    const int64_t lo = 64 * w;
    if (lim <= lo) return 0;
    if (lim - lo >= 64) return ~uint64_t{0};
    return (uint64_t{1} << (lim - lo)) - 1;
}

}  // namespace kmer_lcp

// The cores of one LCP level, deepened in place. Level-1 cores are the
// selected k-mers: label = k-mer code (kmer_bits wide), span = start |
// end << 32 in base coordinates; their bit_rep (0x8000000000000000 | label)
// and bit_size (kmer_bits) are implicit, so only labels and spans are
// stored. All arrays are carved from one slab, and the run-length stream of
// a parse lives in a second one.
struct core_buf {
    int level = 1;
    int kmer_bits = 0;  // label width of the level-1 cores
    int64_t size = 0;
    uint32_t *label = nullptr;
    uint64_t *span = nullptr;
    uint64_t *bit_rep = nullptr;  // levels >= 2
    uint8_t *bit_size = nullptr;  // levels >= 2

    core_buf() = default;
    core_buf(const core_buf &) = delete;
    core_buf &operator=(const core_buf &) = delete;
    uint32_t start(int64_t i) const { return static_cast<uint32_t>(span[i]); }
    uint32_t end(int64_t i) const {
        return static_cast<uint32_t>(span[i] >> 32);
    }
    // The previous level's index of core i's first constituent, after a
    // deepening step: key[i] = b | e << 32 holds core i's constituents as
    // compressed cores [b, e), and compressed core b is the previous level's
    // core b + 1 compressed against core b.
    int64_t first(int64_t i) const {
        return static_cast<uint32_t>(key[i]) + int64_t{1};
    }
    // Room for n cores, dropping the current ones (called before a parse
    // writes its level-1 cores); grows by at least 1.5x so per-read calls
    // stay rare.
    void reserve(int64_t n) {
        if (n <= capacity) return;
        const int64_t c = std::max(n, capacity + capacity / 2);
        const size_t uc = static_cast<size_t>(c);
        // The code array has 64 bytes of padding for the pair-mask loads, and
        // the pair masks hold a bit per core plus two words.
        auto arrays = [&](auto &&f) {
            f(label, uc);
            f(span, uc);
            f(bit_rep, uc);
            f(bit_size, uc);
            f(label_next, uc);
            f(span_next, uc);
            f(bit_rep_next, uc);
            f(bit_size_next, uc);
            f(code, uc + 64);
            f(key, uc);
            f(eq, uc / 64 + 2);
            f(gt, uc / 64 + 2);
            f(lt, uc / 64 + 2);
        };
        auto bytes = [](auto *ptr, size_t len) {
            return parse_slab::aligned(len * sizeof(*ptr));
        };
        size_t total = 0;
        arrays([&](auto *&ptr, size_t len) { total += bytes(ptr, len); });
        std::byte *m = cores_mem.ensure(total);
        arrays([&](auto *&ptr, size_t len) {
            ptr = reinterpret_cast<std::remove_reference_t<decltype(ptr)>>(m);
            m += bytes(ptr, len);
        });
        capacity = c;
    }
    // Run-length scratch for n bases: {bases, positions}.
    std::pair<uint8_t *, uint32_t *> rle_scratch(size_t n) {
        const size_t pos_bytes = parse_slab::aligned(n * sizeof(uint32_t));
        std::byte *m = rle_mem.ensure(pos_bytes + n);
        return {reinterpret_cast<uint8_t *>(m + pos_bytes),
                reinterpret_cast<uint32_t *>(m)};
    }
    // Deepens to lcp_level, like lcptools' lps_deepen.
    void deepen(int lcp_level) {
        while (level < lcp_level) {
            // As lcp_dct: below two cores a level leaves none, and the parse
            // ends.
            if (size < 2) {
                size = 0;
                ++level;
                return;
            }
            size = next_level();
            std::swap(label, label_next);
            std::swap(span, span_next);
            std::swap(bit_rep, bit_rep_next);
            std::swap(bit_size, bit_size_next);
            ++level;
        }
    }
    // Drops the first level-1 core; level-1 bit_rep and bit_size are implicit
    // in the label, so only labels and spans move.
    void drop_first() {
        --size;
        std::memmove(label, label + 1, size * sizeof(*label));
        std::memmove(span, span + 1, size * sizeof(*span));
    }
    void append_to(lcp_read &lr, int64_t b, int64_t e) const {
        lr.cores.insert(lr.cores.end(), label + b, label + e);
        for (int64_t i = b; i < e; ++i) lr.pos.emplace_back(start(i), end(i));
    }
    void append_to(lcp_read &lr) const { append_to(lr, 0, size); }

   private:
    int64_t capacity = 0;
    // The next level, swapped in by deepen.
    uint32_t *label_next = nullptr;
    uint64_t *span_next = nullptr;
    uint64_t *bit_rep_next = nullptr;
    uint8_t *bit_size_next = nullptr;
    uint8_t *code = nullptr;  // compressed cores
    uint64_t *key = nullptr;  // constituents of each next-level core
    uint64_t *eq = nullptr, *gt = nullptr, *lt = nullptr;  // pair masks
    parse_slab cores_mem;
    parse_slab rle_mem;

    // One level: size >= 2 cores -> the next level's cores in the *_next
    // arrays; returns their count. Compressed core j is core j + 1
    // compressed against core j: its label and end are core j + 1's, its
    // start is core j's. The rules follow lcptools' parse3().
    int64_t next_level() {
        using namespace kmer_lcp;
        const int64_t m = size - 1;
        uint8_t *c = code;
        if (level == 1) {
            compress_kmers(label, m, static_cast<uint32_t>(kmer_bits), c);
        } else {
            compress(bit_rep, bit_size, m, c);
        }
        std::memset(c + m, 0, 64);

        const int64_t words = (m + 63) / 64;
        uint64_t *E = eq, *G = gt, *Lt = lt;
        pair_masks(c, words, E, G, Lt);
        E[words] = G[words] = Lt[words] = 0;
        // Pairs (i, i + 1) exist for i <= m - 2.
        const int64_t npairs = std::max<int64_t>(m - 1, 0);
        for (int64_t w = npairs / 64; w < words; w++) {
            const uint64_t keep = below(w, npairs);
            E[w] &= keep;
            G[w] &= keep;
            Lt[w] &= keep;
        }

        int64_t k = 0;
        // New core k over compressed [b, e); bit_rep keeps whole codes from
        // the right while they fit in 64 bits (at most 21 bits for 3 codes).
        auto put = [&](int64_t b, int64_t e) {
            span_next[k] = (span[b] & 0xFFFFFFFFu) | (span[e] & ~0xFFFFFFFFull);
            key[k] = static_cast<uint64_t>(b) | static_cast<uint64_t>(e) << 32;
            uint32_t bits = 0;
            uint64_t rep = 0;
            if (e - b == 3) {
                const uint32_t s0 = code_bits(c[b]), s1 = code_bits(c[b + 1]),
                               s2 = code_bits(c[b + 2]);
                bits = s0 + s1 + s2;
                rep = c[b + 2] | uint64_t{c[b + 1]} << s2 |
                      uint64_t{c[b]} << (s2 + s1);
            } else {
                uint32_t index = 0;
                bool full = false;
                for (int64_t t = e - 1; b <= t; t--) {
                    const uint32_t sz = code_bits(c[t]);
                    bits += sz;
                    if (!full && index + sz <= 64) {
                        rep |= uint64_t{c[t]} << index;
                        index += sz;
                    } else {
                        full = true;
                    }
                }
            }
            bit_rep_next[k] = rep & 0x7FFFFFFFFFFFFFFFull;
            bit_size_next[k] =
                static_cast<uint8_t>(std::min<uint32_t>(bits, 63));
            k++;
        };

        int64_t i2 = m;  // end of the last core; m: none yet
        for (int64_t w = 0; w < words; w++) {
            const uint64_t e0 = E[w], e1 = at(E, w, 1);
            const uint64_t g0 = G[w], g1 = at(G, w, 1), gm = at(G, w, -1);
            const uint64_t l0 = Lt[w], l1 = at(Lt, w, 1), l2 = at(Lt, w, 2);
            const uint64_t rint = ~e0 & e1;
            const uint64_t lmin = g0 & l1;
            uint64_t lmax = l0 & g1 & ~gm & ~l2;
            if (w == 0) lmax &= ~uint64_t{1};
            // Every core needs i + 2 < m, LMAX i + 3 < m.
            uint64_t starts =
                ((rint | lmin) & below(w, m - 2)) | (lmax & below(w, m - 3));
            while (starts) {
                const int64_t i = 64 * w + __builtin_ctzll(starts);
                starts &= starts - 1;
                // The core ends one past the first unequal pair after i + 1:
                // the end of the run for RINT, pair i + 1 itself for
                // LMIN/LMAX.
                int64_t pw = (i + 1) / 64;
                uint64_t v = ~E[pw] & (~uint64_t{0} << ((i + 1) % 64));
                while (v == 0) v = ~E[++pw];
                const int64_t t = 64 * pw + __builtin_ctzll(v) + 1;
                if (t >= m) continue;  // a run reaching the end is no core
                if (i2 < i) put(i2 - 1, i + 1);  // the gap since the last core
                i2 = t + 1;
                put(i, i2);
            }
        }
        labels(key, label, label_next, k);
        return k;
    }
};
