#ifndef _SYNCMER_H_
#define _SYNCMER_H_

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <limits>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSSE3__)
#include <tmmintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

#include "asm.h"
#include "gtl/vector.hpp"
namespace syncmer {
template <class HF, class INT_TYPE>
concept HashFunction = requires(HF h) {
    { HF{INT_TYPE{}} } -> std::same_as<HF>;
    { h(INT_TYPE{}) } -> std::same_as<INT_TYPE>;
};

// MurmurHash3_32 of the 4-byte key `value` (this implementation's body step
// rotates h1 by 15), written branch- and memory-free so batch callers
// vectorize it.
struct MurMur {
    unsigned seed;
    MurMur(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const {
        uint32_t k1 = value * 0xcc9e2d51u;
        k1 = (k1 << 15) | (k1 >> 17);
        k1 *= 0x1b873593u;
        uint32_t h1 = seed ^ k1;
        h1 = (h1 << 15) | (h1 >> 17);
        h1 = h1 * 5 + 0xe6546b64u;
        h1 ^= 4u;  // key length
        h1 ^= h1 >> 16;
        h1 *= 0x85ebca6bu;
        h1 ^= h1 >> 13;
        h1 *= 0xc2b2ae35u;
        h1 ^= h1 >> 16;
        return h1;
    }
};
struct Identity {
    unsigned seed;
    Identity(unsigned seed) : seed{seed} {}
    unsigned operator()(unsigned value) const { return value; }
};

// Per-thread scratch for GenericSubmerRule::select_live, reused across reads.
struct live_scratch {
    // Leading pad of `q` so q[p - 12] is addressable for every p >= 0; the
    // pad stays zero (only [0, n) is ever written).
    static constexpr int64_t kPad = 16;
    // q[p]: the 4 symbols ending at p, 2 bits each, oldest in the high bits
    // (symbols before the stream start read as 0). Any s-mer (s <= 16) or
    // k-mer (k <= 16) ending at p is an OR of q[p], q[p-4], q[p-8], q[p-12].
    gtl::vector<uint8_t> q;
    // Per-position s-mer order key; only the one matching s is used.
    gtl::vector<uint8_t> r8;    // rank, s <= 4
    gtl::vector<uint16_t> r16;  // rank, s <= 8
    gtl::vector<uint32_t> h32;  // hash, s > 8
    gtl::vector<uint8_t> sel;
    void require(size_t n) {
        if (sel.size() < n) {
            q.resize(n + kPad);
            sel.resize(n);
        }
    }
    uint8_t *codes() { return q.data() + kPad; }
};

// r[p] = rank[q[p] & m] for p in [p0, n), with a 256-entry byte table (the
// s <= 4 s-mer rank LuT, zero-padded). A scalar gather costs ~2-3 cycles per
// element; the table fits in registers, so SIMD looks up 16/32 per step.
inline void lookup_rank8(const uint8_t *__restrict q, uint8_t *__restrict r,
                         int64_t p0, int64_t n, uint8_t m,
                         const uint8_t *__restrict rank) {
    int64_t p = p0;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    // Four 64-byte tbl lookups; after the first, indices are shifted down by
    // 64 and tbx leaves lanes whose index is out of range (>= 64) unchanged.
    const uint8x16x4_t t0 = vld1q_u8_x4(rank), t1 = vld1q_u8_x4(rank + 64),
                       t2 = vld1q_u8_x4(rank + 128),
                       t3 = vld1q_u8_x4(rank + 192);
    const uint8x16_t vm = vdupq_n_u8(m), c64 = vdupq_n_u8(64);
    for (; p + 16 <= n; p += 16) {
        uint8x16_t i = vandq_u8(vld1q_u8(q + p), vm);
        uint8x16_t v = vqtbl4q_u8(t0, i);
        i = vsubq_u8(i, c64);
        v = vqtbx4q_u8(v, t1, i);
        i = vsubq_u8(i, c64);
        v = vqtbx4q_u8(v, t2, i);
        i = vsubq_u8(i, c64);
        v = vqtbx4q_u8(v, t3, i);
        vst1q_u8(r + p, v);
    }
#elif defined(__AVX2__)
    // Row hi = index >> 4 selects a 16-byte table row; pshufb looks up the
    // low nibble in every row and the row compare keeps the matching one.
    const __m256i vm = _mm256_set1_epi8(static_cast<char>(m));
    const __m256i nib = _mm256_set1_epi8(0x0F);
    for (; p + 32 <= n; p += 32) {
        const __m256i i = _mm256_and_si256(
            _mm256_loadu_si256(reinterpret_cast<const __m256i *>(q + p)), vm);
        const __m256i lo = _mm256_and_si256(i, nib);
        const __m256i hi = _mm256_and_si256(_mm256_srli_epi16(i, 4), nib);
        __m256i v = _mm256_setzero_si256();
        for (int row = 0; row < 16; ++row) {
            const __m256i t = _mm256_broadcastsi128_si256(_mm_loadu_si128(
                reinterpret_cast<const __m128i *>(rank + 16 * row)));
            const __m256i hit =
                _mm256_cmpeq_epi8(hi, _mm256_set1_epi8(static_cast<char>(row)));
            v = _mm256_or_si256(
                v, _mm256_and_si256(hit, _mm256_shuffle_epi8(t, lo)));
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(r + p), v);
    }
#elif defined(__SSSE3__)
    const __m128i vm = _mm_set1_epi8(static_cast<char>(m));
    const __m128i nib = _mm_set1_epi8(0x0F);
    for (; p + 16 <= n; p += 16) {
        const __m128i i = _mm_and_si128(
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(q + p)), vm);
        const __m128i lo = _mm_and_si128(i, nib);
        const __m128i hi = _mm_and_si128(_mm_srli_epi16(i, 4), nib);
        __m128i v = _mm_setzero_si128();
        for (int row = 0; row < 16; ++row) {
            const __m128i t = _mm_loadu_si128(
                reinterpret_cast<const __m128i *>(rank + 16 * row));
            const __m128i hit =
                _mm_cmpeq_epi8(hi, _mm_set1_epi8(static_cast<char>(row)));
            v = _mm_or_si128(v, _mm_and_si128(hit, _mm_shuffle_epi8(t, lo)));
        }
        _mm_storeu_si128(reinterpret_cast<__m128i *>(r + p), v);
    }
#endif
    for (; p < n; ++p) r[p] = rank[q[p] & m];
}

// Code of the m symbols (m <= 16) ending at p, from 4-symbol codes.
inline uint32_t live_code(const uint8_t *q, int64_t p, uint32_t mask) {
    return ((uint32_t(q[p - 12]) << 24) | (uint32_t(q[p - 8]) << 16) |
            (uint32_t(q[p - 4]) << 8) | q[p]) &
           mask;
}

// Per-ISA vector operations on window keys of type R (the primary template
// has none: lanes == 0 selects the scalar loop). Masks are all-ones lanes;
// all() is both the all-true mask and the largest key (min identity).
// store_sel writes one 0/1 byte per lane.
template <class R>
struct window_simd {
    static constexpr int64_t lanes = 0;
};
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
template <>
struct window_simd<uint8_t> {
    using V = uint8x16_t;
    static constexpr int64_t lanes = 16;
    static V load(const uint8_t *p) { return vld1q_u8(p); }
    static V all() { return vdupq_n_u8(0xFF); }
    static V vmin(V a, V b) { return vminq_u8(a, b); }
    static V le(V a, V b) { return vcleq_u8(a, b); }
    static V lt(V a, V b) { return vcltq_u8(a, b); }
    static V vand(V a, V b) { return vandq_u8(a, b); }
    static V vor(V a, V b) { return vorrq_u8(a, b); }
    static void store_sel(uint8_t *s, V m) {
        vst1q_u8(s, vandq_u8(m, vdupq_n_u8(1)));
    }
};
template <>
struct window_simd<uint16_t> {
    using V = uint16x8_t;
    static constexpr int64_t lanes = 8;
    static V load(const uint16_t *p) { return vld1q_u16(p); }
    static V all() { return vdupq_n_u16(0xFFFF); }
    static V vmin(V a, V b) { return vminq_u16(a, b); }
    static V le(V a, V b) { return vcleq_u16(a, b); }
    static V lt(V a, V b) { return vcltq_u16(a, b); }
    static V vand(V a, V b) { return vandq_u16(a, b); }
    static V vor(V a, V b) { return vorrq_u16(a, b); }
    static void store_sel(uint8_t *s, V m) {
        vst1_u8(s, vand_u8(vmovn_u16(m), vdup_n_u8(1)));
    }
};
template <>
struct window_simd<uint32_t> {
    using V = uint32x4_t;
    static constexpr int64_t lanes = 4;
    static V load(const uint32_t *p) { return vld1q_u32(p); }
    static V all() { return vdupq_n_u32(0xFFFFFFFFu); }
    static V vmin(V a, V b) { return vminq_u32(a, b); }
    static V le(V a, V b) { return vcleq_u32(a, b); }
    static V lt(V a, V b) { return vcltq_u32(a, b); }
    static V vand(V a, V b) { return vandq_u32(a, b); }
    static V vor(V a, V b) { return vorrq_u32(a, b); }
    static void store_sel(uint8_t *s, V m) {
        const uint16x4_t h = vmovn_u32(m);
        const uint8x8_t b =
            vand_u8(vmovn_u16(vcombine_u16(h, h)), vdup_n_u8(1));
        const uint32_t x = vget_lane_u32(vreinterpret_u32_u8(b), 0);
        std::memcpy(s, &x, 4);
    }
};
#elif defined(__AVX2__)
// Unsigned compares from min: a <= b ⇔ min(a, b) == a.
template <>
struct window_simd<uint8_t> {
    using V = __m256i;
    static constexpr int64_t lanes = 32;
    static V load(const uint8_t *p) {
        return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
    }
    static V all() { return _mm256_set1_epi8(-1); }
    static V vmin(V a, V b) { return _mm256_min_epu8(a, b); }
    static V le(V a, V b) {
        return _mm256_cmpeq_epi8(_mm256_min_epu8(a, b), a);
    }
    static V lt(V a, V b) { return _mm256_andnot_si256(le(b, a), all()); }
    static V vand(V a, V b) { return _mm256_and_si256(a, b); }
    static V vor(V a, V b) { return _mm256_or_si256(a, b); }
    static void store_sel(uint8_t *s, V m) {
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(s),
                            _mm256_and_si256(m, _mm256_set1_epi8(1)));
    }
};
template <>
struct window_simd<uint16_t> {
    using V = __m256i;
    static constexpr int64_t lanes = 16;
    static V load(const uint16_t *p) {
        return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
    }
    static V all() { return _mm256_set1_epi16(-1); }
    static V vmin(V a, V b) { return _mm256_min_epu16(a, b); }
    static V le(V a, V b) {
        return _mm256_cmpeq_epi16(_mm256_min_epu16(a, b), a);
    }
    static V lt(V a, V b) { return _mm256_andnot_si256(le(b, a), all()); }
    static V vand(V a, V b) { return _mm256_and_si256(a, b); }
    static V vor(V a, V b) { return _mm256_or_si256(a, b); }
    static void store_sel(uint8_t *s, V m) {
        // packs works per 128-bit lane; qwords 0 and 2 hold lanes 0-7, 8-15.
        const V p = _mm256_permute4x64_epi64(_mm256_packs_epi16(m, m), 0x08);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(s),
            _mm_and_si128(_mm256_castsi256_si128(p), _mm_set1_epi8(1)));
    }
};
template <>
struct window_simd<uint32_t> {
    using V = __m256i;
    static constexpr int64_t lanes = 8;
    static V load(const uint32_t *p) {
        return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
    }
    static V all() { return _mm256_set1_epi32(-1); }
    static V vmin(V a, V b) { return _mm256_min_epu32(a, b); }
    static V le(V a, V b) {
        return _mm256_cmpeq_epi32(_mm256_min_epu32(a, b), a);
    }
    static V lt(V a, V b) { return _mm256_andnot_si256(le(b, a), all()); }
    static V vand(V a, V b) { return _mm256_and_si256(a, b); }
    static V vor(V a, V b) { return _mm256_or_si256(a, b); }
    static void store_sel(uint8_t *s, V m) {
        // Per 128-bit lane: 4 dwords -> 4 bytes in the low dword.
        const V p = _mm256_packs_epi16(_mm256_packs_epi32(m, m), m);
        const uint32_t lo = _mm_cvtsi128_si32(_mm256_castsi256_si128(p));
        const uint32_t hi = _mm_cvtsi128_si32(_mm256_extracti128_si256(p, 1));
        const uint64_t x = (lo | (uint64_t(hi) << 32)) & 0x0101010101010101ULL;
        std::memcpy(s, &x, 8);
    }
};
#elif defined(__SSE2__)
// SSE2 has unsigned byte min only; 16/32-bit keys compare as signed after
// flipping the sign bit.
template <>
struct window_simd<uint8_t> {
    using V = __m128i;
    static constexpr int64_t lanes = 16;
    static V load(const uint8_t *p) {
        return _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
    }
    static V all() { return _mm_set1_epi8(-1); }
    static V vmin(V a, V b) { return _mm_min_epu8(a, b); }
    static V le(V a, V b) { return _mm_cmpeq_epi8(_mm_min_epu8(a, b), a); }
    static V lt(V a, V b) { return _mm_andnot_si128(le(b, a), all()); }
    static V vand(V a, V b) { return _mm_and_si128(a, b); }
    static V vor(V a, V b) { return _mm_or_si128(a, b); }
    static void store_sel(uint8_t *s, V m) {
        _mm_storeu_si128(reinterpret_cast<__m128i *>(s),
                         _mm_and_si128(m, _mm_set1_epi8(1)));
    }
};
template <>
struct window_simd<uint16_t> {
    using V = __m128i;
    static constexpr int64_t lanes = 8;
    static V flip(V a) { return _mm_xor_si128(a, _mm_set1_epi16(-32768)); }
    static V load(const uint16_t *p) {
        return _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
    }
    static V all() { return _mm_set1_epi16(-1); }
    static V vmin(V a, V b) { return flip(_mm_min_epi16(flip(a), flip(b))); }
    static V le(V a, V b) {
        return _mm_andnot_si128(_mm_cmpgt_epi16(flip(a), flip(b)), all());
    }
    static V lt(V a, V b) { return _mm_cmpgt_epi16(flip(b), flip(a)); }
    static V vand(V a, V b) { return _mm_and_si128(a, b); }
    static V vor(V a, V b) { return _mm_or_si128(a, b); }
    static void store_sel(uint8_t *s, V m) {
        _mm_storel_epi64(
            reinterpret_cast<__m128i *>(s),
            _mm_and_si128(_mm_packs_epi16(m, m), _mm_set1_epi8(1)));
    }
};
template <>
struct window_simd<uint32_t> {
    using V = __m128i;
    static constexpr int64_t lanes = 4;
    static V flip(V a) { return _mm_xor_si128(a, _mm_set1_epi32(INT32_MIN)); }
    static V load(const uint32_t *p) {
        return _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
    }
    static V all() { return _mm_set1_epi32(-1); }
    static V vmin(V a, V b) {
        const V gt = _mm_cmpgt_epi32(flip(a), flip(b));
        return _mm_or_si128(_mm_and_si128(gt, b), _mm_andnot_si128(gt, a));
    }
    static V le(V a, V b) {
        return _mm_andnot_si128(_mm_cmpgt_epi32(flip(a), flip(b)), all());
    }
    static V lt(V a, V b) { return _mm_cmpgt_epi32(flip(b), flip(a)); }
    static V vand(V a, V b) { return _mm_and_si128(a, b); }
    static V vor(V a, V b) { return _mm_or_si128(a, b); }
    static void store_sel(uint8_t *s, V m) {
        const V p = _mm_packs_epi16(_mm_packs_epi32(m, m), m);
        const uint32_t x =
            _mm_cvtsi128_si32(_mm_and_si128(p, _mm_set1_epi8(1)));
        std::memcpy(s, &x, 4);
    }
};
#endif

// Window selection over per-position s-mer keys h (rank or hash of the s-mer
// ending at p; comparing keys == comparing hashes). For each e in [e0, n),
// window j of the k-mer ending at e is h[e - j] (0 <= j < w), and µ is the
// smallest j among equal minimal keys; sel[e] = 1 iff selected.
//   open (offset t):  µ == t ⇔ h[t] < every newer window (j < t) and
//                     h[t] <= every older window (j > t)
//   closed:           µ == 0 ⇔ h[0] <= all others;
//                     µ == w-1 ⇔ h[w-1] < all others
// The SIMD loop evaluates a vector of consecutive k-mers per window sweep,
// so w and t stay runtime values (templating the kernel on (w, t) was as
// fast but doubled asm2.cpp compile time); two vectors share each sweep to
// halve the short-loop overhead. The scalar loop finishes the tail.
template <class R>
void select_windows(bool open, int w, int t, const R *__restrict h,
                    uint8_t *__restrict sel, int64_t e0, int64_t n) {
    int64_t e = e0;
    // Windows folded into `lo`: [a0, a1) (open j < t, closed 1..w-2); into
    // `hi`: [b0, b1) (open j > t). The pivot windows are compared after.
    const int a0 = open ? 0 : 1, a1 = open ? t : w - 1;
    const int b0 = open ? t + 1 : 0, b1 = open ? w : 0;
    using S = window_simd<R>;
    if constexpr (S::lanes > 0) {
        using V = typename S::V;
        constexpr int64_t L = S::lanes;
        auto finish = [&](int64_t x, V lo, V hi) {
            if (open) {
                const V q = S::load(h + x - t);
                S::store_sel(sel + x,
                             t == 0 ? S::le(q, hi)
                                    : S::vand(S::lt(q, lo), S::le(q, hi)));
            } else {
                const V f = S::load(h + x), l = S::load(h + x - (w - 1));
                S::store_sel(sel + x, S::vor(S::le(f, S::vmin(lo, l)),
                                             S::lt(l, S::vmin(lo, f))));
            }
        };
        for (; e + 2 * L <= n; e += 2 * L) {
            V lo0 = S::all(), hi0 = S::all(), lo1 = S::all(), hi1 = S::all();
            for (int j = a0; j < a1; ++j) {
                lo0 = S::vmin(lo0, S::load(h + e - j));
                lo1 = S::vmin(lo1, S::load(h + e + L - j));
            }
            for (int j = b0; j < b1; ++j) {
                hi0 = S::vmin(hi0, S::load(h + e - j));
                hi1 = S::vmin(hi1, S::load(h + e + L - j));
            }
            finish(e, lo0, hi0);
            finish(e + L, lo1, hi1);
        }
        for (; e + L <= n; e += L) {
            V lo = S::all(), hi = S::all();
            for (int j = a0; j < a1; ++j) lo = S::vmin(lo, S::load(h + e - j));
            for (int j = b0; j < b1; ++j) hi = S::vmin(hi, S::load(h + e - j));
            finish(e, lo, hi);
        }
    }
    constexpr R kMax = std::numeric_limits<R>::max();
    for (; e < n; ++e) {
        R lo = kMax, hi = kMax;
        for (int j = a0; j < a1; ++j) lo = std::min<R>(lo, h[e - j]);
        for (int j = b0; j < b1; ++j) hi = std::min<R>(hi, h[e - j]);
        if (open) {
            const R q = h[e - t];
            sel[e] = (t == 0 || q < lo) & (q <= hi);
        } else {
            const R f = h[e], l = h[e - (w - 1)];
            const bool first = f <= std::min(lo, l);
            const bool last = l < std::min(lo, f);
            sel[e] = first | last;
        }
    }
}

// Open/closed syncmer rule (Edgar 2021), evaluated live: each s-mer gets an
// order key (the dense rank of hf(s-mer) for s <= 8, the hash itself above)
// and select_live() decides every window of a read at once. Nothing is
// tabulated per k-mer. k <= syncmer_params::max_k (32-bit k-mer codes).
// Precondition: cfg.validate() succeeded.
template <class INT_TYPE, HashFunction<INT_TYPE> HASH_FUNCTION>
struct GenericSubmerRule {
    using int_type = INT_TYPE;
    syncmer_params cfg;
    int_type seed;
    HASH_FUNCTION hf;
    // s bases packed 2 bits each: the s-mer occupies the low 2*s bits.
    int_type s_mask;
    // s <= 8: dense rank of hf(s-mer) for every s-mer code. Ranks compare
    // exactly like the hashes (ties included) in 1-2 bytes.
    gtl::vector<uint8_t> smer_rank8;    // s <= 4, zero-padded to 256 entries
    gtl::vector<uint16_t> smer_rank16;  // 4 < s <= 8

    constexpr int_type k() const { return cfg.k; }
    constexpr int_type s() const { return cfg.smer.s.value(); }
    constexpr int_type t() const { return cfg.smer.t.value(); }
    constexpr auto type() const -> syncmer_params::typeenum {
        return cfg.type();
    }
    GenericSubmerRule(const syncmer_params &cfg, int_type seed)
        : cfg{cfg},
          seed{seed},
          hf{seed},
          s_mask{(int_type)((1ULL << (2 * cfg.smer.s.value())) - 1ULL)} {
        assert(cfg.validate() && "syncmer_params validated before use");
        if (s() <= 8) build_smer_ranks();
    }
    void build_smer_ranks() {
        const size_t count = size_t{1} << (2 * s());
        gtl::vector<int_type> hashes(count);
        for (size_t v = 0; v < count; ++v) hashes[v] = hf((int_type)v);
        gtl::vector<int_type> order(hashes.begin(), hashes.end());
        std::sort(order.begin(), order.end());
        order.erase(std::unique(order.begin(), order.end()), order.end());
        auto rank = [&](size_t v) {
            return std::lower_bound(order.begin(), order.end(), hashes[v]) -
                   order.begin();
        };
        if (s() <= 4) {
            smer_rank8.assign(256, 0);  // lookup_rank8 reads a full byte table
            for (size_t v = 0; v < count; ++v) smer_rank8[v] = rank(v);
        } else {
            smer_rank16.resize(count);
            for (size_t v = 0; v < count; ++v) smer_rank16[v] = rank(v);
        }
    }
    // Selection over a run-length symbol stream (2-bit codes): fills
    // S.codes() (see live_scratch::q) for [0, n) and sets S.sel[e] = 1 iff the
    // k-mer ending at symbol e (symbols [e-k+1, e]) is selected, for e in
    // [k-1, n); other sel entries are left untouched. Window j of that k-mer
    // is the s-mer ending at e-j (j = 0 is the lowest bits of the k-mer code),
    // so µ is the smallest j among equal minimal hashes. Closed syncmer: µ at
    // the first or last window; open syncmer: µ at t.
    // Precondition: S.require(n) or larger.
    void select_live(const uint8_t *rle, int64_t n, live_scratch &S) const {
        const int64_t K = k(), Sz = s();
        uint8_t *__restrict q = S.codes();
        for (int64_t p = 0; p < std::min<int64_t>(n, 3); ++p) {
            uint8_t c = 0;
            for (int64_t d = p - 3; d <= p; ++d)
                c = (c << 2) | (d < 0 ? 0 : rle[d]);
            q[p] = c;
        }
        for (int64_t p = 3; p < n; ++p) {
            q[p] = (rle[p - 3] << 6) | (rle[p - 2] << 4) | (rle[p - 1] << 2) |
                   rle[p];
        }
        const bool open = type() == syncmer_params::typeenum::OpenSyncmer;
        const int w = static_cast<int>(K - Sz + 1);
        const int tt = open ? static_cast<int>(t()) : 0;
        const uint32_t smask = s_mask;
        uint8_t *sel = S.sel.data();
        if (Sz <= 4) {
            if (S.r8.size() < static_cast<size_t>(n)) S.r8.resize(n);
            uint8_t *__restrict r = S.r8.data();
            lookup_rank8(q, r, Sz - 1, n, static_cast<uint8_t>(smask),
                         smer_rank8.data());
            select_windows<uint8_t>(open, w, tt, r, sel, K - 1, n);
        } else if (Sz <= 8) {
            if (S.r16.size() < static_cast<size_t>(n)) S.r16.resize(n);
            uint16_t *__restrict r = S.r16.data();
            const uint16_t *rank = smer_rank16.data();
            for (int64_t p = Sz - 1; p < n; ++p) {
                r[p] = rank[((uint32_t(q[p - 4]) << 8) | q[p]) & smask];
            }
            select_windows<uint16_t>(open, w, tt, r, sel, K - 1, n);
        } else {
            if (S.h32.size() < static_cast<size_t>(n)) S.h32.resize(n);
            uint32_t *__restrict r = S.h32.data();
            for (int64_t p = Sz - 1; p < n; ++p)
                r[p] = hf(live_code(q, p, smask));
            select_windows<uint32_t>(open, w, tt, r, sel, K - 1, n);
        }
    }
};

}  // namespace syncmer

#endif
