#pragma once
#if defined(__SSSE3__)
#include <tmmintrin.h>
#endif

#ifndef LCMER_H
#define LCMER_H
#include <stdint.h>

#include <algorithm>
#include <bit>
#include <cstring>
#include <gtl/vector.hpp>
#include <iterator>
#include <type_traits>

#include "kmer_lcp.h"
#include "nibble2base.hpp"
#include "nt16.h"
#include "syncmer.h"
using SequenceContainerType = CrumbledSubSeq64;
// The fragments of one sequence, in order, wherever they are stored (a node's
// fragment list, a layout's vector).
using fragment_span = std::span<const SequenceContainerType>;

// Makes sure we have enough scratch space
template <class T, double GF = 2, size_t ALIGN_BLOCK_BITS = 5>
struct charsi_base {
    static constexpr size_t align_block_bits = ALIGN_BLOCK_BITS;
    static constexpr size_t align_block_mask = (1 << align_block_bits) - 1;
    T *data;
    size_t capacity;
    charsi_base() : data{NULL}, capacity{0} {}
    // Rule of 3 just delete
    charsi_base(const charsi_base &) = delete;
    charsi_base &operator=(const charsi_base &) = delete;
    // Rule of 5
    charsi_base(charsi_base &&other) noexcept
        : data{other.data}, capacity{other.capacity} {
        other.data = nullptr;
        other.capacity = 0;
    }
    charsi_base &operator=(charsi_base &&other) noexcept {
        if (this != &other) {
            free(data);
            data = other.data;
            capacity = other.capacity;
            other.data = nullptr;
            other.capacity = 0;
        }
        return *this;
    }

    ~charsi_base() { free(data); }
    int require(size_t cap) {
        if (cap > capacity) {
            capacity = cap * GF;
            free(data);
            size_t target_cap = capacity * sizeof(T);
            target_cap = (target_cap | align_block_mask) + 1;
            data = (decltype(data))aligned_alloc(32, target_cap);
        }
        return data == NULL;
    }
};
using charsi = charsi_base<uint8_t, 2.0, 5>;
using pos_scratch = charsi_base<uint32_t, 2.0, 5>;

struct litecmer {
    struct words_t {
        uint32_t hi;
        uint32_t lo;
    };
    words_t words;
    static constexpr uint32_t kPolyBase = 0x9e3779b1u;
    static_assert(sizeof(words_t) == sizeof(uint64_t));
    static_assert(std::is_trivially_copyable_v<words_t>);

    litecmer() : words{0, 0} {}
    litecmer(const litecmer &other) : words(other.words) {}
    litecmer(std::forward_iterator auto it, int n) : words{0, 0} {
        if (n > 1) {
            uint32_t h = 0;
            for (int i = 0; i < n - 1; ++i) {
                h = h * kPolyBase + *(it++);
            }
            words.hi = h;
            words.lo = *it;
        } else {
            words.lo = *it;
        }
    }

    litecmer &operator=(const litecmer &other) {
        words = other.words;
        return *this;
    }

    constexpr uint64_t value() const {
        return (static_cast<uint64_t>(words.hi) << 32) |
               static_cast<uint64_t>(words.lo);
    }

    static uint32_t poly_pow(uint32_t base, int exp) {
        uint32_t out = 1;
        for (int i = 0; i < exp; ++i) {
            out *= base;
        }
        return out;
    }

    // Order-sensitive rolling update for first N-1 cores.
    // Requires state initialized by the iterator constructor.
    void update_lcmer_poly(uint32_t new_core, uint32_t out_core, int n) {
        if (n > 1) {
            const uint32_t kOutMul = poly_pow(kPolyBase, n - 2);
            uint32_t h = words.hi;
            h -= out_core * kOutMul;
            h = h * kPolyBase + words.lo;
            words.hi = h;
        }
        words.lo = new_core;
    }

    void update_lcmer(uint32_t new_core, uint32_t out_core, int n) {
        update_lcmer_poly(new_core, out_core, n);
    }

    constexpr auto operator<=>(const litecmer &other) const {
        return value() <=> other.value();
    };
    constexpr bool operator==(const litecmer &other) const {
        return value() == other.value();
    }
};
using LCmer = litecmer;

template <>
struct std::hash<litecmer> {
    std::size_t operator()(const litecmer &l) const noexcept {
        return l.value();
    }
};

// Fused nt16 one-hot decode + RLE compaction in a single pass.
// Input: packed nt16 nibbles (2 bases per byte, as the BAM and FASTX readers
// produce), length in bases.
// Output: rle_data/rle_pos filled with RLE-compacted 2-bit bases. Runs are
// split on the one-hot nibble and decoded through kOnehotDecode (anything
// but A/C/G/T, e.g. N, decodes as A) identically on every path.
// The SIMD paths store whole vectors past the last written entry: callers
// must provide len + kOnehotRlePad writable entries in rle_data and rle_pos.
inline constexpr int64_t kOnehotRlePad = 32;
inline constexpr uint8_t kOnehotDecode[16] = {0, 0, 1, 0, 2, 0, 0, 0,
                                              3, 0, 0, 0, 0, 0, 0, 0};

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
// Binary reduction: compact 16 one-hot bytes into contiguous packed output.
// Each byte is 0 (unchanged) or (lane_idx << 4) | one_hot_value, one-hot
// {1,2,4,8} in the low nibble, lane 0-15 in the high nibble. Merge adjacent
// pairs in 3 rounds (16→8→4→2) via widen + variable-shift + OR, yielding two
// uint64 of up to 8 packed bytes each; unpack index (>>4) and value (&0xF),
// decode one-hot→2-bit, store via NEON.
static inline int compact16_emit(
    uint8x16_t data_vec,             // 16 one-hot bytes (nt16 nibble values)
    uint8x16_t neq_vec,              // 0xFF where changed, 0x00 where same
    const uint8x16_t &idx_template,  // precomputed (i<<4) for i=0..15
    const uint8x16_t &decode_tbl,    // one-hot→2-bit decode table (16 entries)
    uint8_t *rle_data,               // output: decoded 2-bit values
    uint32_t *rle_pos,               // output: base positions
    uint32_t pos_base  // position offset for this half (chunk_start or +16)
) {
    // Build combined byte per lane: (lane_idx << 4) | value for changed, 0 for
    // unchanged. A changed lane 0 of nibble 0 would pack to 0 as well; nibble
    // 15 decodes the same, so the value packs as that instead.
    const uint8x16_t zero_nibble = vceqq_u8(data_vec, vdupq_n_u8(0));
    const uint8x16_t value =
        vorrq_u8(data_vec, vandq_u8(zero_nibble, vdupq_n_u8(0x0F)));
    uint8x16_t live = vandq_u8(value, neq_vec);
    uint8x16_t combined = vorrq_u8(live, vandq_u8(neq_vec, idx_template));

    // --- Round 1: 16 bytes → 8 uint16 (merge adjacent pairs) ---
    // Deinterleave: even[k] = combined[2k], odd[k] = combined[2k+1]
    uint8x8x2_t uzp = vuzp_u8(vget_low_u8(combined), vget_high_u8(combined));
    uint8x8_t r1e = uzp.val[0];  // first of each pair
    uint8x8_t r1o = uzp.val[1];  // second of each pair
    uint8x8_t r1e_alive = vcgt_u8(r1e, vdup_n_u8(0));
    // Widen to uint16 so two bytes can coexist in one value
    uint16x8_t r1l = vmovl_u8(r1e);
    uint16x8_t r1r = vmovl_u8(r1o);
    // If left is alive, shift right element up by 8 bits (1 byte); else keep at
    // position 0
    int16x8_t r1s = vreinterpretq_s16_u16(
        vshlq_n_u16(vmovl_u8(vshr_n_u8(r1e_alive, 7)), 3));
    uint16x8_t r1m = vorrq_u16(r1l, vshlq_u16(r1r, r1s));
    // Count kept elements per pair (0, 1, or 2)
    uint8x8_t r1o_alive = vcgt_u8(r1o, vdup_n_u8(0));
    uint8x8_t r1cnt = vadd_u8(vshr_n_u8(r1e_alive, 7), vshr_n_u8(r1o_alive, 7));

    // --- Round 2: 8 uint16 → 4 uint32 (merge adjacent pairs of pairs) ---
    uint16x4x2_t uzp16 = vuzp_u16(vget_low_u16(r1m), vget_high_u16(r1m));
    uint32x4_t r2l = vmovl_u16(uzp16.val[0]);  // quads' left halves
    uint32x4_t r2r = vmovl_u16(uzp16.val[1]);  // quads' right halves
    // Left-pair counts: r1cnt at even positions
    uint8x8x2_t cnt_uzp = vuzp_u8(r1cnt, r1cnt);
    uint8x8_t r2lc = cnt_uzp.val[0];
    // Shift right half by (left_count * 8) bits
    uint16x4_t r2lc16 = vget_low_u16(vmovl_u8(r2lc));
    int32x4_t r2s = vreinterpretq_s32_u32(vmulq_n_u32(vmovl_u16(r2lc16), 8));
    uint32x4_t r2m = vorrq_u32(r2l, vshlq_u32(r2r, r2s));
    // Quad counts
    uint8x8_t r2rc = cnt_uzp.val[1];
    uint8x8_t r2cnt = vadd_u8(r2lc, r2rc);

    // --- Round 3: 4 uint32 → 2 uint64 (merge adjacent quads into octets) ---
    uint32x2x2_t uzp32 = vuzp_u32(vget_low_u32(r2m), vget_high_u32(r2m));
    uint64x2_t r3l = vmovl_u32(uzp32.val[0]);
    uint64x2_t r3r = vmovl_u32(uzp32.val[1]);
    int lc_a = vget_lane_u8(r2cnt, 0);  // left quad count of octet 0
    int lc_b = vget_lane_u8(r2cnt, 2);  // left quad count of octet 1
    int rc_a = vget_lane_u8(r2cnt, 1);  // right quad count of octet 0
    int rc_b = vget_lane_u8(r2cnt, 3);  // right quad count of octet 1
    int64x2_t r3s = {(int64_t)lc_a * 8, (int64_t)lc_b * 8};
    uint64x2_t r3m = vorrq_u64(r3l, vshlq_u64(r3r, r3s));

    int n_lo = lc_a + rc_a;  // packed byte count in r3m[0]
    int n_hi = lc_b + rc_b;  // packed byte count in r3m[1]
    int n_total = n_lo + n_hi;
    if (n_total == 0) return 0;

    // --- Unpack: extract values + indices from packed bytes, store via NEON
    // ---
    uint8x8_t blo = vreinterpret_u8_u64(vget_low_u64(r3m));
    uint8x8_t bhi = vreinterpret_u8_u64(vget_high_u64(r3m));

    // Values: low nibble is one-hot → decode to 2-bit via table lookup
    vst1_u8(rle_data, vqtbl1_u8(decode_tbl, vand_u8(blo, vdup_n_u8(0x0F))));
    if (n_hi > 0)
        vst1_u8(rle_data + n_lo,
                vqtbl1_u8(decode_tbl, vand_u8(bhi, vdup_n_u8(0x0F))));

    // Positions: high nibble is lane index → widen to uint32, add pos_base
    uint32x4_t cs = vdupq_n_u32(pos_base);
    uint16x8_t ilo16 = vmovl_u8(vshr_n_u8(blo, 4));
    vst1q_u32(rle_pos, vaddq_u32(vmovl_u16(vget_low_u16(ilo16)), cs));
    vst1q_u32(rle_pos + 4, vaddq_u32(vmovl_u16(vget_high_u16(ilo16)), cs));
    if (n_hi > 0) {
        uint16x8_t ihi16 = vmovl_u8(vshr_n_u8(bhi, 4));
        vst1q_u32(rle_pos + n_lo,
                  vaddq_u32(vmovl_u16(vget_low_u16(ihi16)), cs));
        vst1q_u32(rle_pos + n_lo + 4,
                  vaddq_u32(vmovl_u16(vget_high_u16(ihi16)), cs));
    }

    return n_total;
}
#elif defined(__SSSE3__)
// Left-pack table: for an 8-bit lane mask m, idx[m] lists the set lanes in
// ascending order, padded with 0x80 (pshufb writes 0 for those lanes).
struct onehot_left_pack_table {
    uint8_t idx[256][8];
    constexpr onehot_left_pack_table() : idx{} {
        for (int m = 0; m < 256; ++m) {
            int n = 0;
            for (int lane = 0; lane < 8; ++lane) {
                if (m & (1 << lane)) idx[m][n++] = static_cast<uint8_t>(lane);
            }
            while (n < 8) idx[m][n++] = 0x80;
        }
    }
};
inline constexpr onehot_left_pack_table kOnehotLeftPack{};

// Emits the run starts among 16 bases without a per-run loop. d: one-hot
// bytes; mask bit i set = lane i starts a run. Each 8-lane half is packed
// by one pshufb whose control bytes are also the lanes' indices, so the
// same control, widened, gives the positions. Stores are unconditional
// (up to 16 entries past the returned count).
static inline int onehot_emit16(__m128i d, unsigned mask,
                                const __m128i &decode_tbl, uint8_t *rle_data,
                                uint32_t *rle_pos, uint32_t pos_base) {
    const unsigned m_lo = mask & 0xFF, m_hi = mask >> 8;
    const int c_lo = std::popcount(m_lo);
    uint64_t lo_idx, hi_idx;
    std::memcpy(&lo_idx, kOnehotLeftPack.idx[m_lo], 8);
    std::memcpy(&hi_idx, kOnehotLeftPack.idx[m_hi], 8);
    // High half lanes are 8..15; 0x80 padding stays >= 0x80 after the add.
    const __m128i ctrl =
        _mm_set_epi64x(static_cast<long long>(hi_idx + 0x0808080808080808ULL),
                       static_cast<long long>(lo_idx));
    const __m128i codes =
        _mm_shuffle_epi8(decode_tbl, _mm_and_si128(d, _mm_set1_epi8(0x0F)));
    const __m128i packed = _mm_shuffle_epi8(codes, ctrl);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(rle_data), packed);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(rle_data + c_lo),
                     _mm_srli_si128(packed, 8));
#if defined(__AVX2__)
    const __m256i base = _mm256_set1_epi32(static_cast<int>(pos_base));
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(rle_pos),
                        _mm256_add_epi32(_mm256_cvtepu8_epi32(ctrl), base));
    _mm256_storeu_si256(
        reinterpret_cast<__m256i *>(rle_pos + c_lo),
        _mm256_add_epi32(_mm256_cvtepu8_epi32(_mm_srli_si128(ctrl, 8)), base));
#else
    const __m128i zero = _mm_setzero_si128();
    const __m128i base = _mm_set1_epi32(static_cast<int>(pos_base));
    const __m128i w_lo = _mm_unpacklo_epi8(ctrl, zero);  // 8 x u16
    const __m128i w_hi = _mm_unpackhi_epi8(ctrl, zero);
    auto store4 = [&](uint32_t *dst, __m128i v) {
        _mm_storeu_si128(reinterpret_cast<__m128i *>(dst),
                         _mm_add_epi32(v, base));
    };
    store4(rle_pos, _mm_unpacklo_epi16(w_lo, zero));
    store4(rle_pos + 4, _mm_unpackhi_epi16(w_lo, zero));
    store4(rle_pos + c_lo, _mm_unpacklo_epi16(w_hi, zero));
    store4(rle_pos + c_lo + 4, _mm_unpackhi_epi16(w_hi, zero));
#endif
    return c_lo + std::popcount(m_hi);
}
#endif

// Scalar RLE of bases [from, len): `prev` is the one-hot nibble of base
// from - 1. Shared by every path for the tail and as the portable fallback.
static inline int64_t onehot_to_rle_scalar(const uint8_t *data, int64_t len,
                                           int64_t from, uint8_t prev,
                                           uint8_t *rle_data, uint32_t *rle_pos,
                                           int64_t out) {
    for (int64_t i = from; i < len; ++i) {
        const uint8_t nib =
            (i & 1) ? (data[i >> 1] & 0x0F) : (data[i >> 1] >> 4);
        if (nib != prev) {
            rle_data[out] = kOnehotDecode[nib];
            rle_pos[out] = static_cast<uint32_t>(i);
            ++out;
        }
        prev = nib;
    }
    return out;
}

inline int64_t onehot_to_rle(const uint8_t *data, int64_t len,
                             uint8_t *rle_data, uint32_t *rle_pos) {
    if (len <= 0) return 0;

    // Base 0 always starts a run.
    rle_data[0] = kOnehotDecode[data[0] >> 4];
    rle_pos[0] = 0;
    int64_t out = 1;
    int64_t byte_idx = 0;

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    const uint8x16_t lo_mask = vdupq_n_u8(0x0F);
    const uint8x16_t decode_tbl = vld1q_u8(kOnehotDecode);

    // Lane index template: high nibble = lane index, used to tag each base with
    // its position
    static const uint8_t idx_arr[16] = {0x00, 0x10, 0x20, 0x30, 0x40, 0x50,
                                        0x60, 0x70, 0x80, 0x90, 0xA0, 0xB0,
                                        0xC0, 0xD0, 0xE0, 0xF0};
    const uint8x16_t idx_template = vld1q_u8(idx_arr);

    // Track last one-hot value for cross-chunk comparison.
    // One-hot format: A=1, C=2, G=4, T=8.
    uint8x16_t prev_last = vdupq_n_u8(data[0] >> 4);

    // Each iteration processes 16 input bytes → 32 bases.
    for (; byte_idx + 16 <= len / 2; byte_idx += 16) {
        // --- Phase 1: Unpack nt16 nibbles to one-hot bytes (no decode) ---
        uint8x16_t v = vld1q_u8(data + byte_idx);
        uint8x16_t hi = vshrq_n_u8(v, 4);      // high nibbles as bytes
        uint8x16_t lo = vandq_u8(v, lo_mask);  // low nibbles as bytes
        uint8x16x2_t z = vzipq_u8(hi, lo);
        uint8x16_t d0 = z.val[0];  // one-hot bases [0..15]
        uint8x16_t d1 = z.val[1];  // one-hot bases [16..31]

        // --- Phase 2: Build previous-base vectors for RLE comparison ---
        uint8x16_t p0 = vextq_u8(prev_last, d0, 15);
        uint8x16_t p1 = vextq_u8(d0, d1, 15);
        prev_last = d1;

        // --- Phase 3: Compare each base with its predecessor ---
        uint8x16_t neq0 = vmvnq_u8(vceqq_u8(d0, p0));
        uint8x16_t neq1 = vmvnq_u8(vceqq_u8(d1, p1));

        // Base 0 was already emitted; mark as unchanged.
        if (byte_idx == 0) neq0 = vsetq_lane_u8(0, neq0, 0);

        // Fast path: skip if all 32 bases match their predecessors.
        if (vmaxvq_u8(vorrq_u8(neq0, neq1)) == 0) continue;

        int64_t chunk_start = byte_idx * 2;

        // --- Phase 4: Binary reduction emit (no scalar ctz loop) ---
        if (vmaxvq_u8(neq0) != 0) {
            out += compact16_emit(d0, neq0, idx_template, decode_tbl,
                                  rle_data + out, rle_pos + out,
                                  static_cast<uint32_t>(chunk_start));
        }
        if (vmaxvq_u8(neq1) != 0) {
            out += compact16_emit(d1, neq1, idx_template, decode_tbl,
                                  rle_data + out, rle_pos + out,
                                  static_cast<uint32_t>(chunk_start + 16));
        }
    }
#elif defined(__SSSE3__)
    const __m128i lo_mask = _mm_set1_epi8(0x0F);
    const __m128i decode_tbl =
        _mm_loadu_si128(reinterpret_cast<const __m128i *>(kOnehotDecode));
    // Last one-hot base of the previous block, in lane 15.
    __m128i prev_last = _mm_set1_epi8(static_cast<char>(data[0] >> 4));

    // Each iteration processes 16 input bytes → 32 bases.
    for (; byte_idx + 16 <= len / 2; byte_idx += 16) {
        // Unpack nt16 nibbles to one-hot bytes, bases [0..15] and [16..31].
        const __m128i v =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(data + byte_idx));
        const __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), lo_mask);
        const __m128i lo = _mm_and_si128(v, lo_mask);
        const __m128i d0 = _mm_unpacklo_epi8(hi, lo);
        const __m128i d1 = _mm_unpackhi_epi8(hi, lo);
        // Predecessor of each base; a lane starts a run if it differs.
        const __m128i p0 = _mm_alignr_epi8(d0, prev_last, 15);
        const __m128i p1 = _mm_alignr_epi8(d1, d0, 15);
        prev_last = d1;
        unsigned m0 = ~_mm_movemask_epi8(_mm_cmpeq_epi8(d0, p0)) & 0xFFFFu;
        const unsigned m1 =
            ~_mm_movemask_epi8(_mm_cmpeq_epi8(d1, p1)) & 0xFFFFu;
        if (byte_idx == 0) m0 &= ~1u;  // base 0 was already emitted

        const uint32_t chunk_start = static_cast<uint32_t>(byte_idx * 2);
        out += onehot_emit16(d0, m0, decode_tbl, rle_data + out, rle_pos + out,
                             chunk_start);
        out += onehot_emit16(d1, m1, decode_tbl, rle_data + out, rle_pos + out,
                             chunk_start + 16);
    }
#endif
    // Remaining bases (all of them without SIMD). The predecessor of base
    // 2 * byte_idx is the low nibble of the previous byte.
    const int64_t from = std::max<int64_t>(1, 2 * byte_idx);
    const uint8_t prev =
        (from & 1) ? (data[0] >> 4) : (data[byte_idx - 1] & 0x0F);
    return onehot_to_rle_scalar(data, len, from, prev, rle_data, rle_pos, out);
}

// Packs 64 selection bytes (each 0 or 1) into a 64-bit mask, bit i = s[i].
inline uint64_t sel_mask64(const uint8_t *s) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    // Shift byte i left by i % 8, then three pairwise-add levels sum each run
    // of 8 bytes into one byte: byte j of the result holds bits 8j..8j+7.
    const uint8x16x4_t v = vld1q_u8_x4(s);
    static constexpr int8_t kShift[16] = {0, 1, 2, 3, 4, 5, 6, 7,
                                          0, 1, 2, 3, 4, 5, 6, 7};
    const int8x16_t sh = vld1q_s8(kShift);
    const uint8x16_t p =
        vpaddq_u8(vshlq_u8(v.val[0], sh), vshlq_u8(v.val[1], sh));
    const uint8x16_t q =
        vpaddq_u8(vshlq_u8(v.val[2], sh), vshlq_u8(v.val[3], sh));
    const uint8x16_t r = vpaddq_u8(vpaddq_u8(p, q), vdupq_n_u8(0));
    return vgetq_lane_u64(vreinterpretq_u64_u8(r), 0);
#elif defined(__AVX2__)
    // Shifting each 16-bit lane left by 7 moves bit 0 of every byte to its
    // bit 7, which movemask collects.
    auto half = [](const uint8_t *p) {
        const __m256i v =
            _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        return static_cast<uint32_t>(
            _mm256_movemask_epi8(_mm256_slli_epi16(v, 7)));
    };
    return half(s) | (uint64_t{half(s + 32)} << 32);
#elif defined(__SSE2__)
    uint64_t mask = 0;
    for (int b = 0; b < 4; ++b) {
        const __m128i v =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 16 * b));
        mask |= uint64_t{static_cast<uint32_t>(
                    _mm_movemask_epi8(_mm_slli_epi16(v, 7)))}
                << (16 * b);
    }
    return mask;
#else
    // 8 bytes of 0/1 → 8 bits via one multiply per word.
    uint64_t mask = 0;
    for (int b = 0; b < 8; ++b) {
        uint64_t w;
        std::memcpy(&w, s + 8 * b, sizeof w);
        mask |= ((w * 0x0102040810204080ULL) >> 56) << (8 * b);
    }
    return mask;
#endif
}

// Level-1 cores of an RLE stream into str (resets it): the selected k-mers.
inline uint64_t parse1_rle_syncmers(const uint8_t *rle_data,
                                    const uint32_t *rle_pos, int64_t rle_len,
                                    int64_t seq_len, core_buf *str,
                                    const auto &syncmer_functor) {
    const auto k = syncmer_functor.k();
    str->level = 1;
    str->kmer_bits = 2 * static_cast<int>(k);
    str->size = 0;
    if (rle_len < k) return 0;

    // Selection for every window at once, then emit only the selected ones:
    // iterate set bits of a 64-window mask (sel_mask64); each k-mer code is 4
    // loads of 4-symbol codes. A core spans its RLE-space k-mer (runs
    // [e-k+1, e]) in original base coordinates: start of the first run up to
    // the start of the next run.
    thread_local syncmer::live_scratch scratch;
    constexpr int64_t kBlock = 64;
    scratch.require(static_cast<size_t>(rle_len + kBlock));
    uint8_t *sel = scratch.sel.data();
    std::memset(sel, 0, k - 1);
    std::memset(sel + rle_len, 0, kBlock);
    syncmer_functor.select_live(rle_data, rle_len, scratch);
    const uint8_t *q = scratch.codes();
    const uint32_t kmask = static_cast<uint32_t>((1ULL << (2 * k)) - 1);
    // Room for exactly the selected k-mers (a fraction of the windows).
    int64_t selected = 0;
    for (int64_t base = 0; base < rle_len; base += kBlock)
        selected += std::popcount(sel_mask64(sel + base));
    str->reserve(selected);
    int64_t core_index = 0;
    uint32_t *label = str->label;
    uint64_t *span = str->span;
    for (int64_t base = 0; base < rle_len; base += kBlock) {
        uint64_t mask = sel_mask64(sel + base);
        while (mask) {
            const int64_t e = base + std::countr_zero(mask);
            mask &= mask - 1;
            const uint32_t si = rle_pos[e - k + 1];
            const uint32_t ei = (e + 1 < rle_len)
                                    ? rle_pos[e + 1]
                                    : static_cast<uint32_t>(seq_len);
            label[core_index] = syncmer::live_code(q, e, kmask);
            span[core_index] = si | static_cast<uint64_t>(ei) << 32;
            core_index++;
        }
    }
    str->size = core_index;
    return core_index;
}

// Level-1 cores of bases [from, to) of the concatenated fragments, with
// positions in whole-sequence coordinates: the RLE stream of the range is
// windowed with the same code as the read path. A run is NOT restarted when
// the next fragment begins with the run's base; the range's first run starts
// at `from`, and a k-mer ending in its last run ends at `to`.
inline uint64_t parse1_CrumbledSubSeqs_syncmers(fragment_span sseqs,
                                                size_t from, size_t to,
                                                core_buf *str,
                                                const auto &syncmer_functor) {
    if (sseqs.empty()) {
        fprintf(stderr, "Error: Empty sequence list\n");
        str->level = 1;
        str->size = 0;
        return 0;
    }
    const auto [rle_data, rle_pos] = str->rle_scratch(to - from);
    int64_t rle_len = 0;
    size_t l_so_far = 0;
    uint8_t prev_base = 255;
    for (const auto &s : sseqs) {
        if (l_so_far >= to) break;
        // [lo, hi): the part of this fragment inside the range.
        const size_t lo = std::max(from, l_so_far) - l_so_far;
        const size_t hi = std::min(to, l_so_far + s.len) - l_so_far;
        for (size_t i = lo; i < hi; ++i) {
            const uint8_t b = s.code(i);
            if (b != prev_base) {
                rle_data[rle_len] = b;
                rle_pos[rle_len] = static_cast<uint32_t>(l_so_far + i);
                ++rle_len;
                prev_base = b;
            }
        }
        l_so_far += s.len;
    }
    return parse1_rle_syncmers(rle_data, rle_pos, rle_len,
                               static_cast<int64_t>(to), str, syncmer_functor);
}

// Level-1 cores of the whole concatenated fragments.
inline void init_lps_crumbled(core_buf *str, fragment_span sseqs,
                              const auto &syncmer_functor) {
    size_t len = 0;
    for (const auto &s : sseqs) len += s.len;
    parse1_CrumbledSubSeqs_syncmers(sseqs, 0, len, str, syncmer_functor);
}

// Bases per window of lps_parse_crumbled. Windows bound the per-thread parse
// buffers (tens of bytes per base) independently of the sequence length,
// except where an exact periodic stretch longer than the window forces the
// window to double until it covers the stretch; the whole-sequence parse
// holds that stretch anyway. The window is large against the re-parsed
// overlap (the context of a few final-level cores), so re-parsing stays a
// small fraction of the work.
inline constexpr uint64_t kLpsWindowBases = uint64_t{1} << 17;

// Window start of a core that no later window can be anchored on.
inline constexpr uint64_t kLpsNoReq = UINT64_MAX;

// A window's cores, the caller's scratch, and the per-core req of
// lps_deepen_window with that function's scratch, reused across the windows
// of one lps_parse_crumbled.
struct lps_window {
    core_buf *cores;
    gtl::vector<uint64_t> req, next_req;
};

// One deepening level on a window's cores. `w.cores` holds a run of
// consecutive whole-parse cores; `left_open`: the run does not start the
// sequence. core_buf::deepen gives lps_deepen1's cores for the same input
// (DCT, then parse3's rules), so the argument below is parse3's.
// - Right end: DCT compresses a core against its left neighbour, and parse3
//   decides at position j from positions j-1..j+3 and the equal-core run
//   after j. A decision that would need cores past the end is not taken
//   (no LMAX in the last 3 positions, no RINT on a run that reaches the
//   end) and nothing is taken after it, so the outputs of a prefix are a
//   prefix of the outputs.
// - Left end: position 0 is not compressed and is never a constituent, and
//   only an LMAX at position 1 (parse3's begin) would read it; an SSEQ core
//   is emitted only between two special cores. So the window can miss the
//   whole parse's outputs before its first special core but emits no core
//   the whole parse lacks: its outputs are again a run of consecutive
//   whole-parse cores, and none is trimmed.
// - req[i] is a window start that guarantees core i: a window starting
//   there, or earlier, emits it if it reaches far enough right. Output i-1,
//   if special, is decided at its first constituent first(i-1); if an SSEQ
//   core, output i is special and decided after it. A window holding input
//   core first(i-1) - 2 decides both at position 2 or later, where every
//   decision is the whole parse's, so req[i] is that core's req.
// Returns false when the parse ends (fewer than two cores).
inline bool lps_deepen_window(lps_window &w, bool left_open) {
    core_buf *s = w.cores;
    const int64_t n = s->size;
    s->deepen(s->level + 1);
    // Below two cores, a level leaves none.
    if (n < 2) return false;
    const int64_t m = s->size;
    w.next_req.resize(m);
    for (int64_t i = 0; i < m; ++i) {
        uint64_t r = 0;  // a window starting at the sequence start emits it
        if (i > 0 && s->first(i - 1) >= 2) {
            r = w.req[s->first(i - 1) - 2];
        } else if (left_open) {
            r = kLpsNoReq;
        }
        w.next_req[i] = r;
    }
    std::swap(w.req, w.next_req);
    return true;
}

// Cores of bases [from, to) of `total`, deepened to `lcp_level` and trimmed
// to those equal to the whole-sequence parse; req as in lps_deepen_window.
// A level-1 core (an RLE k-mer) is exact unless its first run may be cut by
// `from` or it ends at `to`; a window starting before it keeps it.
inline void lps_parse_window(lps_window &w, fragment_span sseqs, uint64_t from,
                             uint64_t to, uint64_t total,
                             const auto &syncmer_functor, int lcp_level) {
    const bool left_open = from > 0;
    const bool right_open = to < total;
    core_buf *win = w.cores;
    parse1_CrumbledSubSeqs_syncmers(sseqs, from, to, win, syncmer_functor);
    if (right_open && win->size > 0 && win->end(win->size - 1) == to)
        --win->size;
    if (left_open && win->size > 0 && win->start(0) == from) win->drop_first();
    w.req.resize(win->size);
    for (int64_t i = 0; i < win->size; ++i)
        w.req[i] = std::max<uint64_t>(win->start(i), 1) - 1;
    while (win->level < lcp_level && lps_deepen_window(w, left_open)) {
    }
}

// Index of the first core after the anchor core (`label`, `pos`: start, end)
// among a window's cores, or -1 when the window does not hold the anchor.
inline int64_t lps_after_anchor(const core_buf *win, uint32_t label,
                                std::pair<uint32_t, uint32_t> pos) {
    int64_t i = 0;
    while (i < win->size && win->start(i) < pos.first) ++i;
    if (i < win->size && win->end(i) == pos.second && win->label[i] == label)
        return i + 1;
    return -1;
}

// Appends to `out` the level-`lcp_level` cores of the concatenated fragments,
// equal to init_lps_crumbled + deepen on the whole sequence, parsed in
// windows of `window` bases so `win`, the caller's scratch, does not grow with
// the sequence. Each later window starts at the req of the last core
// appended to `out`, so it parses that anchor core again and the two runs of
// cores join there. A window that cannot reach past its anchor doubles until
// it can.
inline void lps_parse_crumbled(lcp_read &out, core_buf *win,
                               fragment_span sseqs, const auto &syncmer_functor,
                               int lcp_level,
                               uint64_t window = kLpsWindowBases) {
    uint64_t total = 0;
    for (const auto &s : sseqs) total += s.len;
    lps_window w{win, {}, {}};
    const size_t out0 = out.cores.size();
    uint64_t from = 0;
    uint64_t extent = window;
    for (;;) {
        const uint64_t to = std::min(total, from + extent);
        lps_parse_window(w, sseqs, from, to, total, syncmer_functor, lcp_level);
        int64_t i = 0;
        if (out.cores.size() > out0)
            i = lps_after_anchor(win, out.cores.back(), out.pos.back());
        if (i < 0 && to < total) {
            extent *= 2;
            continue;
        }
        if (i < 0) {
            fprintf(stderr, "lps_parse_crumbled: windows do not join\n");
            abort();
        }
        if (to == total) {
            win->append_to(out, i, win->size);
            break;
        }
        int64_t a = win->size - 1;
        while (a >= i && w.req[a] == kLpsNoReq) --a;
        if (a < i) {
            extent *= 2;
            continue;
        }
        win->append_to(out, i, a + 1);
        // From the same start the window would repeat, holding nothing past a.
        if (w.req[a] == from) {
            extent *= 2;
        } else {
            extent = window;
        }
        from = w.req[a];
    }
}

// Complement of a BAM nibble (A<->T, C<->G, anything else N).
constexpr uint8_t nibble_complement(uint8_t n) {
    constexpr uint8_t comp[16] = {15, 8,  4,  15, 2,  15, 15, 15,
                                  1,  15, 15, 15, 15, 15, 15, 15};
    return comp[n];
}

// Nonzero iff nibble n is not A, C, G or T, the nibbles with exactly one bit
// set. Plain bit operations, so a scan of a whole read vectorizes.
constexpr uint8_t non_acgt_bits(uint8_t n) {
    return static_cast<uint8_t>((n & (n - 1)) | (n == 0));
}

// Whether any of len packed nibbles is not A, C, G or T.
inline bool has_non_acgt(const uint8_t *nibbles, int64_t len) {
    uint8_t other = 0;
    for (int64_t i = 0; i < len / 2; ++i) {
        other |=
            non_acgt_bits(nibbles[i] >> 4) | non_acgt_bits(nibbles[i] & 0x0F);
    }
    if (len % 2 == 1) other |= non_acgt_bits(nibbles[len / 2] >> 4);
    return other != 0;
}

// Turns the run-length stream onehot_to_rle wrote for `nibbles` (len bases,
// rle_len runs) in place into the one it writes for their reverse complement,
// and returns its run count. The runs come in reverse order, and the run that
// ended at e starts at len - e. An A, C, G or T run's base b complements to
// 3 - b. Every other nibble complements to N, which decodes as A, so such a
// run's base is 0 and neighbouring runs of two such nibbles become one run.
inline int64_t rle_reverse_complement(const uint8_t *nibbles, int64_t len,
                                      uint8_t *rle_data, uint32_t *rle_pos,
                                      int64_t rle_len) {
    if (rle_len <= 0) return 0;
    std::reverse(rle_data, rle_data + rle_len);
    for (int64_t j = 0; j < rle_len; ++j) rle_data[j] ^= 3;
    // Forward run i ends where run i + 1 starts: once rle_pos[1, rle_len) is
    // reversed, reverse run j >= 1 starts at len - rle_pos[j].
    std::reverse(rle_pos + 1, rle_pos + rle_len);
    for (int64_t j = 1; j < rle_len; ++j) {
        rle_pos[j] = static_cast<uint32_t>(len) - rle_pos[j];
    }
    if (!has_non_acgt(nibbles, len)) return rle_len;

    // Some run is not A, C, G or T: its 3 - b is wrong, and it merges with a
    // neighbouring such run.
    int64_t out = 0;
    bool prev_other = false;
    for (int64_t j = 0; j < rle_len; ++j) {
        // Reverse base rle_pos[j] complements forward base `at`, which is the
        // high nibble of its byte when even.
        const uint32_t at = static_cast<uint32_t>(len) - 1 - rle_pos[j];
        const uint8_t nib = (nibbles[at >> 1] >> ((~at & 1u) << 2)) & 0x0F;
        const bool other = non_acgt_bits(nib) != 0;
        if (other && prev_other) continue;
        rle_data[out] = rle_data[j];
        if (other) rle_data[out] = kOnehotDecode[nibble_complement(nib)];
        rle_pos[out] = rle_pos[j];
        ++out;
        prev_other = other;
    }
    return out;
}

// Cores of the run-length stream rle_data/rle_pos (rle_len runs of a len-base
// sequence) at lcp_level, appended to lr.
inline void parse_cores_from_rle(const uint8_t *rle_data,
                                 const uint32_t *rle_pos, int64_t rle_len,
                                 size_t len, int lcp_level, core_buf *str,
                                 const auto &syncmer_functor, lcp_read &lr) {
    parse1_rle_syncmers(rle_data, rle_pos, rle_len, static_cast<int64_t>(len),
                        str, syncmer_functor);
    str->deepen(lcp_level);
    str->append_to(lr);
}

// Cores of len packed nibbles, appended to lr. The run-length stream stays in
// chs/pos_buf; returns its run count.
inline int64_t parse_cores_from_nibbles(const uint8_t *data, size_t len,
                                        int lcp_level, core_buf *str,
                                        const auto &syncmer_functor,
                                        charsi &chs, lcp_read &lr,
                                        pos_scratch &pos_buf) {
    chs.require(len + kOnehotRlePad);
    pos_buf.require(len + kOnehotRlePad);
    const int64_t rle_len = onehot_to_rle(data, len, chs.data, pos_buf.data);
    parse_cores_from_rle(chs.data, pos_buf.data, rle_len, len, lcp_level, str,
                         syncmer_functor, lr);
    lr.data = data;
    return rle_len;
}
#endif
