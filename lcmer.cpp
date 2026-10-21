#include "lcmer.h"

#include <bit>
#include <cstdint>

#include "asm.h"
#include "ext/lcptools/lcptools_ho.h"
#include "fastx.h"
#include "gtl/phmap.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "nibble2base.hpp"
#include "syncmer.h"

uint8_t *onehotto2bitDNA(uint8_t *data, size_t len) {
    uint8_t *ret = (uint8_t *)malloc(len);
    for (size_t i = 0; i < len / 2; ++i) {
        ret[2 * i] = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        ret[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
    return ret;
}
void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &chs) {
    chs.require(len);
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    static const uint8_t decode4_tbl_arr[16] = {
        0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0};
    const uint8x16_t tbl = vld1q_u8(decode4_tbl_arr);
    const uint8x16_t mask = vdupq_n_u8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len / 2; i += 16) {
        uint8x16_t v = vld1q_u8(data + i);
        uint8x16_t hi = vqtbl1q_u8(tbl, vshrq_n_u8(v, 4));
        uint8x16_t lo = vqtbl1q_u8(tbl, vandq_u8(v, mask));
        uint8x16x2_t interleaved = vzipq_u8(hi, lo);
        vst1q_u8(chs.data + 2 * i,      interleaved.val[0]);
        vst1q_u8(chs.data + 2 * i + 16, interleaved.val[1]);
    }
    for (; i < len / 2; ++i) {
        chs.data[2 * i]     = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        chs.data[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
#elif defined(__SSSE3__)
    static const uint8_t decode4_tbl_arr[16] = {
        0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0};
    const __m128i tbl = _mm_loadu_si128(
        reinterpret_cast<const __m128i *>(decode4_tbl_arr));
    const __m128i mask = _mm_set1_epi8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len / 2; i += 16) {
        __m128i v = _mm_loadu_si128(
            reinterpret_cast<const __m128i *>(data + i));
        __m128i hi = _mm_shuffle_epi8(tbl,
            _mm_and_si128(_mm_srli_epi16(v, 4), mask));
        __m128i lo = _mm_shuffle_epi8(tbl, _mm_and_si128(v, mask));
        __m128i interleaved_lo = _mm_unpacklo_epi8(hi, lo);
        __m128i interleaved_hi = _mm_unpackhi_epi8(hi, lo);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + 2 * i), interleaved_lo);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + 2 * i + 16), interleaved_hi);
    }
    for (; i < len / 2; ++i) {
        chs.data[2 * i]     = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        chs.data[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
#else
    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * i]     = 27 - std::countl_zero(static_cast<unsigned>(data[i]));
        chs.data[2 * i + 1] = std::countr_zero(static_cast<unsigned>(data[i]));
    }
#endif
    if (len % 2 == 1) {
        chs.data[len - 1] = 27 - std::countl_zero(static_cast<unsigned>(data[len / 2]));
    }
}
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &chs) {
    chs.require(len);
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    static const uint8_t rc_tbl_arr[16] = {
        0,3,2,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};
    const uint8x16_t tbl = vld1q_u8(rc_tbl_arr);
    const uint8x16_t mask = vdupq_n_u8(0x0F);
    size_t half = len / 2;
    size_t i = 0;
    for (; i + 16 <= half; i += 16) {
        uint8x16_t v = vld1q_u8(data + i);
        uint8x16_t lo = vqtbl1q_u8(tbl, vandq_u8(v, mask));
        uint8x16_t hi = vqtbl1q_u8(tbl, vshrq_n_u8(v, 4));
        uint8x16x2_t interleaved = vzipq_u8(lo, hi);
        uint8x16_t r0 = vrev64q_u8(interleaved.val[1]);
        r0 = vcombine_u8(vget_high_u8(r0), vget_low_u8(r0));
        uint8x16_t r1 = vrev64q_u8(interleaved.val[0]);
        r1 = vcombine_u8(vget_high_u8(r1), vget_low_u8(r1));
        size_t out_pos = len - 2 * (i + 16);
        vst1q_u8(chs.data + out_pos,      r0);
        vst1q_u8(chs.data + out_pos + 16, r1);
    }
    size_t mid = half - 1;
    for (; i < half; ++i) {
        chs.data[2 * (mid - i)]     = std::countl_zero(static_cast<unsigned>(data[i])) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - std::countr_zero(static_cast<unsigned>(data[i]));
    }
#elif defined(__SSSE3__)
    static const uint8_t rc_tbl_arr[16] = {
        0,3,2,0, 1,0,0,0, 0,0,0,0, 0,0,0,0};
    const __m128i tbl = _mm_loadu_si128(
        reinterpret_cast<const __m128i *>(rc_tbl_arr));
    const __m128i mask = _mm_set1_epi8(0x0F);
    // Byte-reverse pattern for _mm_shuffle_epi8: reverses 16 bytes
    const __m128i rev_pattern = _mm_set_epi8(
        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    size_t half = len / 2;
    size_t i = 0;
    for (; i + 16 <= half; i += 16) {
        __m128i v = _mm_loadu_si128(
            reinterpret_cast<const __m128i *>(data + i));
        __m128i lo = _mm_shuffle_epi8(tbl, _mm_and_si128(v, mask));
        __m128i hi = _mm_shuffle_epi8(tbl,
            _mm_and_si128(_mm_srli_epi16(v, 4), mask));
        // Interleave lo nibble first, then hi nibble (reversed order for the
        // complement).
        __m128i interleaved_lo = _mm_unpacklo_epi8(lo, hi);
        __m128i interleaved_hi = _mm_unpackhi_epi8(lo, hi);
        // Reverse both 16-byte vectors
        __m128i r0 = _mm_shuffle_epi8(interleaved_hi, rev_pattern);
        __m128i r1 = _mm_shuffle_epi8(interleaved_lo, rev_pattern);
        size_t out_pos = len - 2 * (i + 16);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + out_pos), r0);
        _mm_storeu_si128(
            reinterpret_cast<__m128i *>(chs.data + out_pos + 16), r1);
    }
    size_t mid = half - 1;
    for (; i < half; ++i) {
        chs.data[2 * (mid - i)]     = std::countl_zero(static_cast<unsigned>(data[i])) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - std::countr_zero(static_cast<unsigned>(data[i]));
    }
#else
    size_t mid = len / 2 - 1;
    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * (mid - i)]     = std::countl_zero(static_cast<unsigned>(data[i])) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - std::countr_zero(static_cast<unsigned>(data[i]));
    }
#endif
    if (len % 2 == 1) {
        chs.data[0] = std::countl_zero(static_cast<unsigned>(data[len / 2])) - 24;
    }
}

#ifdef LASM_UNIT_TESTS
#include "tests/lcmer.test.inc"
#endif
