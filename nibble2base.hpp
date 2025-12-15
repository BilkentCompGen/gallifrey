#ifndef NIBBLE2BASE_HPP
#define NIBBLE2BASE_HPP
// #include <sys/_types/_u_char.h>
#include <algorithm>
#include "htslib/hts.h"

typedef unsigned int u_int;
typedef unsigned char u_char;
typedef unsigned short u_short;
#include <htslib/sam.h>
#include <sys/sysctl.h>
#include <sys/types.h>

#include <cstdint>
enum class CPUType { X86, ARM, UNKNOWN };
#define BUILDING_SIMD_NIBBLE2BASE
static inline void nibble2base_default(uint8_t *nib, char *seq, int len) {
    static const char code2base[513] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len / 2;
    seq[0] = 0;

    for (i = 0; i < len2; i++)
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i * 2], &code2base[(size_t)nib[i] * 2], 2);

    if ((i *= 2) < len) seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

#if defined HAVE_ATTRIBUTE_CONSTRUCTOR &&                                   \
    ((defined __x86_64__ && defined HAVE_X86INTRIN_H && HAVE_X86INTRIN_H && \
      defined HAVE_ATTRIBUTE_TARGET_SSSE3 &&                                \
      defined HAVE_BUILTIN_CPU_SUPPORT_SSSE3) ||                            \
     (defined __ARM_NEON))
#define BUILDING_SIMD_NIBBLE2BASE
#endif

static inline void nibble2base(uint8_t *nib, char *seq, int len) {
#ifdef BUILDING_SIMD_NIBBLE2BASE
    extern void (*htslibx_nibble2base)(uint8_t *nib, char *seq, int len);
    htslibx_nibble2base(nib, seq, len);
#else
    nibble2base_default(nib, seq, len);
#endif
}

// These must be defined before the first system include to ensure that legacy
// BSD types needed by <sys/sysctl.h> remain defined when _XOPEN_SOURCE is set.
#if defined __APPLE__
#define _DARWIN_C_SOURCE
#elif defined __NetBSD__
#define _NETBSD_SOURCE
#endif

#include "htslib/sam.h"
// #include "sam_internal.h"

#if defined __x86_64__ && defined HAVE_X86INTRIN_H && HAVE_X86INTRIN_H
#include <x86intrin.h>
#elif defined __ARM_NEON
#include <arm_neon.h>
#endif

#if defined __arm__ || defined __aarch64__

#ifdef HAVE_OPENBSD
/*
 * Extra check for elf_aux_info() on configure-less OpenBSD builds.  Once
 * version 7.5 has dropped off support, this can be changed to an assumption
 * that the function exists in the Makefile-generated config.h.
 */
#include <sys/param.h>
#if OpenBSD >= 202409
#define HAVE_ELF_AUX_INFO
#endif
#endif

#if defined HAVE_GETAUXVAL || defined HAVE_ELF_AUX_INFO
#include <sys/auxv.h>
#elif defined __APPLE__
#include <sys/sysctl.h>
#include <sys/types.h>
#elif defined __NetBSD__
#include <stddef.h>
#include <sys/param.h>
#include <sys/sysctl.h>
#ifdef __aarch64__
#include <aarch64/armreg.h>
#else
#include <arm/armreg.h>
#endif
#elif defined _WIN32
#include <processthreadsapi.h>
#endif

#endif
#ifdef BUILDING_SIMD_NIBBLE2BASE

inline void (*htslibx_nibble2base)(uint8_t *nib, char *seq,
                                   int len) = nibble2base_default;

#if defined __x86_64__
#include <emmintrin.h>
#include <tmmintrin.h>
/*
 * Convert a nibble encoded BAM sequence to a string of bases.
 *
 * Using SSSE3 instructions, 16 codepoints that hold 2 bases each can be
 * unpacked into 32 indexes from 0-15. Using the pshufb instruction these can
 * be converted to the IUPAC characters.
 * It falls back on the nibble2base_default function for the remainder.
 */

__attribute__((target("ssse3"))) static void nibble2base_ssse3(uint8_t *nib,
                                                               char *seq,
                                                               int len) {
    const char *seq_end_ptr = seq + len;
    char *seq_cursor = seq;
    uint8_t *nibble_cursor = nib;
    const char *seq_vec_end_ptr = seq_end_ptr - (2 * sizeof(__m128i) - 1);
    __m128i nuc_lookup_vec = _mm_lddqu_si128((__m128i *)seq_nt16_str);
    /* Nucleotides are encoded 4-bits per nucleotide and stored in 8-bit bytes
       as follows: |AB|CD|EF|GH|. The 4-bit codes (going from 0-15) can be used
       together with the pshufb instruction as a lookup table. The most
       efficient way is to use bitwise AND and shift to create two vectors. One
       with all the upper codes (|A|C|E|G|) and one with the lower codes
       (|B|D|F|H|). The lookup can then be performed and the resulting vectors
       can be interleaved again using the unpack instructions. */
    while (seq_cursor < seq_vec_end_ptr) {
        __m128i encoded = _mm_lddqu_si128((__m128i *)nibble_cursor);
        __m128i encoded_upper = _mm_srli_epi64(encoded, 4);
        encoded_upper = _mm_and_si128(encoded_upper, _mm_set1_epi8(15));
        __m128i encoded_lower = _mm_and_si128(encoded, _mm_set1_epi8(15));
        __m128i nucs_upper = _mm_shuffle_epi8(nuc_lookup_vec, encoded_upper);
        __m128i nucs_lower = _mm_shuffle_epi8(nuc_lookup_vec, encoded_lower);
        __m128i first_nucleotides = _mm_unpacklo_epi8(nucs_upper, nucs_lower);
        __m128i second_nucleotides = _mm_unpackhi_epi8(nucs_upper, nucs_lower);
        _mm_storeu_si128((__m128i *)seq_cursor, first_nucleotides);
        _mm_storeu_si128((__m128i *)(seq_cursor + sizeof(__m128i)),
                         second_nucleotides);
        nibble_cursor += sizeof(__m128i);
        seq_cursor += 2 * sizeof(__m128i);
    }
    nibble2base_default(nibble_cursor, seq_cursor, seq_end_ptr - seq_cursor);
}

__attribute__((constructor)) static CPUType nibble2base_resolve(void) {
    if (__builtin_cpu_supports("ssse3")) {
        htslibx_nibble2base = nibble2base_ssse3;
        return CPUType::X86;
    }
    return CPUType::UNKNOWN;
}

#elif defined __ARM_NEON

static void nibble2base_neon(uint8_t *nib, char *seq0, int len) {
    uint8x16_t low_nibbles_mask = vdupq_n_u8(0x0f);
    uint8x16_t nuc_lookup_vec = vld1q_u8((const uint8_t *)seq_nt16_str);
#ifndef __aarch64__
    uint8x8x2_t nuc_lookup_vec2 = {
        {vget_low_u8(nuc_lookup_vec), vget_high_u8(nuc_lookup_vec)}};
#endif

    uint8_t *seq = (uint8_t *)seq0;
    int blocks;

    for (blocks = len / 32; blocks > 0; --blocks) {
        uint8x16_t encoded = vld1q_u8(nib);
        nib += 16;

        /* Translate the high and low nibbles to nucleotide letters separately,
           then interleave them back together via vzipq for writing. */

        uint8x16_t high_nibbles = vshrq_n_u8(encoded, 4);
        uint8x16_t low_nibbles = vandq_u8(encoded, low_nibbles_mask);

#ifdef __aarch64__
        uint8x16_t high_nucleotides = vqtbl1q_u8(nuc_lookup_vec, high_nibbles);
        uint8x16_t low_nucleotides = vqtbl1q_u8(nuc_lookup_vec, low_nibbles);
#else
        uint8x8_t high_low =
            vtbl2_u8(nuc_lookup_vec2, vget_low_u8(high_nibbles));
        uint8x8_t high_high =
            vtbl2_u8(nuc_lookup_vec2, vget_high_u8(high_nibbles));
        uint8x16_t high_nucleotides = vcombine_u8(high_low, high_high);

        uint8x8_t low_low = vtbl2_u8(nuc_lookup_vec2, vget_low_u8(low_nibbles));
        uint8x8_t low_high =
            vtbl2_u8(nuc_lookup_vec2, vget_high_u8(low_nibbles));
        uint8x16_t low_nucleotides = vcombine_u8(low_low, low_high);
#endif

#ifdef __aarch64__
        vst1q_u8_x2(seq, vzipq_u8(high_nucleotides, low_nucleotides));
#else
        // Avoid vst1q_u8_x2 as GCC erroneously omits it on 32-bit ARM
        uint8x16x2_t nucleotides = {{high_nucleotides, low_nucleotides}};
        vst2q_u8(seq, nucleotides);
#endif
        seq += 32;
    }

    if (len % 32 != 0) nibble2base_default(nib, (char *)seq, len % 32);
}

static __attribute__((constructor)) CPUType nibble2base_resolve(void) {
#ifdef __ARM_NEON__
    htslibx_nibble2base = nibble2base_neon;
    return CPUType::ARM;
#endif
    return CPUType::UNKNOWN;
}

#endif
#endif  // BUILDING_SIMD_NIBBLE2BASE

struct bam_nibble_data {
    size_t len : 63;
    bool skip_first : 1;
    uint8_t *data;

    bam_nibble_data(const uint8_t *nibbles, size_t from, size_t to)
        : len((to - from)),
          skip_first(from % 2 == 1),
          data(new uint8_t[len / 2 + 1]) {  //(skip_first?1:0)]) {
        memcpy(data, nibbles + (from / 2), (len / 2) + 1);
    }
    bam_nibble_data() : len(0), skip_first(false), data(nullptr) {}
    bam_nibble_data(bam_nibble_data &&other)
        : len(other.len), skip_first(other.skip_first), data(other.data) {
        other.data = nullptr;
    }
    bam_nibble_data(const bam_nibble_data &other) = delete;
    bam_nibble_data &operator=(const bam_nibble_data &other) = delete;
    ~bam_nibble_data() { delete[] data; }
    bam_nibble_data &operator=(bam_nibble_data &&other) {
        if (this != &other) {
            delete[] data;
            data = other.data;
            skip_first = other.skip_first;
            len = other.len;
            other.data = nullptr;
        }
        return *this;
    }
    char *convert_to_ascii(char *buffer) const {
        nibble2base(data, buffer, len + (skip_first ? 1 : 0));
        if (skip_first) {
            buffer[0] = 'N';
        }
        buffer[len + (skip_first ? 1 : 0)] = 0;
        return buffer + len + (skip_first ? 1 : 0);
    }
};

// Uses naive conversions for now can be optimized later if necessary
struct CrumbledSubSeq {
    size_t len : 62;
    uint8_t skip_first_n : 2;
    uint8_t *data;
    struct from_ascii {};
    struct from_nibbles {};
    static constexpr size_t storage_size(size_t bases) noexcept {
        return bases / 4 + 4;
    }
    uint8_t *nibble2crumb(uint8_t *dest, const uint8_t *src, size_t from,
                          size_t to) {
        size_t _from = from & (~3);
        size_t _to = to & (~3);
        if (_to != to) {
            _to += 4;
        }
        for (size_t i = _from; i < _to; i += 4) {
            size_t write_idx = (i - _from) / 4;
            dest[write_idx] = seq_nt16_int[bam_seqi(src, i)] |
                              seq_nt16_int[bam_seqi(src, i + 1)] << 2 |
                              seq_nt16_int[bam_seqi(src, i + 2)] << 4 |
                              seq_nt16_int[bam_seqi(src, i + 3)] << 6;
            // dest[write_idx] = 0xff;
        }
        skip_first_n = from - _from;

        return dest;
    }
    uint8_t *chars2crumb(uint8_t *dest, const char *src, size_t from,
                         size_t to) {
        for (size_t i = from; i < to; i++) {
            size_t write_index = (i - from) / 4;
            size_t write_crumb = (i - from) % 4;
            if (write_crumb == 0) {
                dest[write_index] = 0;
            }
            switch (src[i]) {
                case 'a':
                case 'A':
                    break;
                case 'c':
                case 'C':
                    dest[write_index] |= 1 << (2 * write_crumb);
                    break;
                case 'g':
                case 'G':
                    dest[write_index] |= 2 << (2 * write_crumb);
                    break;
                case 'U':
                case 'u':
                case 'T':
                case 't':
                    dest[write_index] |= 3 << (2 * write_crumb);
                    break;
                default:
                    break;
            }
        }
        return dest;
    }
    char *crumb2base(char *dest, const uint8_t *src, size_t len) const {
        static constexpr const char *lut = "ACGT";
        for (size_t i = skip_first_n; i < len + skip_first_n; ++i) {
            dest[i - skip_first_n] = lut[(src[i / 4] >> (2 * (i % 4))) & 3];
        }
        return dest;
    }
    CrumbledSubSeq(const char *nibbles, size_t from, size_t to,
                   from_ascii) noexcept
        : len((to - from)),
          skip_first_n(0),
          data(new uint8_t[len / 4 + 4]) {  //(skip_first?1:0)]) {
        chars2crumb(data, nibbles, from, to);
    }
    CrumbledSubSeq(const uint8_t *nibbles, size_t from, size_t to,
                   from_nibbles) noexcept
        : len((to - from)),
          skip_first_n(from & 3),
          data(new uint8_t[len / 4 + 4]) {  //(skip_first?1:0)]) {
        nibble2crumb(data, nibbles, from, to);
    }
    CrumbledSubSeq(const uint8_t *nibbles, size_t from, size_t to) noexcept
        : len((to - from)),
          skip_first_n(from & 3),
          data(new uint8_t[len / 4 + 4]) {  //(skip_first?1:0)]) {
        nibble2crumb(data, nibbles, from, to);
    }
    CrumbledSubSeq(size_t length, uint8_t skip) noexcept
        : len(length),
          skip_first_n(skip),
          data(new uint8_t[storage_size(length)]()) {}
    CrumbledSubSeq() noexcept : CrumbledSubSeq(0, 0) {}
    CrumbledSubSeq(CrumbledSubSeq &&other) noexcept
        : len(other.len), skip_first_n(other.skip_first_n), data(other.data) {
        other.data = nullptr;
    }
    CrumbledSubSeq(const CrumbledSubSeq &other) = delete;
    CrumbledSubSeq &operator=(const CrumbledSubSeq &other) = delete;
    ~CrumbledSubSeq() { delete[] data; }
    CrumbledSubSeq &operator=(CrumbledSubSeq &&other) noexcept {
        if (this != &other) {
            delete[] data;
            data = other.data;
            skip_first_n = other.skip_first_n;
            len = other.len;
            other.data = nullptr;
        }
        return *this;
    }
    char *convert_to_ascii(char *buffer) const {
        crumb2base(buffer, data, len);
        buffer[len] = '\0';
        return buffer;
    }
    [[nodiscard]] CrumbledSubSeq reverse_complement() const {
        CrumbledSubSeq rc{len, 0};
        for (size_t i = 0; i < len; ++i) {
            const size_t src_pos = len - 1 - i + skip_first_n;
            const size_t src_byte = src_pos / 4;
            const size_t src_shift = (src_pos % 4) * 2;
            const uint8_t code = (data[src_byte] >> src_shift) & 0x3;
            const uint8_t complement = static_cast<uint8_t>(code ^ 0x3);
            const size_t dst_byte = i / 4;
            const size_t dst_shift = (i % 4) * 2;
            rc.data[dst_byte] |= static_cast<uint8_t>(complement << dst_shift);
        }
        return rc;
    }
};

// https://github.com/ruanjue/wtdbg2/blob/b77c5657c8095412317e4a20fe3668f5bde6b1ac/dna.h#L108
template <class T>
T reverse_complement_kmer(T kmer, int K) {
    uint32_t shift_count = (sizeof(T) * 8) - 2 * K;
    kmer = ~kmer;
    kmer = ((kmer & 0x3333333333333333LLU) << 2) |
           ((kmer & 0xCCCCCCCCCCCCCCCCLLU) >> 2);
    kmer = ((kmer & 0x0F0F0F0F0F0F0F0FLLU) << 4) |
           ((kmer & 0xF0F0F0F0F0F0F0F0LLU) >> 4);
    if constexpr (sizeof(T) == 8) {
        kmer = __builtin_bswap64(kmer);
    } else if constexpr (sizeof(T) == 4) {
        kmer = __builtin_bswap32(kmer);
    } else if constexpr (sizeof(T) == 2) {
        kmer = __builtin_bswap16(kmer);
    }
    kmer = kmer >> shift_count;
    return kmer;
}
// 64-bit packed variant: 32 bases per word, 2 bits per base.
struct CrumbledSubSeq64 {
    size_t len : 59;
    uint8_t skip_first_n : 6;  // offset within first word in bases (0-31)
    uint64_t *data;
    struct from_ascii {};
    struct from_nibbles {};

    static constexpr size_t bases_per_word = 32;
    static constexpr size_t storage_size(size_t bases) noexcept {
        return bases / bases_per_word + 1;  // extra word for padding
    }

    static inline uint8_t base_to_code(char c) {
        switch (c) {
            case 'a':
            case 'A':
                return 0;
            case 'c':
            case 'C':
                return 1;
            case 'g':
            case 'G':
                return 2;
            case 'u':
            case 'U':
            case 't':
            case 'T':
                return 3;
            default:
                return 0;
        }
    }

    template <class FetchCode>
    void pack_range(uint64_t *dest, size_t aligned_from, size_t to,
                    FetchCode &&fetch) {
        const size_t total_bases = to - aligned_from;
        for (size_t i = 0; i < total_bases; ++i) {
            const uint64_t code = fetch(aligned_from + i) & 0x3;
            const size_t bit_index = i * 2;
            dest[bit_index >> 6] |= code << (bit_index & 63);
        }
    }

    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to,
                     from_nibbles) noexcept
        : len(to - from),
          skip_first_n(static_cast<uint8_t>(from & 31)),
          data(nullptr) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        data = new uint64_t[storage_size(total_bases)]();
        pack_range(data, aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(seq_nt16_int[bam_seqi(nibbles, idx)]);
        });
    }

    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to) noexcept
        : CrumbledSubSeq64(nibbles, from, to, from_nibbles{}) {}

    CrumbledSubSeq64(const char *bases, size_t from, size_t to,
                     from_ascii) noexcept
        : len(to - from),
          skip_first_n(static_cast<uint8_t>(from & 31)),
          data(nullptr) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        data = new uint64_t[storage_size(total_bases)]();
        pack_range(data, aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(base_to_code(bases[idx]));
        });
    }

    CrumbledSubSeq64(size_t length, uint8_t skip_bases) noexcept
        : len(length),
          skip_first_n(static_cast<uint8_t>(skip_bases & 63)),
          data(new uint64_t[storage_size(length + skip_bases)]()) {}
    CrumbledSubSeq64() noexcept : CrumbledSubSeq64(0, 0) {}
    CrumbledSubSeq64(CrumbledSubSeq64 &&other) noexcept
        : len(other.len),
          skip_first_n(other.skip_first_n),
          data(other.data) {
        other.data = nullptr;
    }
    CrumbledSubSeq64(const CrumbledSubSeq64 &other) = delete;
    CrumbledSubSeq64 &operator=(const CrumbledSubSeq64 &other) = delete;
    ~CrumbledSubSeq64() { delete[] data; }
    CrumbledSubSeq64 &operator=(CrumbledSubSeq64 &&other) noexcept {
        if (this != &other) {
            delete[] data;
            data = other.data;
            skip_first_n = other.skip_first_n;
            len = other.len;
            other.data = nullptr;
        }
        return *this;
    }

    char *convert_to_ascii(char *buffer) const {
        static constexpr const char *lut = "ACGT";
        size_t base_index = 0;


        auto scalar_emit = [&](size_t count) {
            const size_t bit_base = static_cast<size_t>(skip_first_n) * 2;
            for (size_t i = 0; i < count; ++i) {
                const size_t bit_index = bit_base + (base_index + i) * 2;
                const uint64_t word = data[bit_index >> 6];
                buffer[base_index + i] =
                    lut[(word >> (bit_index & 63)) & 0x3];
            }
            base_index += count;
        };

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
        // Align to word boundary so the NEON path works for any skip_first_n.
        if (len >= bases_per_word) {
            const size_t total = len + skip_first_n;
            const size_t blocks =
                (total + bases_per_word - 1) / bases_per_word;
            const uint64_t *word_ptr = data;
            auto *out = reinterpret_cast<uint8_t *>(buffer);

            static const uint8x8_t ascii_table =
                {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
            const uint8x8_t mask = vdup_n_u8(0x3);

            for (size_t b = 0; b < blocks; ++b) {
                const uint8x8_t packed =
                    vld1_u8(reinterpret_cast<const uint8_t *>(word_ptr));

                const uint8x8_t c0 = vand_u8(packed, mask);
                const uint8x8_t c1 = vand_u8(vshr_n_u8(packed, 2), mask);
                const uint8x8_t c2 = vand_u8(vshr_n_u8(packed, 4), mask);
                const uint8x8_t c3 = vand_u8(vshr_n_u8(packed, 6), mask);

                uint8x8x4_t expanded;
                expanded.val[0] = vtbl1_u8(ascii_table, c0);
                expanded.val[1] = vtbl1_u8(ascii_table, c1);
                expanded.val[2] = vtbl1_u8(ascii_table, c2);
                expanded.val[3] = vtbl1_u8(ascii_table, c3);

                vst4_u8(out, expanded);

                word_ptr++;
                out += bases_per_word;
            }
            buffer[total] = '\0';
            return buffer + skip_first_n;
        }
#elif defined(__SSSE3__) || defined(__AVX2__) // Not tested ever
        // SSSE3 path mirrors the NEON strategy.
        if (len >= bases_per_word) {
            size_t remaining = len;
            const size_t rem = (skip_first_n + base_index) & (bases_per_word - 1);
            const size_t prefix =
                std::min(remaining, (bases_per_word - rem) & (bases_per_word - 1));
            if (prefix) {
                scalar_emit(prefix);
                remaining -= prefix;
            }

            const uint64_t *word_ptr =
                data + ((bit_base + base_index * 2) >> 6);
            auto *out = reinterpret_cast<uint8_t *>(buffer + base_index);

            const __m128i lut = _mm_setr_epi8('A', 'C', 'G', 'T', 'A', 'C', 'G',
                                              'T', 'A', 'C', 'G', 'T', 'A', 'C',
                                              'G', 'T');
            const __m128i mask = _mm_set1_epi8(0x3);

            while (remaining > 0) {
                // Load 8 bytes -> 32 bases spread across bytes.
                __m128i packed =
                    _mm_loadl_epi64(reinterpret_cast<const __m128i *>(word_ptr));

                __m128i c0 = _mm_and_si128(packed, mask);
                __m128i c1 = _mm_and_si128(_mm_srli_epi16(packed, 2), mask);
                __m128i c2 = _mm_and_si128(_mm_srli_epi16(packed, 4), mask);
                __m128i c3 = _mm_and_si128(_mm_srli_epi16(packed, 6), mask);

                c0 = _mm_shuffle_epi8(lut, c0);
                c1 = _mm_shuffle_epi8(lut, c1);
                c2 = _mm_shuffle_epi8(lut, c2);
                c3 = _mm_shuffle_epi8(lut, c3);

                // Interleave to restore base order.
                __m128i t0lo = _mm_unpacklo_epi8(c0, c1);
                __m128i t0hi = _mm_unpackhi_epi8(c0, c1);
                __m128i t1lo = _mm_unpacklo_epi8(c2, c3);
                __m128i t1hi = _mm_unpackhi_epi8(c2, c3);

                __m128i o0 = _mm_unpacklo_epi16(t0lo, t1lo);
                __m128i o1 = _mm_unpackhi_epi16(t0lo, t1lo);
                __m128i o2 = _mm_unpacklo_epi16(t0hi, t1hi);
                __m128i o3 = _mm_unpackhi_epi16(t0hi, t1hi);

                _mm_storeu_si128(reinterpret_cast<__m128i *>(out), o0);
                _mm_storeu_si128(reinterpret_cast<__m128i *>(out + 16), o1);
                _mm_storeu_si128(reinterpret_cast<__m128i *>(out + 32), o2);
                _mm_storeu_si128(reinterpret_cast<__m128i *>(out + 48), o3);

                size_t step = std::min(remaining, bases_per_word);
                word_ptr++;
                out += step;
                base_index += step;
                remaining -= step;
            }
            buffer[len] = '\0';
            return buffer + skip_first_n;
        }
#endif
        scalar_emit(len);
        buffer[len] = '\0';
        return buffer;
    }

    [[nodiscard]] CrumbledSubSeq64 reverse_complement_old() const {
        CrumbledSubSeq64 rc{len, 0};
        for (size_t i = 0; i < len; ++i) {
            const size_t src_bit = (skip_first_n + (len - 1 - i)) * 2;
            const uint64_t code =
                (data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
            const uint64_t complement = code ^ 0x3;
            const size_t dst_bit = i * 2;
            rc.data[dst_bit >> 6] |= complement << (dst_bit & 63);
        }
        return rc;
    }
    [[nodiscard]] CrumbledSubSeq64 reverse_complement() const {
        CrumbledSubSeq64 rc{len, 31};
        if(len < 32){
            // Fix the skip_first
            uint64_t word = (data[0] >> skip_first_n * 2);
            if(storage_size(len+skip_first_n) > 1){
                word |= (data[1] << (64 - skip_first_n * 2));
            }

            rc.data[0] = reverse_complement_kmer(word, len);
            rc.skip_first_n = 0;
            return rc;
        }
        size_t num_blocks = storage_size(len + skip_first_n);
        for(size_t i = 0 ; i < num_blocks; ++i){
            rc.data[num_blocks - i - 1] = reverse_complement_kmer(data[i], 32);
        }
        rc.skip_first_n =  (num_blocks * 32 - len) - skip_first_n;
        return rc;
    }
};
#endif
