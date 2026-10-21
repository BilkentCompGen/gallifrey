#ifndef NIBBLE2BASE_HPP
#define NIBBLE2BASE_HPP

#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <gtl/vector.hpp>
#include <iterator>
#include <memory>
#include <span>
#include <vector>

#include "htslib/hts.h"
#include "kmer_bits.h"

static const uint8_t decode4_tbl[16] = {
    /* 0 */ 0,
    /* 1 */ 0,
    /* 2 */ 1,
    /* 3 */ 0,
    /* 4 */ 2,
    /* 5 */ 0,
    /* 6 */ 0,
    /* 7 */ 0,
    /* 8 */ 3,
    /* 9 */ 0,
    /* A */ 0,
    /* B */ 0,
    /* C */ 0,
    /* D */ 0,
    /* E */ 0,
    /* F */ 0};

#if (defined(__ARM_NEON) || defined(__ARM_NEON__))

#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>

inline void onehotto2bitDNA_neon(
    const uint8_t *__restrict data, size_t len, uint8_t *__restrict low,
    uint8_t *__restrict high) {
    const uint8x16_t tbl = vld1q_u8(decode4_tbl);
    const uint8x16_t mask = vdupq_n_u8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        uint8x16_t v = vld1q_u8(data + i);
        uint8x16_t lo = vandq_u8(v, mask);
        uint8x16_t hi = vshrq_n_u8(v, 4);
        uint8x16_t dlo = vqtbl1q_u8(tbl, lo);
        uint8x16_t dhi = vqtbl1q_u8(tbl, hi);
        vst1q_u8(low + i, dlo);
        vst1q_u8(high + i, dhi);
    }
    for (; i < len; ++i) {
        uint8_t b = data[i];
        low[i] = decode4_tbl[b & 0x0F];
        high[i] = decode4_tbl[b >> 4];
    }
}

#elif defined(__SSSE3__)

#include <tmmintrin.h>
#include <stddef.h>
#include <stdint.h>

inline void onehotto2bitDNA_sse(
    const uint8_t *__restrict data, size_t len, uint8_t *__restrict low,
    uint8_t *__restrict high) {
    const __m128i tbl = _mm_loadu_si128(reinterpret_cast<const __m128i *>(decode4_tbl));
    const __m128i mask = _mm_set1_epi8(0x0F);
    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i *>(data + i));
        __m128i lo = _mm_and_si128(v, mask);
        __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), mask);
        __m128i dlo = _mm_shuffle_epi8(tbl, lo);
        __m128i dhi = _mm_shuffle_epi8(tbl, hi);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(low + i), dlo);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(high + i), dhi);
    }
    for (; i < len; ++i) {
        uint8_t b = data[i];
        low[i] = decode4_tbl[b & 0x0F];
        high[i] = decode4_tbl[b >> 4];
    }
}

#endif

inline void onehotto2bitDNA_vector(const uint8_t *__restrict data, size_t len,
                                   uint8_t *__restrict low,
                                   uint8_t *__restrict high) {
#if (defined(__ARM_NEON) || defined(__ARM_NEON__))
    return onehotto2bitDNA_neon(data, len, low, high);
#elif defined(__SSSE3__)
    return onehotto2bitDNA_sse(data, len, low, high);
#else
    for (size_t i = 0; i < len; ++i) {
        uint8_t b = data[i];
        low[i] = decode4_tbl[b & 0x0F];
        high[i] = decode4_tbl[b >> 4];
    }
#endif
}

// Bump allocator for 2-bit words: zeroed storage carved from large blocks,
// all freed together by release() or the destructor. Not thread-safe; each
// inserting thread gets its own (dbg_index_base::fragment_arena). Aligned to
// the x86-64 cache line so that arenas of concurrent workers, kept side by
// side, never share a line.
class alignas(64) word_arena {
   public:
    word_arena() = default;
    // Moves hand the blocks over; the words stay where they are, so
    // fragments borrowing them stay valid under the new owner.
    word_arena(word_arena &&other) noexcept { take(other); }
    word_arena &operator=(word_arena &&other) noexcept {
        if (this != &other) take(other);
        return *this;
    }
    [[nodiscard]] uint64_t *allocate(size_t words) {
        if (words > left_) {
            left_ = std::max(words, kBlockWords);
            blocks_.push_back(std::make_unique<uint64_t[]>(left_));
            cur_ = blocks_.back().get();
        }
        uint64_t *out = cur_;
        cur_ += words;
        left_ -= words;
        return out;
    }
    // Frees every block; storage handed out before is invalid afterwards.
    void release() noexcept {
        std::vector<std::unique_ptr<uint64_t[]>>().swap(blocks_);
        cur_ = nullptr;
        left_ = 0;
    }

   private:
    void take(word_arena &other) noexcept {
        blocks_ = std::move(other.blocks_);
        cur_ = std::exchange(other.cur_, nullptr);
        left_ = std::exchange(other.left_, 0);
        other.blocks_.clear();
    }
    // 1 MiB blocks: thousands of fragments per malloc, while the unused tail
    // of one block per thread stays negligible against a graph.
    static constexpr size_t kBlockWords = size_t{1} << 17;
    std::vector<std::unique_ptr<uint64_t[]>> blocks_;
    uint64_t *cur_ = nullptr;
    size_t left_ = 0;
};

// 64-bit packed variant: 32 bases per word, 2 bits per base.
struct CrumbledSubSeq64 {
    size_t len : 57;
    int skip_first_n : 6;  // offset within first word in bases (0-31)
    // The words belong to a word_arena, not to this fragment: destruction
    // and assignment leave them alone. Moves carry the flag along.
    size_t borrowed : 1 = 0;
    uint64_t *data;
    uint64_t origin_node;
    struct from_ascii {};
    struct from_nibbles {};

    static constexpr size_t bases_per_word = 32;
    static constexpr size_t storage_size(size_t bases) noexcept {
        return bases / bases_per_word + 1;  // extra word for padding
    }
    template <std::forward_iterator It>
    static size_t combined_length(It first, It last) {
        size_t total = 0;
        for (auto it = first; it != last; ++it) {
            total += static_cast<size_t>((*it).len);
        }
        return total;
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
                     uint64_t origin_node, from_nibbles) noexcept
        : len(to - from),
          skip_first_n(static_cast<uint8_t>(from & 31)),
          data(nullptr),
          origin_node(origin_node) {
        const size_t aligned_from = from & ~static_cast<size_t>(31);
        const size_t total_bases = (to - aligned_from);
        data = new uint64_t[storage_size(total_bases)]();
        pack_range(data, aligned_from, to, [&](size_t idx) {
            return static_cast<uint64_t>(seq_nt16_int[bam_seqi(nibbles, idx)]);
        });
    }

    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin_node) noexcept
        : CrumbledSubSeq64(nibbles, from, to, origin_node, from_nibbles{}) {}

    CrumbledSubSeq64(const char *bases, size_t from, size_t to,
                     uint64_t origin_node, from_ascii) noexcept
        : len(to - from),
          skip_first_n(static_cast<uint8_t>(from & 31)),
          data(nullptr),
          origin_node(origin_node) {
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
          data(new uint64_t[storage_size(length + skip_bases)]()),
          origin_node(0) {}
    // As above, with the words from `arena` (the fragment borrows them).
    CrumbledSubSeq64(size_t length, uint8_t skip_bases,
                     word_arena &arena) noexcept
        : len(length),
          skip_first_n(static_cast<uint8_t>(skip_bases & 63)),
          borrowed(1),
          data(arena.allocate(storage_size(length + skip_bases))),
          origin_node(0) {}
    CrumbledSubSeq64(const uint8_t *nibbles, size_t from, size_t to,
                     uint64_t origin_node, word_arena &arena) noexcept
        : CrumbledSubSeq64(to - from, static_cast<uint8_t>(from & 31), arena) {
        this->origin_node = origin_node;
        pack_range(data, from & ~static_cast<size_t>(31), to, [&](size_t idx) {
            return static_cast<uint64_t>(seq_nt16_int[bam_seqi(nibbles, idx)]);
        });
    }
    CrumbledSubSeq64() noexcept : CrumbledSubSeq64(0, 0) {}
    CrumbledSubSeq64(CrumbledSubSeq64 &&other) noexcept
        : len(other.len),
          skip_first_n(other.skip_first_n),
          borrowed(other.borrowed),
          data(other.data),
          origin_node(other.origin_node) {
        other.data = nullptr;
        other.origin_node = 0;
    }

    template <std::forward_iterator It>
    CrumbledSubSeq64(It first, It last)
        : len(combined_length(first, last)),
          skip_first_n(0),
          data(new uint64_t[storage_size(len + 32)]()),
          origin_node(0) {
        size_t dest_pos = 0;
        for (auto it = first; it != last; ++it) {
            const auto &frag = *it;
            for (size_t b = 0; b < frag.len; ++b) {
                const size_t src_bit =
                    (static_cast<size_t>(frag.skip_first_n) + b) * 2;
                const uint64_t code =
                    (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                const size_t dst_bit = dest_pos * 2;
                data[dst_bit >> 6] |= code << (dst_bit & 63);
                ++dest_pos;
            }
            origin_node = frag.origin_node;
        }
    }

    CrumbledSubSeq64(std::span<const CrumbledSubSeq64> subseqs,
                     size_t start_idx, size_t from, size_t to)
        : len{to - from},
          skip_first_n{0},
          data{new uint64_t[storage_size(to - from + 32)]()},
          origin_node{0} {
        size_t sub_len = (from < to) ? (to - from) : 0;
        if (sub_len == 0 || subseqs.empty()) {
            return;
        }

        size_t cumulative = 0;
        size_t dest_pos = 0;

        for (size_t i = start_idx; i < subseqs.size() && dest_pos < sub_len; ++i) {
            const auto &frag = subseqs[i];
            size_t frag_start = cumulative;
            size_t frag_end = cumulative + frag.len;

            // Compute the overlap of [frag_start, frag_end) with [from, to)
            size_t overlap_start = std::max(frag_start, from);
            size_t overlap_end = std::min(frag_end, to);

            if (overlap_start < overlap_end) {
                // Offset within this fragment
                size_t local_from = overlap_start - frag_start;
                size_t copy_len = overlap_end - overlap_start;

                // Copy bases from frag to result
                for (size_t b = 0; b < copy_len; ++b) {
                    size_t src_bit = (frag.skip_first_n + local_from + b) * 2;
                    uint64_t code =
                        (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                    size_t dst_bit = dest_pos * 2;
                    data[dst_bit >> 6] |= code << (dst_bit & 63);
                    ++dest_pos;
                }
                origin_node = frag.origin_node;
            }
            cumulative += frag.len;
        }
    };
    CrumbledSubSeq64(std::span<const CrumbledSubSeq64> subseqs, size_t from,
                     size_t to)
        : len{to - from},
          skip_first_n{0},
          data{new uint64_t[storage_size(to - from + 32)]()},
          origin_node{0} {
        size_t sub_len = (from < to) ? (to - from) : 0;
        if (sub_len == 0 || subseqs.empty()) {
            return;
        }

        size_t cumulative = 0;
        size_t dest_pos = 0;

        for (size_t i = 0; i < subseqs.size() && dest_pos < sub_len; ++i) {
            const auto &frag = subseqs[i];
            size_t frag_start = cumulative;
            size_t frag_end = cumulative + frag.len;

            // Compute the overlap of [frag_start, frag_end) with [from, to)
            size_t overlap_start = std::max(frag_start, from);
            size_t overlap_end = std::min(frag_end, to);

            if (overlap_start < overlap_end) {
                // Offset within this fragment
                size_t local_from = overlap_start - frag_start;
                size_t copy_len = overlap_end - overlap_start;

                // Copy bases from frag to result
                for (size_t b = 0; b < copy_len; ++b) {
                    size_t src_bit = (frag.skip_first_n + local_from + b) * 2;
                    uint64_t code =
                        (frag.data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
                    size_t dst_bit = dest_pos * 2;
                    data[dst_bit >> 6] |= code << (dst_bit & 63);
                    ++dest_pos;
                }
                origin_node = frag.origin_node;
            }
            cumulative += frag.len;
        }
    };
    CrumbledSubSeq64(size_t len, int skip_first_n, uint64_t *data,
                     uint64_t origin_node) noexcept
        : len{len},
          skip_first_n{skip_first_n},
          data(data),
          origin_node(origin_node) {}
    static CrumbledSubSeq64 copy_of(const CrumbledSubSeq64 &other) {
        const size_t words = storage_size(other.len + other.skip_first_n);
        uint64_t *copied_data = new uint64_t[words]();
        std::memcpy(copied_data, other.data, words * sizeof(uint64_t));
        return CrumbledSubSeq64(other.len, other.skip_first_n, copied_data,
                                other.origin_node);
    }
    CrumbledSubSeq64(const CrumbledSubSeq64 &other) = delete;
    CrumbledSubSeq64 &operator=(const CrumbledSubSeq64 &other) = delete;
    ~CrumbledSubSeq64() {
        if (!borrowed) delete[] data;
    }
    CrumbledSubSeq64 &operator=(CrumbledSubSeq64 &&other) noexcept {
        if (this != &other) {
            if (!borrowed) delete[] data;
            data = other.data;
            skip_first_n = other.skip_first_n;
            borrowed = other.borrowed;
            len = other.len;
            origin_node = other.origin_node;
            other.data = nullptr;
            other.origin_node = 0;
        }
        return *this;
    }

    // 2-bit code (A=0, C=1, G=2, T=3) of base i.
    [[nodiscard]] uint8_t code(size_t i) const noexcept {
        const size_t bit = (static_cast<size_t>(skip_first_n) + i) * 2;
        return (data[bit >> 6] >> (bit & 63)) & 0x3;
    }

    // Base codes of [from, to), one byte per base: out[0] is base `from`.
    void unpack_bytes(size_t from, size_t to, uint8_t *out) const noexcept {
        size_t i = from;
        for (; i + bases_per_word <= to; i += bases_per_word)
            unpack_word(word_at(i), out + (i - from));
        for (; i < to; ++i) out[i - from] = code(i);
    }

    // The same codes last-first (out[0] is base to - 1), not complemented.
    void unpack_bytes_reversed(size_t from, size_t to,
                               uint8_t *out) const noexcept {
        size_t i = to;
        for (; i - from >= bases_per_word; i -= bases_per_word) {
            unpack_word_reversed(word_at(i - bases_per_word), out);
            out += bases_per_word;
        }
        for (; i > from; --i) *out++ = code(i - 1);
    }

    // The 32 bases from base i as one word (base i in the low bits); needs
    // i + 32 <= len. The next word is read only when the window spans it,
    // so it lies within storage_size(len + skip_first_n).
    [[nodiscard]] uint64_t word_at(size_t i) const noexcept {
        const size_t pos = static_cast<size_t>(skip_first_n) + i;
        const size_t w = pos / bases_per_word;
        const unsigned shift = (pos % bases_per_word) * 2;
        if (shift == 0) return data[w];
        return (data[w] >> shift) | (data[w + 1] << (64 - shift));
    }

    // The 32 codes of one word, one byte each. Each byte holds 4 bases low
    // bits first: splitting every byte into its low and high nibble and
    // interleaving them gives one nibble (2 bases) per byte in base order;
    // splitting and interleaving the nibbles' 2-bit halves gives the codes.
    static void unpack_word(uint64_t word, uint8_t *out) noexcept {
#if (defined(__ARM_NEON) || defined(__ARM_NEON__))
        const uint8x16x2_t codes = unpack_word_neon(word);
        vst1q_u8(out, codes.val[0]);
        vst1q_u8(out + 16, codes.val[1]);
#elif defined(__SSSE3__)
        __m128i lo, hi;
        unpack_word_sse(word, lo, hi);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out), lo);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out + 16), hi);
#else
        for (size_t k = 0; k < bases_per_word; ++k)
            out[k] = (word >> (2 * k)) & 0x3;
#endif
    }

    // unpack_word with the 32 bytes in reverse order.
    static void unpack_word_reversed(uint64_t word, uint8_t *out) noexcept {
#if (defined(__ARM_NEON) || defined(__ARM_NEON__))
        const uint8x16x2_t codes = unpack_word_neon(word);
        const uint8x16_t lo = vrev64q_u8(codes.val[0]);
        const uint8x16_t hi = vrev64q_u8(codes.val[1]);
        vst1q_u8(out, vcombine_u8(vget_high_u8(hi), vget_low_u8(hi)));
        vst1q_u8(out + 16, vcombine_u8(vget_high_u8(lo), vget_low_u8(lo)));
#elif defined(__SSSE3__)
        __m128i lo, hi;
        unpack_word_sse(word, lo, hi);
        const __m128i rev =
            _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out),
                         _mm_shuffle_epi8(hi, rev));
        _mm_storeu_si128(reinterpret_cast<__m128i *>(out + 16),
                         _mm_shuffle_epi8(lo, rev));
#else
        for (size_t k = 0; k < bases_per_word; ++k)
            out[k] = (word >> (2 * (bases_per_word - 1 - k))) & 0x3;
#endif
    }

#if (defined(__ARM_NEON) || defined(__ARM_NEON__))
    // Codes of bases 0-15 (val[0]) and 16-31 (val[1]) of one word.
    static uint8x16x2_t unpack_word_neon(uint64_t word) noexcept {
        const uint8x8_t bytes = vcreate_u8(word);
        const uint8x8x2_t nib =
            vzip_u8(vand_u8(bytes, vdup_n_u8(0x0F)), vshr_n_u8(bytes, 4));
        const uint8x16_t nibbles = vcombine_u8(nib.val[0], nib.val[1]);
        return vzipq_u8(vandq_u8(nibbles, vdupq_n_u8(0x3)),
                        vshrq_n_u8(nibbles, 2));
    }
#elif defined(__SSSE3__)
    // Codes of bases 0-15 (lo) and 16-31 (hi) of one word.
    static void unpack_word_sse(uint64_t word, __m128i &lo,
                                __m128i &hi) noexcept {
        const __m128i bytes = _mm_cvtsi64_si128(static_cast<long long>(word));
        const __m128i nibble_mask = _mm_set1_epi8(0x0F);
        const __m128i nibbles = _mm_unpacklo_epi8(
            _mm_and_si128(bytes, nibble_mask),
            _mm_and_si128(_mm_srli_epi16(bytes, 4), nibble_mask));
        const __m128i code_mask = _mm_set1_epi8(0x3);
        const __m128i low2 = _mm_and_si128(nibbles, code_mask);
        const __m128i high2 =
            _mm_and_si128(_mm_srli_epi16(nibbles, 2), code_mask);
        lo = _mm_unpacklo_epi8(low2, high2);
        hi = _mm_unpackhi_epi8(low2, high2);
    }
#endif

    // The bases as letters, then '\0' at buffer[len]; needs len + 1 bytes.
    // The codes are unpacked into the buffer and turned into letters in place.
    char *convert_to_ascii(char *buffer) const {
        uint8_t *codes = reinterpret_cast<uint8_t *>(buffer);
        unpack_bytes(0, len, codes);
        for (size_t i = 0; i < len; ++i) buffer[i] = kBaseLetters[codes[i]];
        buffer[len] = '\0';
        return buffer;
    }

    [[nodiscard]] CrumbledSubSeq64 reverse_complement() const {
        CrumbledSubSeq64 rc{len, 0};
        for (size_t i = 0; i < len; ++i) {
            const size_t src_bit = (skip_first_n + (len - 1 - i)) * 2;
            const uint64_t code = (data[src_bit >> 6] >> (src_bit & 63)) & 0x3;
            const uint64_t complement = code ^ 0x3;
            const size_t dst_bit = i * 2;
            rc.data[dst_bit >> 6] |= complement << (dst_bit & 63);
        }
        return rc;
    }
    [[nodiscard]] CrumbledSubSeq64 reverse_complement_old() const {
        assert(skip_first_n >= 0);
        CrumbledSubSeq64 rc{len, 31};
        if (len < 32) {
            // Fix the skip_first
            uint64_t word = (data[0] >> skip_first_n * 2);
            if (storage_size(len + skip_first_n) > 1) {
                word |= (data[1] << (64 - skip_first_n * 2));
            }

            rc.data[0] = reverse_complement_kmer(word, len);
            rc.skip_first_n = 0;
            return rc;
        }
        size_t num_blocks = storage_size(len + skip_first_n);
        for (size_t i = 0; i < num_blocks - 1; ++i) {
            rc.data[num_blocks - i - 1] = reverse_complement_kmer(data[i], 32);
        }
        if ((len + skip_first_n) % 32 != 0)
            rc.data[0] = reverse_complement_kmer(data[num_blocks - 1],
                                                 (len + skip_first_n) % 32);
        rc.skip_first_n = (num_blocks * 32 - len) - skip_first_n;
        return rc;
    }
};
static_assert(sizeof(CrumbledSubSeq64) == 24,
              "the borrowed flag shares the length word");
struct seqveciterator {
    size_t idx = 0;
    size_t bidx = 0;
    size_t bases_before = 0;
    void move_only(const gtl::vector<CrumbledSubSeq64> &sseqs,
                   uint64_t until) {
        size_t new_idx = idx;
        size_t new_bidx = bidx;
        while (new_idx < sseqs.size() &&
               bases_before + (sseqs[new_idx].len - new_bidx) < until) {
            bases_before += sseqs[new_idx].len - new_bidx;
            new_bidx = 0;
            new_idx++;
        }
        new_bidx += until - bases_before;
        idx = new_idx;
        bidx = new_bidx;
        bases_before = until;
    }
    CrumbledSubSeq64 get_and_move(
        const gtl::vector<CrumbledSubSeq64> &sseqs, uint64_t until) {
        // Save bases_before before the loop mutates it, so we can
        // compute the full extraction range for the constructor.
        const size_t saved_bases_before = bases_before;
        size_t new_idx = idx;
        size_t new_bidx = bidx;
        while (new_idx < sseqs.size() &&
               bases_before + (sseqs[new_idx].len - new_bidx) < until) {
            bases_before += sseqs[new_idx].len - new_bidx;
            new_bidx = 0;
            new_idx++;
        }
        new_bidx += until - bases_before;
        // 'from' and 'to' are positions relative to the concatenation
        // starting at start_idx (idx).  We start at bidx within the
        // first fragment and need (until - saved_bases_before) bases.
        CrumbledSubSeq64 seq{sseqs, idx, bidx,
                                  bidx + (until - saved_bases_before)};
        idx = new_idx;
        bidx = new_bidx;
        bases_before = until;
        return seq;
    }
};

#endif
