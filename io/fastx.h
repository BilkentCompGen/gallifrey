#ifndef _GALLIFREY_FASTX_H_
#define _GALLIFREY_FASTX_H_

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "../nt16.h"
#include "input_stream.h"

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#elif defined(__SSSE3__)
#include <tmmintrin.h>
#endif

// Native FASTA/FASTQ input path.
//
// Reading FASTA/FASTQ through a generic sequence record would pay for a copy
// per line and a quality array the assembler never looks at. The assembler
// only needs the packed 4-bit sequence and its length, so this reader scans
// blocks with memchr and packs bases straight into the record buffer in the
// nt16 nibble layout of the BAM SEQ field (A=1, C=2, G=4, T=8, other=15, high
// nibble first). Everything downstream -- onehot_to_rle(), CrumbledSubSeq64's
// from_nibbles ctor -- is therefore byte-identical to the BAM path.

// Recycled sequence record of every reader. Buffers grow only.
struct seq_record {
    uint8_t *nibbles = nullptr;
    uint32_t len = 0;  // bases
    uint32_t cap = 0;  // allocated bytes

    void reserve_bases(uint32_t bases) {
        // +8 bytes of slack so nibble writers/readers that round up to a word
        // never touch unallocated memory.
        const uint32_t need = bases / 2 + 8;
        if (need <= cap) return;
        uint32_t ncap = cap ? cap : 8192;
        while (ncap < need) ncap += ncap / 2;
        uint8_t *p = static_cast<uint8_t *>(realloc(nibbles, ncap));
        if (!p) {
            fprintf(stderr, "seq_record: out of memory (%u bytes)\n", ncap);
            abort();
        }
        nibbles = p;
        cap = ncap;
    }
};

inline seq_record *seq_record_init() { return new seq_record(); }

inline void seq_record_destroy(seq_record *r) {
    if (!r) return;
    free(r->nibbles);
    delete r;
}

// `Source` hands out the bytes in blocks: bool next(const uint8_t *&, size_t &)
// and const std::string &error() (empty after a clean end), as input_stream.
template <class Source>
class fastx_reader {
    // kseq record grammar, driven incrementally so a record may straddle any
    // number of blocks: header line -> sequence lines -> ('+' line -> quality
    // bytes) for FASTQ. A record ends at a line starting with '>' or '@', or
    // once `len` quality bytes have been consumed.
    enum class state { start, header, seq, plus, qual };

    Source &src_;
    const uint8_t *buf_ = nullptr;  // current block of src_
    size_t pos_ = 0;
    size_t end_ = 0;
    state st_ = state::start;
    bool line_start_ = true;
    bool eof_ = false;
    bool failed_ = false;
    uint32_t qual_got_ = 0;
    std::string name_;
    bool name_done_ = false;
    // Shape of the record being read: its sequence lines, whether it has a
    // '+' line, whether its quality spans a line break, and whether the input
    // ended before its quality was complete.
    uint32_t seq_lines_ = 0;
    bool fastq_ = false;
    bool qual_wrapped_ = false;
    bool truncated_ = false;

    // Takes the next block. Returns false at EOF or on error.
    bool refill() {
        if (eof_) return false;
        const uint8_t *p = nullptr;
        size_t n = 0;
        if (!src_.next(p, n)) {
            eof_ = true;
            if (!src_.error().empty()) failed_ = true;
            return false;
        }
        buf_ = p;
        pos_ = 0;
        end_ = n;
        return true;
    }

    static inline void pack_one(uint8_t *dst, uint32_t &len, uint8_t c) {
        const uint8_t nib = seq_nt16_table[c];
        if (len & 1) {
            dst[len >> 1] = static_cast<uint8_t>(dst[len >> 1] | nib);
        } else {
            dst[len >> 1] = static_cast<uint8_t>(nib << 4);
        }
        ++len;
    }

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    using vec_t = uint8x16_t;
    using half_t = uint8x8_t;
    static inline vec_t vload(const uint8_t *p) { return vld1q_u8(p); }
    // 0xFF in every lane holding one of ACGTacgt, 0x00 elsewhere.
    static inline vec_t acgt_mask(vec_t v) {
        const vec_t lower = vorrq_u8(v, vdupq_n_u8(0x20));
        return vorrq_u8(vorrq_u8(vceqq_u8(lower, vdupq_n_u8('a')),
                                 vceqq_u8(lower, vdupq_n_u8('c'))),
                        vorrq_u8(vceqq_u8(lower, vdupq_n_u8('g')),
                                 vceqq_u8(lower, vdupq_n_u8('t'))));
    }
    static inline bool all_acgt(vec_t v) {
        return vminvq_u8(acgt_mask(v)) == 0xFF;
    }
    static inline bool all_acgt(vec_t a, vec_t b) {
        return vminvq_u8(vandq_u8(acgt_mask(a), acgt_mask(b))) == 0xFF;
    }
    // 16 ASCII ACGT bases -> 8 packed nt16 bytes. (c >> 1) & 3 maps A/C/T/G
    // (either case) to 0/1/2/3; the table turns that into the one-hot nibble.
    static inline half_t pack16(vec_t v) {
        static const uint8_t code2nib[16] = {1, 2, 8, 4, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0};
        const vec_t nib = vqtbl1q_u8(vld1q_u8(code2nib),
                                     vandq_u8(vshrq_n_u8(v, 1), vdupq_n_u8(3)));
        const vec_t hi = vuzp1q_u8(nib, nib);  // even bases
        const vec_t lo = vuzp2q_u8(nib, nib);  // odd bases
        return vget_low_u8(vorrq_u8(vshlq_n_u8(hi, 4), lo));
    }
    static inline void store8(uint8_t *dst, half_t v) { vst1_u8(dst, v); }
    static inline void store16(uint8_t *dst, half_t a, half_t b) {
        vst1q_u8(dst, vcombine_u8(a, b));
    }
#define GALLIFREY_FASTX_SIMD 1
#elif defined(__SSSE3__)
    using vec_t = __m128i;
    using half_t = __m128i;  // low 8 bytes significant
    static inline vec_t vload(const uint8_t *p) {
        return _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
    }
    static inline vec_t acgt_mask(vec_t v) {
        const vec_t lower = _mm_or_si128(v, _mm_set1_epi8(0x20));
        return _mm_or_si128(
            _mm_or_si128(_mm_cmpeq_epi8(lower, _mm_set1_epi8('a')),
                         _mm_cmpeq_epi8(lower, _mm_set1_epi8('c'))),
            _mm_or_si128(_mm_cmpeq_epi8(lower, _mm_set1_epi8('g')),
                         _mm_cmpeq_epi8(lower, _mm_set1_epi8('t'))));
    }
    static inline bool all_acgt(vec_t v) {
        return _mm_movemask_epi8(acgt_mask(v)) == 0xFFFF;
    }
    static inline bool all_acgt(vec_t a, vec_t b) {
        return _mm_movemask_epi8(_mm_and_si128(acgt_mask(a), acgt_mask(b))) ==
               0xFFFF;
    }
    static inline half_t pack16(vec_t v) {
        static const uint8_t code2nib[16] = {1, 2, 8, 4, 0, 0, 0, 0,
                                             0, 0, 0, 0, 0, 0, 0, 0};
        const vec_t lut =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(code2nib));
        // Gather even / odd lanes into the low 8 bytes.
        const vec_t even_idx = _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, -1, -1,
                                             -1, -1, -1, -1, -1, -1);
        const vec_t odd_idx = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1, -1,
                                            -1, -1, -1, -1, -1, -1);
        const vec_t nib = _mm_shuffle_epi8(
            lut, _mm_and_si128(_mm_srli_epi16(v, 1), _mm_set1_epi8(3)));
        const vec_t hi = _mm_shuffle_epi8(nib, even_idx);
        const vec_t lo = _mm_shuffle_epi8(nib, odd_idx);
        // Nibbles are <= 8, so the 16-bit shift cannot carry across bytes.
        return _mm_or_si128(_mm_slli_epi16(hi, 4), lo);
    }
    static inline void store8(uint8_t *dst, half_t v) {
        _mm_storel_epi64(reinterpret_cast<__m128i *>(dst), v);
    }
    static inline void store16(uint8_t *dst, half_t a, half_t b) {
        _mm_storeu_si128(reinterpret_cast<__m128i *>(dst),
                         _mm_unpacklo_epi64(a, b));
    }
#define GALLIFREY_FASTX_SIMD 1
#endif

    // Packs [p, q) of ASCII bases onto the end of `rec`, dropping any byte at
    // or below ' ' (CR, LF, spaces, tabs). This is the hot loop of the reader:
    // it, not the I/O, bounds throughput. Plain ACGT runs go 32 bases per
    // iteration through SIMD; any window containing whitespace, N or IUPAC
    // codes is packed one base at a time and the vector loop resumes after it.
    static void pack_bases(seq_record &rec, const uint8_t *p,
                           const uint8_t *q) {
        rec.reserve_bases(rec.len + static_cast<uint32_t>(q - p));
        uint8_t *dst = rec.nibbles;
        uint32_t len = rec.len;
#ifdef GALLIFREY_FASTX_SIMD
        // Vector stores need an even base offset (whole output bytes); pack
        // one base scalar to get there.
        while ((len & 1) && p != q) {
            if (*p > ' ') pack_one(dst, len, *p);
            ++p;
        }
        auto scalar16 = [&]() {
            const uint8_t *stop = p + 16;
            for (; p != stop; ++p) {
                if (*p > ' ') pack_one(dst, len, *p);
            }
        };
        while ((len & 1) == 0 && q - p >= 32) {
            const vec_t a = vload(p);
            const vec_t b = vload(p + 16);
            if (all_acgt(a, b)) {
                store16(dst + (len >> 1), pack16(a), pack16(b));
                len += 32;
                p += 32;
            } else if (all_acgt(a)) {
                store8(dst + (len >> 1), pack16(a));
                len += 16;
                p += 16;
                scalar16();
            } else {
                scalar16();
            }
        }
        while ((len & 1) == 0 && q - p >= 16) {
            const vec_t a = vload(p);
            if (all_acgt(a)) {
                store8(dst + (len >> 1), pack16(a));
                len += 16;
                p += 16;
            } else {
                scalar16();
            }
        }
#endif
        for (; p != q; ++p) {
            const uint8_t c = *p;
            if (c <= ' ') continue;
            pack_one(dst, len, c);
        }
        rec.len = len;
    }

   public:
    // Reads FASTA/FASTQ from `src` (opened by the caller and not owned).
    // Record boundaries are handled incrementally, so any block size of `src`
    // parses identically.
    explicit fastx_reader(Source &src) : src_(src) {}
    fastx_reader(const fastx_reader &) = delete;
    fastx_reader &operator=(const fastx_reader &) = delete;

    bool failed() const { return failed_; }

    // Why the input could not be read or decompressed; empty when it was
    // malformed FASTA/FASTQ, or read cleanly.
    const std::string &stream_error() const { return src_.error(); }

    // Name of the record returned by the last successful next(): the header
    // line truncated at the first whitespace.
    const std::string &name() const { return name_; }

    // Shape of the record returned by the last successful next(): FASTQ
    // (it has a '+' line), and FASTQ in exactly four lines whose quality
    // was read to its full length.
    bool fastq_record() const { return fastq_; }
    bool four_line_fastq_record() const {
        return fastq_ && seq_lines_ == 1 && !qual_wrapped_ && !truncated_;
    }
    // Whether that record is FASTQ whose input ended before its quality was
    // complete.
    bool truncated_record() const { return truncated_; }

    // Fills `rec` with the next record's packed sequence. Returns false at EOF.
    bool next(seq_record &rec) {
        rec.len = 0;
        while (true) {
            if (pos_ == end_ && !refill()) {
                // End of stream: emit whatever the last record accumulated,
                // unless the stream broke off in it.
                truncated_ = st_ == state::plus ||
                             (st_ == state::qual && qual_got_ < rec.len);
                st_ = state::start;
                if (failed_ || rec.len == 0) return false;
                if (rec.len & 1) rec.nibbles[rec.len >> 1] &= 0xF0;
                return true;
            }
            const uint8_t *base = buf_;
            switch (st_) {
                case state::start: {
                    // Skip blank lines, then consume the record marker.
                    while (pos_ < end_ && base[pos_] <= ' ') ++pos_;
                    if (pos_ == end_) break;
                    const uint8_t c = base[pos_];
                    if (c != '>' && c != '@') {
                        failed_ = true;
                        eof_ = true;
                        return false;
                    }
                    ++pos_;
                    name_.clear();
                    name_done_ = false;
                    seq_lines_ = 0;
                    fastq_ = false;
                    qual_wrapped_ = false;
                    truncated_ = false;
                    st_ = state::header;
                    break;
                }
                case state::header: {
                    const uint8_t *nl = static_cast<const uint8_t *>(
                        memchr(base + pos_, '\n', end_ - pos_));
                    const size_t stop =
                        nl ? static_cast<size_t>(nl - base) : end_;
                    if (!name_done_) {
                        // Name is the header up to the first whitespace; the
                        // comment after it is discarded.
                        size_t w = pos_;
                        while (w < stop && base[w] > ' ') ++w;
                        name_.append(
                            reinterpret_cast<const char *>(base + pos_),
                            w - pos_);
                        if (w < stop) name_done_ = true;
                    }
                    if (!nl) {
                        pos_ = end_;
                        break;
                    }
                    pos_ = stop + 1;
                    st_ = state::seq;
                    line_start_ = true;
                    break;
                }
                case state::plus: {
                    const uint8_t *nl = static_cast<const uint8_t *>(
                        memchr(base + pos_, '\n', end_ - pos_));
                    if (!nl) {
                        pos_ = end_;
                        break;
                    }
                    pos_ = static_cast<size_t>(nl - base) + 1;
                    st_ = state::qual;
                    qual_got_ = 0;
                    break;
                }
                case state::seq: {
                    if (line_start_) {
                        const uint8_t c = base[pos_];
                        if (c == '>' || c == '@') {
                            // Next record's header; leave it unconsumed.
                            st_ = state::start;
                            if (rec.len & 1) rec.nibbles[rec.len >> 1] &= 0xF0;
                            return true;
                        }
                        if (c == '+') {
                            st_ = state::plus;
                            fastq_ = true;
                            break;
                        }
                        ++seq_lines_;
                    }
                    const uint8_t *nl = static_cast<const uint8_t *>(
                        memchr(base + pos_, '\n', end_ - pos_));
                    const size_t stop =
                        nl ? static_cast<size_t>(nl - base) : end_;
                    pack_bases(rec, base + pos_, base + stop);
                    pos_ = nl ? stop + 1 : end_;
                    line_start_ = (nl != nullptr);
                    break;
                }
                case state::qual: {
                    // Quality is never used; only its length matters (it must
                    // equal the sequence length, which is what makes a quality
                    // line starting with '@' unambiguous). Skip it in spans:
                    // one memchr finds the next newline, everything before it
                    // counts. Single-line quality is one memchr over len
                    // bytes; multi-line quality still works. '\r' never counts.
                    while (pos_ < end_ && qual_got_ < rec.len) {
                        const size_t need = rec.len - qual_got_;
                        const size_t span =
                            end_ - pos_ < need ? end_ - pos_ : need;
                        const uint8_t *p = base + pos_;
                        const uint8_t *nl =
                            static_cast<const uint8_t *>(memchr(p, '\n', span));
                        if (!nl) {
                            // A '\r' ending the span is a CRLF split across
                            // blocks (or the need-cut landing on one).
                            qual_got_ += static_cast<uint32_t>(
                                p[span - 1] == '\r' ? span - 1 : span);
                            pos_ += span;
                            continue;
                        }
                        size_t n = static_cast<size_t>(nl - p);
                        if (n > 0 && nl[-1] == '\r') --n;
                        qual_got_ += static_cast<uint32_t>(n);
                        pos_ = static_cast<size_t>(nl - base) + 1;
                        if (qual_got_ < rec.len) qual_wrapped_ = true;
                    }
                    if (qual_got_ == rec.len) {
                        st_ = state::start;
                        if (rec.len & 1) rec.nibbles[rec.len >> 1] &= 0xF0;
                        return true;
                    }
                    break;
                }
            }
        }
    }
};

#endif
