#ifndef _LASM_BAM_READER_H_
#define _LASM_BAM_READER_H_

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "fastx.h"
#include "input_stream.h"

// BAM input. Only the sequence of each record is kept; flags are not looked
// at, so secondary and supplementary records are read like any other. Records
// come out as seq_record in the packed nt16 layout of the BAM SEQ field.

// Contiguous byte ranges over an input_stream's blocks. A range inside the
// current block points into it; only a range straddling blocks is copied.
class byte_cursor {
    input_stream &in_;
    const uint8_t *p_ = nullptr;
    size_t left_ = 0;
    std::vector<uint8_t> carry_;

    bool refill() { return in_.next(p_, left_); }

   public:
    explicit byte_cursor(input_stream &in) : in_(in) {}

    // The next `n` (> 0) bytes, valid until the next call. nullptr when the
    // stream ends first; `got` is then the count that was left.
    const uint8_t *take(size_t n, size_t &got) {
        if (left_ >= n) {
            const uint8_t *r = p_;
            p_ += n;
            left_ -= n;
            got = n;
            return r;
        }
        carry_.resize(n);
        got = left_;
        if (left_ > 0) memcpy(carry_.data(), p_, left_);
        left_ = 0;
        while (got < n) {
            if (!refill()) return nullptr;
            const size_t c = std::min(n - got, left_);
            memcpy(carry_.data() + got, p_, c);
            p_ += c;
            left_ -= c;
            got += c;
        }
        return carry_.data();
    }

    // Skips `n` bytes; false if the stream ends first.
    bool skip(size_t n) {
        while (n > left_) {
            n -= left_;
            left_ = 0;
            if (!refill()) return false;
        }
        p_ += n;
        left_ -= n;
        return true;
    }
};

class bam_reader {
    // Fixed part of a BAM record after block_size (SAM spec 4.2).
    static constexpr size_t kBamCoreSize = 32;
    static constexpr size_t kBamReadNameLenAt = 8;
    static constexpr size_t kBamCigarCountAt = 12;
    static constexpr size_t kBamSeqLenAt = 16;
    // block_size is an int32 in the spec; htslib refuses larger records.
    static constexpr uint32_t kBamMaxBlockSize = INT32_MAX - kBamCoreSize;

    input_stream &in_;
    byte_cursor cur_{in_};
    std::string error_;

    bool fail(const char *what) {
        error_ = in_.error();
        if (error_.empty()) error_ = what;
        return false;
    }

   public:
    // Reads the BAM file `in` (opened by the caller and not owned) from its
    // start.
    explicit bam_reader(input_stream &in) : in_(in) {}
    bam_reader(const bam_reader &) = delete;
    bam_reader &operator=(const bam_reader &) = delete;

    // Reads past the header; call once before next().
    bool read_header() {
        size_t got = 0;
        const uint8_t *magic = cur_.take(4, got);
        if (!magic || memcmp(magic, "BAM\1", 4) != 0) {
            return fail("not a BAM file");
        }
        const uint8_t *l_text = cur_.take(4, got);
        if (!l_text || !cur_.skip(load_le32(l_text))) {
            return fail("truncated BAM header");
        }
        const uint8_t *n_ref_p = cur_.take(4, got);
        if (!n_ref_p) return fail("truncated BAM header");
        const uint32_t n_ref = load_le32(n_ref_p);
        for (uint32_t i = 0; i < n_ref; ++i) {
            const uint8_t *l_name = cur_.take(4, got);
            // Name, then its l_ref.
            if (!l_name || !cur_.skip(load_le32(l_name) + size_t{4})) {
                return fail("truncated BAM header");
            }
        }
        return true;
    }

    // Fills `rec` with the next record's sequence. False at the end of the
    // input or on error.
    bool next(seq_record &rec) {
        if (!error_.empty()) return false;
        size_t got = 0;
        const uint8_t *size_p = cur_.take(4, got);
        if (!size_p) {
            if (got == 0 && in_.error().empty()) return false;
            return fail("truncated BAM record");
        }
        const uint32_t block_size = load_le32(size_p);
        if (block_size < kBamCoreSize || block_size > kBamMaxBlockSize) {
            return fail("corrupt BAM record: bad block_size");
        }
        const uint8_t *r = cur_.take(block_size, got);
        if (!r) return fail("truncated BAM record");
        const size_t l_read_name = r[kBamReadNameLenAt];
        const size_t n_cigar_op = load_le16(r + kBamCigarCountAt);
        const int32_t l_seq = static_cast<int32_t>(load_le32(r + kBamSeqLenAt));
        const size_t seq_at = kBamCoreSize + l_read_name + 4 * n_cigar_op;
        if (l_seq < 0) return fail("corrupt BAM record: negative l_seq");
        const size_t bases = static_cast<size_t>(l_seq);
        const size_t seq_bytes = (bases + 1) / 2;
        // Sequence, then a quality byte per base.
        if (seq_at + seq_bytes + bases > block_size) {
            return fail("corrupt BAM record: sequence past the record end");
        }
        rec.reserve_bases(static_cast<uint32_t>(bases));
        memcpy(rec.nibbles, r + seq_at, seq_bytes);
        rec.len = static_cast<uint32_t>(bases);
        return true;
    }

    bool failed() const { return !error_.empty(); }
    const std::string &error() const { return error_; }
};

#endif
