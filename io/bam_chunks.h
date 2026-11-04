#ifndef _LASM_IO_BAM_CHUNKS_H_
#define _LASM_IO_BAM_CHUNKS_H_

#include <cstdint>
#include <cstring>
#include <string_view>
#include <vector>

#include "bam_reader.h"
#include "chunks.h"

// Record starts of unaligned BAM. A record start is an unaligned BAM record
// (refID, pos, mate refID and mate pos -1, no CIGAR, unmapped flag, printable
// NUL-terminated name, quality of phred values or missing, auxiliary fields
// that end at the record's end) followed by two more such records (or the end
// of the stream). The first record after the header must pass, else the file
// is not cut (aligned BAM, unusual input).
//
// A chunk's reader stops at the next chunk's first record. BAM records carry
// their length, so it ends cleanly only if a record ends exactly there: a cut
// inside a record fails the run ("truncated BAM record") instead of changing
// it.

// Records longer than this are never taken as cut starts (they only leave
// a chunk longer).
inline constexpr uint32_t kMaxCutRecord = uint32_t{1} << 26;

// Whether the auxiliary fields in [p, end) parse and end exactly at `end`.
inline bool aux_fields_fit(const uint8_t *p, const uint8_t *end) {
    while (p < end) {
        if (end - p < 3) return false;
        const bool tag_ok = ((p[0] | 0x20) >= 'a' && (p[0] | 0x20) <= 'z') &&
                            (((p[1] | 0x20) >= 'a' && (p[1] | 0x20) <= 'z') ||
                             (p[1] >= '0' && p[1] <= '9'));
        if (!tag_ok) return false;
        const uint8_t type = p[2];
        p += 3;
        size_t size = 0;
        switch (type) {
            case 'A':
            case 'c':
            case 'C':
                size = 1;
                break;
            case 's':
            case 'S':
                size = 2;
                break;
            case 'i':
            case 'I':
            case 'f':
                size = 4;
                break;
            case 'Z':
            case 'H': {
                const void *nul = memchr(p, 0, static_cast<size_t>(end - p));
                if (!nul) return false;
                p = static_cast<const uint8_t *>(nul) + 1;
                continue;
            }
            case 'B': {
                if (end - p < 5) return false;
                const uint8_t sub = p[0];
                size_t elem = 0;
                if (sub == 'c' || sub == 'C') elem = 1;
                if (sub == 's' || sub == 'S') elem = 2;
                if (sub == 'i' || sub == 'I' || sub == 'f') elem = 4;
                if (elem == 0) return false;
                const uint64_t count = load_le32(p + 1);
                p += 5;
                if (count * elem > static_cast<uint64_t>(end - p)) return false;
                p += count * elem;
                continue;
            }
            default:
                return false;
        }
        if (static_cast<size_t>(end - p) < size) return false;
        p += size;
    }
    return p == end;
}

// Whether an unaligned BAM record starts at `at` of `w` (see the top of the
// file); `next` gets the offset after it. `at_end`: `w` ends with the stream.
inline verdict bam_record_at(std::string_view w, size_t at, bool at_end,
                             size_t &next) {
    const auto *d = reinterpret_cast<const uint8_t *>(w.data());
    const size_t n = w.size();
    auto more = [&] {
        if (at_end) return verdict::no;
        return verdict::unknown;
    };
    if (n - at < 4 + kBamCoreSize) return more();
    const uint8_t *r = d + at + 4;
    const uint32_t block_size = load_le32(d + at);
    if (load_le32(r) != 0xFFFFFFFFu || load_le32(r + 4) != 0xFFFFFFFFu ||
        load_le32(r + 20) != 0xFFFFFFFFu || load_le32(r + 24) != 0xFFFFFFFFu ||
        load_le32(r + 28) != 0) {
        return verdict::no;
    }
    const size_t l_read_name = r[8];
    const uint32_t n_cigar_op = load_le16(r + 12);
    const uint32_t flag = load_le16(r + 14);
    const uint32_t l_seq = load_le32(r + 16);
    if (l_read_name < 1 || n_cigar_op != 0 || (flag & 4) == 0 ||
        (flag & ~0xFFFu) != 0 || l_seq > kMaxCutRecord ||
        block_size < kBamCoreSize || block_size > kMaxCutRecord) {
        return verdict::no;
    }
    const uint64_t fixed = kBamCoreSize + l_read_name +
                           (uint64_t{l_seq} + 1) / 2 + uint64_t{l_seq};
    if (fixed > block_size) return verdict::no;
    if (n - at - 4 < block_size) return more();
    const uint8_t *name = r + kBamCoreSize;
    if (name[l_read_name - 1] != 0) return verdict::no;
    for (size_t i = 0; i + 1 < l_read_name; ++i) {
        if (name[i] < '!' || name[i] > '~') return verdict::no;
    }
    const uint8_t *qual = name + l_read_name + (l_seq + 1) / 2;
    if (l_seq > 0 && qual[0] != 0xFF) {
        for (uint32_t i = 0; i < l_seq; ++i) {
            if (qual[i] > 93) return verdict::no;
        }
    }
    if (!aux_fields_fit(qual + l_seq, r + block_size)) return verdict::no;
    next = at + 4 + block_size;
    return verdict::yes;
}

// Whether a record start is at `at`: three records in a row, or fewer that
// end with the stream.
inline verdict bam_start_at(std::string_view w, size_t at, bool at_end) {
    size_t pos = at;
    for (int i = 0; i < 3; ++i) {
        if (pos == w.size() && at_end && i > 0) return verdict::yes;
        size_t next = 0;
        const verdict v = bam_record_at(w, pos, at_end, next);
        if (v != verdict::yes) return v;
        pos = next;
    }
    return verdict::yes;
}

// The first record start in the window, growing it while the search is
// before the part's end `limit`.
inline bool first_bam_start(chunk_window &w, uint64_t limit, size_t &found) {
    size_t at = 0;
    while (true) {
        const std::string_view b = w.bytes();
        const auto *d = reinterpret_cast<const uint8_t *>(b.data());
        for (; at + 12 <= b.size(); ++at) {
            if (load_le32(d + at + 4) != 0xFFFFFFFFu ||
                load_le32(d + at + 8) != 0xFFFFFFFFu) {
                continue;
            }
            const verdict v = bam_start_at(b, at, w.at_end());
            if (v == verdict::yes) {
                found = at;
                return true;
            }
            if (v == verdict::unknown) break;
        }
        if (!w.can_extend() || w.position(at).at >= limit) return false;
        w.extend();
    }
}

// Decompressed size of the BAM header at the start of `b`; 0 when `b` does
// not hold all of it.
inline size_t bam_header_size(std::string_view b) {
    const auto *d = reinterpret_cast<const uint8_t *>(b.data());
    if (b.size() < 12 || memcmp(d, "BAM\1", 4) != 0) return 0;
    uint64_t at = 8 + uint64_t{load_le32(d + 4)};
    if (at + 4 > b.size()) return 0;
    const uint32_t n_ref = load_le32(d + at);
    at += 4;
    for (uint32_t i = 0; i < n_ref; ++i) {
        if (at + 4 > b.size()) return 0;
        at += 4 + uint64_t{load_le32(d + at)} + 4;
    }
    if (at > b.size()) return 0;
    return static_cast<size_t>(at);
}

// Cuts the records of the BGZF-compressed BAM file `fd` (`size` bytes) into
// up to `parts` chunks. Empty when it is not cut: no end-of-file marker, a
// header or first record the checks do not accept, or no second record start.
inline std::vector<chunk> cut_bam_file(int fd, uint64_t size, size_t parts) {
    if (parts < 2 || !has_bgzf_eof(fd, size)) return {};
    chunk_window w(fd, size, true, 0);
    size_t header = 0;
    while ((header = bam_header_size(w.bytes())) == 0) {
        if (!w.extend()) return {};
    }
    verdict v = bam_start_at(w.bytes(), header, w.at_end());
    while (v == verdict::unknown && w.extend()) {
        v = bam_start_at(w.bytes(), header, w.at_end());
    }
    if (v != verdict::yes || w.failed()) return {};
    return cut_at_record_starts(fd, size, true, parts, w.position(header),
                                first_bam_start);
}

#endif
