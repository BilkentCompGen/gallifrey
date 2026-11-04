#ifndef _LASM_IO_FASTX_CHUNKS_H_
#define _LASM_IO_FASTX_CHUNKS_H_

#include <cstdint>
#include <string_view>
#include <vector>

#include "chunks.h"

// Record starts of FASTA/FASTQ, uncompressed or bgzipped (the same starts,
// looked for in the decompressed stream).
//
// FASTA: a '>' at a line start always begins a record (sequence lines cannot
// hold '>'). A cut goes only before one whose previous record has a base,
// because a chunk's reader drops an empty last record.
// FASTQ: a '@' line may be a quality line. A cut goes only before a '@' line
// that, with the next three lines, is a four-line record (sequence of bases,
// '+' line, quality as long as the sequence) and is followed by a second such
// record or the end of the file. In four-line FASTQ a quality line starting
// with '@' cannot pass: the line after it is the next header, which a sequence
// line may not start with. Wrapped FASTQ can fake both records, so the file's
// first record must be four-line FASTQ, and the readers of the chunks check
// that every record they read is (record_fits_chunk).

enum class fastx_kind { unknown, fasta, fastq };

// A line of a window: [begin, end) without its '\n' and a '\r' before it;
// `next` is where the following line starts.
struct text_line {
    size_t begin;
    size_t end;
    size_t next;
};

// The line starting at `at` in `w`. The last line of the file (`at_eof`: `w`
// ends with the file) needs no '\n'. False when `w` ends inside the line or
// holds no line at `at`.
inline bool line_at(std::string_view w, size_t at, bool at_eof,
                    text_line &line) {
    if (at >= w.size()) return false;
    const size_t nl = w.find('\n', at);
    if (nl == std::string_view::npos && !at_eof) return false;
    line.begin = at;
    if (nl == std::string_view::npos) {
        line.end = w.size();
        line.next = w.size();
    } else {
        line.end = nl;
        line.next = nl + 1;
    }
    if (line.end > at && w[line.end - 1] == '\r') --line.end;
    return true;
}

// Whether the four lines at `at` form one four-line FASTQ record as the
// reader parses it: a '@' header, a sequence of bytes above ' ' that does not
// start with '@', '>' or '+', a '+' line, and a quality line as long as the
// sequence. `next` gets the start of the line after the record.
inline verdict fastq_record_at(std::string_view w, size_t at, bool at_eof,
                               size_t &next) {
    text_line lines[4];
    size_t p = at;
    for (text_line &line : lines) {
        if (!line_at(w, p, at_eof, line)) {
            if (at_eof) return verdict::no;
            return verdict::unknown;
        }
        p = line.next;
    }
    if (w[at] != '@') return verdict::no;
    const std::string_view seq =
        w.substr(lines[1].begin, lines[1].end - lines[1].begin);
    if (!seq.empty() &&
        (seq.front() == '@' || seq.front() == '>' || seq.front() == '+')) {
        return verdict::no;
    }
    for (const char c : seq) {
        if (static_cast<unsigned char>(c) <= ' ') return verdict::no;
    }
    if (lines[2].end == lines[2].begin || w[lines[2].begin] != '+') {
        return verdict::no;
    }
    if (lines[3].end - lines[3].begin != seq.size()) return verdict::no;
    next = p;
    return verdict::yes;
}

// Whether a FASTQ record starts at `at`: a four-line record followed by the
// end of the file or by a second four-line record.
inline verdict fastq_start_at(std::string_view w, size_t at, bool at_eof) {
    size_t next = 0;
    const verdict first = fastq_record_at(w, at, at_eof, next);
    if (first != verdict::yes) return first;
    if (next >= w.size()) {
        if (at_eof) return verdict::yes;
        return verdict::unknown;
    }
    size_t after = 0;
    return fastq_record_at(w, next, at_eof, after);
}

// Whether the '>' line at `at` starts a record of a chunk that the record
// before it ends: that record must hold a base, since the reader drops an
// empty record when its chunk ends there (as at the end of the file) but
// returns it when the next header follows. So the line before must be a
// sequence line with a byte above ' ' and lie wholly in `w` (from byte 1).
inline verdict fasta_start_at(std::string_view w, size_t at) {
    if (at < 2) return verdict::no;
    const size_t nl = w.rfind('\n', at - 2);
    if (nl == std::string_view::npos) return verdict::no;
    const std::string_view prev = w.substr(nl + 1, at - 1 - (nl + 1));
    if (prev.empty() || prev.front() == '>' || prev.front() == '@') {
        return verdict::no;
    }
    for (const char c : prev) {
        if (static_cast<unsigned char>(c) > ' ') return verdict::yes;
    }
    return verdict::no;
}

// The first record start at a line start after byte 0 of `w` (which may lie
// inside a line); `at_eof`: `w` ends with the file. When `w` ends before the
// answer, `found` gets where the search stopped.
inline verdict find_record_start(std::string_view w, fastx_kind kind,
                                 bool at_eof, size_t &found) {
    char marker = '@';
    if (kind == fastx_kind::fasta) marker = '>';
    size_t pos = 0;
    while (true) {
        const size_t nl = w.find('\n', pos);
        if (nl == std::string_view::npos || nl + 1 >= w.size()) {
            found = w.size();
            if (at_eof) return verdict::no;
            return verdict::unknown;
        }
        const size_t at = nl + 1;
        if (w[at] == marker) {
            verdict v = verdict::yes;
            if (kind == fastx_kind::fastq) v = fastq_start_at(w, at, at_eof);
            if (kind == fastx_kind::fasta) v = fasta_start_at(w, at);
            if (v == verdict::yes) {
                found = at;
                return verdict::yes;
            }
            if (v == verdict::unknown) {
                found = at;
                return verdict::unknown;
            }
        }
        pos = at;
    }
}

// The kind of the first record of a stream whose first bytes are `w`
// (`at_eof`: all of it): `yes` with `kind` set when it is '>' or a four-line
// FASTQ record followed by another one (or the end), `no` when it is
// neither, `unknown` when `w` ends before the answer.
inline verdict first_record_kind_in(std::string_view w, bool at_eof,
                                    fastx_kind &kind) {
    size_t at = 0;
    while (at < w.size() && static_cast<unsigned char>(w[at]) <= ' ') ++at;
    if (at == w.size()) {
        if (at_eof) return verdict::no;
        return verdict::unknown;
    }
    if (w[at] == '>') {
        kind = fastx_kind::fasta;
        return verdict::yes;
    }
    if (w[at] != '@') return verdict::no;
    const verdict v = fastq_start_at(w, at, at_eof);
    if (v == verdict::yes) kind = fastx_kind::fastq;
    return v;
}

// The first record start at a line start after byte 0 of the window (which
// may lie inside a line), growing the window while the search is before the
// part's end `limit`.
inline bool first_fastx_start(chunk_window &w, fastx_kind kind, uint64_t limit,
                              size_t &found) {
    while (true) {
        const verdict v = find_record_start(w.bytes(), kind, w.at_end(), found);
        if (v != verdict::unknown) return v == verdict::yes;
        if (!w.can_extend() || w.position(found).at >= limit) return false;
        w.extend();
    }
}

// Cuts the FASTA/FASTQ file `fd` (`size` bytes; `bgzf`: bgzipped) into up to
// `parts` chunks and sets `kind`. Empty when it is not cut: no BGZF
// end-of-file marker, a first record the checks do not accept (wrapped
// FASTQ), or no second record start.
inline std::vector<chunk> cut_fastx_file(int fd, uint64_t size, bool bgzf,
                                         size_t parts, fastx_kind &kind) {
    kind = fastx_kind::unknown;
    if (parts < 2 || size == 0 || (bgzf && !has_bgzf_eof(fd, size))) return {};
    chunk_window w(fd, size, bgzf, 0);
    verdict v = verdict::unknown;
    while (v == verdict::unknown && w.extend()) {
        v = first_record_kind_in(w.bytes(), w.at_end(), kind);
    }
    if (v != verdict::yes || w.failed()) {
        kind = fastx_kind::unknown;
        return {};
    }
    const fastx_kind k = kind;
    return cut_at_record_starts(
        fd, size, bgzf, parts, {0, 0},
        [k](chunk_window &cw, uint64_t limit, size_t &found) {
            return first_fastx_start(cw, k, limit, found);
        });
}

// Whether a record read from a chunk of a cut file has the shape the cut
// relied on: four-line FASTQ in a FASTQ file, no FASTQ in a FASTA file. Any
// other record means a cut may sit inside a record, so the run must stop.
template <class Reader>
bool record_fits_chunk(fastx_kind kind, const Reader &reader) {
    if (kind == fastx_kind::fastq) return reader.four_line_fastq_record();
    return !reader.fastq_record();
}

#endif
