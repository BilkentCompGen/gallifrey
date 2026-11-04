#ifndef _LASM_IO_CHUNKS_H_
#define _LASM_IO_CHUNKS_H_

#include <libdeflate.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "input_stream.h"

// Reading one input file in chunks: the file is cut into chunks that each
// begin at a record start, and every chunk is read on its own. A chunk begins
// only where the sequential reader would begin a record and ends where the
// next chunk begins, so the chunks together hold every record exactly once.
// Where the record starts are depends on the format (fastx_chunks.h,
// bam_chunks.h); this file holds what all formats share: chunk positions, the
// source that reads a chunk (uncompressed or BGZF), and the cut loop.
//
// In a BGZF file a position is the file offset of a block and an offset into
// its data; near each part start the cut goes at the first record start after
// the first block start, where a block start is a BGZF header whose block
// size leads to two more headers (or the end of the file). In an
// uncompressed file a position is a file offset (and offset 0).

// The answer of a check over a window of the file; `unknown` when the window
// ends before the answer.
enum class verdict { yes, no, unknown };

// Reads up to `n` bytes at `offset`, retrying short reads. Returns the bytes
// read (fewer than `n` only at end of file) or -1 on error.
inline int64_t pread_full(int fd, uint8_t *buf, size_t n, uint64_t offset) {
    size_t got = 0;
    while (got < n) {
        const ssize_t r =
            ::pread(fd, buf + got, n - got, static_cast<off_t>(offset + got));
        if (r < 0 && errno == EINTR) continue;
        if (r < 0) return -1;
        if (r == 0) break;
        got += static_cast<size_t>(r);
    }
    return static_cast<int64_t>(got);
}

// A position in a file's (decompressed) stream.
struct chunk_position {
    uint64_t at;      // file offset (of the BGZF block)
    uint32_t offset;  // offset into the BGZF block's data; 0 uncompressed
};

inline bool operator<(const chunk_position &a, const chunk_position &b) {
    if (a.at != b.at) return a.at < b.at;
    return a.offset < b.offset;
}

// Stream bytes [begin, end); the last chunk of a file ends at {size, 0}.
struct chunk {
    chunk_position begin;
    chunk_position end;
};

// The stream bytes of one chunk of an open file, read with pread(2) and, for
// BGZF, inflated and CRC-checked on the calling thread: the Source of a
// fastx_reader or bam_reader. The file must not change while it is read.
class chunk_source {
    // Bytes read at once: compressed BGZF blocks in windows of several (fewer
    // reads), uncompressed text in pieces the reader parses as it goes. A
    // source holds one such buffer, plus a 64 KiB block and an inflater for
    // BGZF: about 1 MiB per reader of uncompressed chunks, 4.1 MiB of BGZF.
    static constexpr size_t kBgzfReadSize = size_t{4} << 20;
    static constexpr size_t kPlainReadSize = size_t{1} << 20;
    int fd_;
    uint64_t file_size_;
    bool bgzf_;
    uint64_t at_;    // file offset of the next read (or BGZF block)
    uint32_t skip_;  // bytes of the first block before the chunk
    chunk_position end_;
    bool done_ = false;
    uint64_t last_at_ = 0;
    size_t read_size_;
    std::unique_ptr<uint8_t[]> buf_;
    uint64_t buf_at_ = 0;
    size_t buf_size_ = 0;
    std::unique_ptr<uint8_t[]> data_;
    std::unique_ptr<libdeflate_decompressor,
                    decltype(&libdeflate_free_decompressor)>
        inflater_{nullptr, &libdeflate_free_decompressor};
    std::string error_;

    bool fail(const char *what) {
        error_ = what;
        return false;
    }

    // File bytes [at, at + n) (n at most a block), or nullptr past the end
    // of the file or on a read error.
    const uint8_t *compressed(uint64_t at, size_t n) {
        if (at >= buf_at_ && at + n <= buf_at_ + buf_size_) {
            return buf_.get() + (at - buf_at_);
        }
        if (at + n > file_size_) return nullptr;
        const size_t want = static_cast<size_t>(
            std::min<uint64_t>(read_size_, file_size_ - at));
        const int64_t got = pread_full(fd_, buf_.get(), want, at);
        if (got < 0) {
            error_ = std::string("read error: ") + strerror(errno);
            return nullptr;
        }
        if (got < static_cast<int64_t>(n)) return nullptr;
        buf_at_ = at;
        buf_size_ = static_cast<size_t>(got);
        return buf_.get();
    }

    // Inflates the block at at_ into data_; its data and block sizes go to
    // `size` and `comp_size`.
    bool inflate_block(uint32_t &size, uint32_t &comp_size) {
        const uint8_t *h = compressed(at_, kBgzfHeaderSize);
        if (!h || !is_bgzf_header(h, kBgzfHeaderSize)) {
            if (!error_.empty()) return false;
            return fail("not a BGZF block");
        }
        comp_size = load_le16(h + 16) + uint32_t{1};
        if (comp_size < kBgzfHeaderSize + kBgzfFooterSize) {
            return fail("corrupt BGZF block size");
        }
        const uint8_t *c = compressed(at_, comp_size);
        if (!c) {
            if (!error_.empty()) return false;
            return fail("truncated BGZF block");
        }
        size = load_le32(c + comp_size - 4);
        if (size > kBgzfMaxBlockSize) {
            return fail("corrupt BGZF block: uncompressed size over 64 KiB");
        }
        if (libdeflate_deflate_decompress(
                inflater_.get(), c + kBgzfHeaderSize,
                comp_size - kBgzfHeaderSize - kBgzfFooterSize, data_.get(),
                size, nullptr) != LIBDEFLATE_SUCCESS) {
            return fail("corrupt BGZF block: inflate failed");
        }
        const uint32_t crc = load_le32(c + comp_size - kBgzfFooterSize);
        if (libdeflate_crc32(0, data_.get(), size) != crc) {
            return fail("corrupt BGZF block: CRC32 mismatch");
        }
        return true;
    }

    bool next_uncompressed(const uint8_t *&p, size_t &n) {
        if (at_ >= end_.at) return false;
        const size_t want =
            static_cast<size_t>(std::min<uint64_t>(read_size_, end_.at - at_));
        const int64_t got = pread_full(fd_, buf_.get(), want, at_);
        if (got < 0) {
            error_ = std::string("read error: ") + strerror(errno);
            return false;
        }
        if (got == 0) return fail("the file became shorter while it was read");
        last_at_ = at_;
        at_ += static_cast<uint64_t>(got);
        p = buf_.get();
        n = static_cast<size_t>(got);
        return true;
    }

    bool next_bgzf(const uint8_t *&p, size_t &n) {
        while (!done_) {
            if (at_ == end_.at && end_.offset == 0 && skip_ == 0) break;
            if (at_ >= file_size_) return fail("chunk ends past the file");
            uint32_t size = 0;
            uint32_t comp_size = 0;
            if (!inflate_block(size, comp_size)) return false;
            last_at_ = at_;
            const uint32_t begin = skip_;
            skip_ = 0;
            uint32_t end = size;
            if (at_ == end_.at) {
                end = end_.offset;
                done_ = true;
            }
            at_ += comp_size;
            if (begin > end || end > size) {
                return fail("chunk outside its BGZF block");
            }
            if (end > begin) {
                p = data_.get() + begin;
                n = end - begin;
                return true;
            }
        }
        done_ = true;
        return false;
    }

   public:
    chunk_source(int fd, uint64_t file_size, bool bgzf, chunk c)
        : fd_(fd),
          file_size_(file_size),
          bgzf_(bgzf),
          at_(c.begin.at),
          skip_(c.begin.offset),
          end_(c.end) {
        read_size_ = kPlainReadSize;
        if (bgzf_) {
            read_size_ = kBgzfReadSize;
            data_ = std::make_unique<uint8_t[]>(kBgzfMaxBlockSize);
            inflater_.reset(libdeflate_alloc_decompressor());
            if (!inflater_) fail("out of memory");
        }
        buf_ = std::make_unique<uint8_t[]>(read_size_);
    }
    chunk_source(const chunk_source &) = delete;
    chunk_source &operator=(const chunk_source &) = delete;

    // The next bytes of the chunk; false at its end or on error.
    bool next(const uint8_t *&p, size_t &n) {
        if (!error_.empty()) return false;
        if (bgzf_) return next_bgzf(p, n);
        return next_uncompressed(p, n);
    }

    // File offset of the read (or BGZF block) the last next() gave bytes of.
    uint64_t last_at() const { return last_at_; }

    // Why the chunk could not be read; empty after a clean end.
    const std::string &error() const { return error_; }
};

// Most stream bytes looked at for one cut.
inline constexpr size_t kMaxCutScan = size_t{1} << 27;

// Stream bytes from file offset `from` (a BGZF block start in a BGZF file)
// on, gathered for the record-start checks, with the file offset each piece
// came from.
class chunk_window {
    chunk_source src_;
    bool bgzf_;
    std::string bytes_;
    struct piece {
        size_t at;  // offset in bytes_
        uint64_t from;
    };
    std::vector<piece> pieces_;
    bool at_end_ = false;

   public:
    chunk_window(int fd, uint64_t size, bool bgzf, uint64_t from)
        : src_(fd, size, bgzf, {{from, 0}, {size, 0}}), bgzf_(bgzf) {}

    // At least doubles bytes(), up to the end of the file or kMaxCutScan, so
    // that checks repeated over the growing window take linear time. False
    // when the window cannot change any more (see can_extend).
    bool extend() {
        if (!can_extend()) return false;
        const size_t want = std::max<size_t>(2 * bytes_.size(), 1);
        while (bytes_.size() < want && bytes_.size() < kMaxCutScan) {
            const uint8_t *p = nullptr;
            size_t n = 0;
            if (!src_.next(p, n)) {
                at_end_ = true;
                break;
            }
            pieces_.push_back({bytes_.size(), src_.last_at()});
            bytes_.append(reinterpret_cast<const char *>(p), n);
        }
        return true;
    }
    std::string_view bytes() const { return bytes_; }
    bool at_end() const { return at_end_; }
    bool failed() const { return !src_.error().empty(); }
    // Whether extend() may still change the window.
    bool can_extend() const {
        return !at_end_ && !failed() && bytes_.size() < kMaxCutScan;
    }

    // The stream position of offset `i` of bytes().
    chunk_position position(size_t i) const {
        size_t k = pieces_.size() - 1;
        while (pieces_[k].at > i) --k;
        const size_t in = i - pieces_[k].at;
        if (bgzf_) return {pieces_[k].from, static_cast<uint32_t>(in)};
        return {pieces_[k].from + in, 0};
    }
};

// Whether a BGZF block starts at file offset `at`: its header and size lead
// to two more headers, or to the end of the file.
inline bool bgzf_block_at(int fd, uint64_t size, uint64_t at) {
    uint64_t pos = at;
    for (int i = 0; i < 3; ++i) {
        if (pos == size && i > 0) return true;
        uint8_t h[kBgzfHeaderSize];
        if (pos + kBgzfHeaderSize > size ||
            pread_full(fd, h, kBgzfHeaderSize, pos) !=
                static_cast<int64_t>(kBgzfHeaderSize) ||
            !is_bgzf_header(h, kBgzfHeaderSize)) {
            return false;
        }
        pos += load_le16(h + 16) + uint64_t{1};
    }
    return true;
}

// The first BGZF block start at or after `target` (within one block size).
inline bool first_block_start(int fd, uint64_t size, uint64_t target,
                              uint64_t &found) {
    const size_t span = static_cast<size_t>(
        std::min<uint64_t>(kBgzfMaxBlockSize + kBgzfHeaderSize, size - target));
    std::vector<uint8_t> w(span);
    const int64_t got = pread_full(fd, w.data(), span, target);
    if (got < 0) return false;
    for (size_t i = 0; i + kBgzfHeaderSize <= static_cast<size_t>(got); ++i) {
        if (w[i] == 0x1f && is_bgzf_header(w.data() + i, got - i) &&
            bgzf_block_at(fd, size, target + i)) {
            found = target + i;
            return true;
        }
    }
    return false;
}

// Whether the file ends with the BGZF end-of-file marker.
inline bool has_bgzf_eof(int fd, uint64_t size) {
    uint8_t tail[sizeof(kBgzfEofBlock)];
    return size >= sizeof(kBgzfEofBlock) &&
           pread_full(fd, tail, sizeof(tail), size - sizeof(tail)) ==
               static_cast<int64_t>(sizeof(tail)) &&
           memcmp(tail, kBgzfEofBlock, sizeof(tail)) == 0;
}

// Where part k of `parts` equal parts of `size` bytes begins.
inline uint64_t part_start(uint64_t size, size_t parts, size_t k) {
    return size / parts * k + size % parts * k / parts;
}

// Cuts the file `fd` (`size` bytes) into up to `parts` chunks: the first
// begins at `first`, each other one at the record start that
// `find_start(window, limit, found)` finds in a window opened near a part
// start (at the first BGZF block start after it, or one byte before it
// uncompressed, so that a record there starts a line of the window). The
// search may give up once it passes the next part start `limit`, where the
// next part's search begins. Empty when fewer than two chunks are found.
template <class FindStart>
std::vector<chunk> cut_at_record_starts(int fd, uint64_t size, bool bgzf,
                                        size_t parts, chunk_position first,
                                        FindStart &&find_start) {
    std::vector<chunk_position> starts{first};
    for (size_t k = 1; k < parts; ++k) {
        const uint64_t target = part_start(size, parts, k);
        const uint64_t limit = part_start(size, parts, k + 1);
        if (target <= starts.back().at) continue;
        uint64_t from = target - 1;
        if (bgzf && !first_block_start(fd, size, target, from)) continue;
        chunk_window w(fd, size, bgzf, from);
        size_t found = 0;
        if (!w.extend() || !find_start(w, limit, found)) continue;
        const chunk_position p = w.position(found);
        if (starts.back() < p) starts.push_back(p);
    }
    std::vector<chunk> chunks;
    if (starts.size() < 2) return chunks;
    for (size_t i = 0; i < starts.size(); ++i) {
        chunk_position end{size, 0};
        if (i + 1 < starts.size()) end = starts[i + 1];
        chunks.push_back({starts[i], end});
    }
    return chunks;
}

#endif
