#ifndef _LASM_INPUT_STREAM_H_
#define _LASM_INPUT_STREAM_H_

#include <blockingconcurrentqueue.h>
#include <fcntl.h>
#include <libdeflate.h>
#include <unistd.h>
#include <zlib.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <expected>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

// Decompressed byte stream of one input file: uncompressed, gzip, or BGZF.
// The compression is detected from the first bytes, not the file name.
//
// BGZF (SAM spec 4.1) is a series of gzip members of at most 64 KiB each, with
// the compressed size in a header subfield and the uncompressed size in the
// footer. Its blocks are decompressed in parallel: one thread reads blocks, the
// io threads inflate them with libdeflate, and the consumer takes them back in
// file order. Plain gzip has no block structure and streams through zlib on
// the calling thread.

inline uint32_t load_le16(const uint8_t *p) {
    return static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8);
}

inline uint32_t load_le32(const uint8_t *p) {
    return load_le16(p) | (load_le16(p + 2) << 16);
}

// Reads up to `n` bytes, retrying short reads. Returns the bytes read (fewer
// than `n` only at end of file) or -1 on error.
inline int64_t read_full(int fd, uint8_t *buf, size_t n) {
    size_t got = 0;
    while (got < n) {
        const ssize_t r = ::read(fd, buf + got, n - got);
        if (r < 0 && errno == EINTR) continue;
        if (r < 0) return -1;
        if (r == 0) break;
        got += static_cast<size_t>(r);
    }
    return static_cast<int64_t>(got);
}

// Sequential bytes of a file, one read(2) per MiB.
class file_chunks {
    static constexpr size_t kChunkSize = 1u << 20;
    int fd_ = -1;
    std::unique_ptr<uint8_t[]> buf_ = std::make_unique<uint8_t[]>(kChunkSize);
    size_t pos_ = 0;
    size_t end_ = 0;

   public:
    void open(int fd) { fd_ = fd; }

    // Same contract as read_full.
    int64_t read(uint8_t *dst, size_t n) {
        size_t got = 0;
        while (got < n) {
            if (pos_ == end_) {
                const int64_t r = read_full(fd_, buf_.get(), kChunkSize);
                if (r < 0) return -1;
                if (r == 0) break;
                pos_ = 0;
                end_ = static_cast<size_t>(r);
            }
            const size_t c = std::min(n - got, end_ - pos_);
            memcpy(dst + got, buf_.get() + pos_, c);
            pos_ += c;
            got += c;
        }
        return static_cast<int64_t>(got);
    }
};

inline constexpr size_t kBgzfHeaderSize = 18;
inline constexpr size_t kBgzfFooterSize = 8;
inline constexpr size_t kBgzfMaxBlockSize = 65536;
// The empty block every BGZF writer appends; its absence means truncation.
inline constexpr uint8_t kBgzfEofBlock[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

inline bool is_gzip_magic(const uint8_t *h, size_t n) {
    return n >= 2 && h[0] == 0x1f && h[1] == 0x8b;
}

// gzip member with FEXTRA holding exactly the 6-byte BC subfield, as htslib
// requires of a BGZF block.
inline bool is_bgzf_header(const uint8_t *h, size_t n) {
    return n >= kBgzfHeaderSize && is_gzip_magic(h, n) && h[2] == 8 &&
           (h[3] & 4) != 0 && load_le16(h + 10) == 6 && h[12] == 'B' &&
           h[13] == 'C' && load_le16(h + 14) == 2;
}

// One BGZF block: the compressed member as read, then its contents.
struct bgzf_block {
    uint8_t comp[kBgzfMaxBlockSize];
    uint8_t data[kBgzfMaxBlockSize];
    uint32_t comp_size = 0;  // whole member, header and footer included
    uint32_t data_size = 0;
    bool end = false;             // no block: the file ended cleanly
    const char *error = nullptr;  // no block: reading or inflating failed
    bool eof_marker = false;      // this block is the 28-byte EOF marker
};

// Reads the next block. At end of file, marks the block `end` when the
// previous block was the EOF marker, else reports truncation.
inline void bgzf_read_block(file_chunks &src, bgzf_block &b,
                            bool after_eof_marker) {
    b.end = false;
    b.error = nullptr;
    b.eof_marker = false;
    const int64_t got = src.read(b.comp, kBgzfHeaderSize);
    if (got < 0) {
        b.error = "read error";
        return;
    }
    if (got == 0) {
        if (after_eof_marker) {
            b.end = true;
        } else {
            b.error = "truncated file: no BGZF end-of-file marker";
        }
        return;
    }
    if (!is_bgzf_header(b.comp, static_cast<size_t>(got))) {
        if (static_cast<size_t>(got) < kBgzfHeaderSize) {
            b.error = "truncated BGZF block";
        } else {
            b.error = "not a BGZF block";
        }
        return;
    }
    const uint32_t size = load_le16(b.comp + 16) + 1;
    if (size < kBgzfHeaderSize + kBgzfFooterSize) {
        b.error = "corrupt BGZF block size";
        return;
    }
    const size_t rest = size - kBgzfHeaderSize;
    if (src.read(b.comp + kBgzfHeaderSize, rest) !=
        static_cast<int64_t>(rest)) {
        b.error = "truncated BGZF block";
        return;
    }
    b.comp_size = size;
    b.data_size = load_le32(b.comp + size - 4);
    if (b.data_size > kBgzfMaxBlockSize) {
        b.error = "corrupt BGZF block: uncompressed size over 64 KiB";
        return;
    }
    b.eof_marker = size == sizeof(kBgzfEofBlock) &&
                   memcmp(b.comp, kBgzfEofBlock, sizeof(kBgzfEofBlock)) == 0;
}

// Inflates a block read by bgzf_read_block to its stated size and checks the
// CRC32 of the result.
inline void bgzf_inflate_block(libdeflate_decompressor *d, bgzf_block &b) {
    const uint8_t *payload = b.comp + kBgzfHeaderSize;
    const size_t payload_size = b.comp_size - kBgzfHeaderSize - kBgzfFooterSize;
    if (libdeflate_deflate_decompress(d, payload, payload_size, b.data,
                                      b.data_size,
                                      nullptr) != LIBDEFLATE_SUCCESS) {
        b.error = "corrupt BGZF block: inflate failed";
        return;
    }
    const uint32_t crc = load_le32(b.comp + b.comp_size - kBgzfFooterSize);
    if (libdeflate_crc32(0, b.data, b.data_size) != crc) {
        b.error = "corrupt BGZF block: CRC32 mismatch";
    }
}

// BGZF blocks in file order. With `threads` > 0, a reader thread and
// `threads` inflate threads run ahead through a ring of blocks; with 0, blocks
// are read and inflated on the calling thread.
class bgzf_stream {
    // Blocks in flight per inflate thread: enough to keep every thread busy
    // while the consumer works through a block.
    static constexpr size_t kBlocksPerThread = 8;
    static constexpr uint32_t kStopWorker = UINT32_MAX;

    struct slot {
        bgzf_block block;
        moodycamel::LightweightSemaphore ready{0, /*maxSpins=*/0};
    };

    file_chunks src_;
    std::string *error_ = nullptr;
    bool done_ = false;

    libdeflate_decompressor *inline_decompressor_ = nullptr;
    std::unique_ptr<bgzf_block> inline_block_;
    bool after_eof_marker_ = false;

    std::unique_ptr<slot[]> ring_;
    size_t ring_size_ = 0;
    size_t next_ = 0;    // ring position the consumer reads next
    bool held_ = false;  // the consumer holds ring_[next_]
    moodycamel::LightweightSemaphore free_slots_{0, /*maxSpins=*/0};
    moodycamel::BlockingConcurrentQueue<uint32_t> to_inflate_;
    std::atomic<bool> stop_{false};
    std::thread reader_;
    std::vector<std::thread> inflaters_;
    int inflater_count_ = 0;

    void read_ahead() {
        bool after_eof_marker = false;
        for (size_t seq = 0;; ++seq) {
            free_slots_.wait();
            if (stop_.load()) break;
            const size_t pos = seq % ring_size_;
            bgzf_block &b = ring_[pos].block;
            bgzf_read_block(src_, b, after_eof_marker);
            if (b.end || b.error) {
                ring_[pos].ready.signal();
                break;
            }
            after_eof_marker = b.eof_marker;
            to_inflate_.enqueue(static_cast<uint32_t>(pos));
        }
        for (int i = 0; i < inflater_count_; ++i) {
            to_inflate_.enqueue(kStopWorker);
        }
    }

    void inflate_loop() {
        libdeflate_decompressor *d = libdeflate_alloc_decompressor();
        uint32_t pos = 0;
        while (true) {
            to_inflate_.wait_dequeue(pos);
            if (pos == kStopWorker) break;
            bgzf_block &b = ring_[pos].block;
            if (d) {
                bgzf_inflate_block(d, b);
            } else {
                b.error = "out of memory";
            }
            ring_[pos].ready.signal();
        }
        libdeflate_free_decompressor(d);
    }

    // The next block in file order: read and inflated here without a ring,
    // else waited for in the ring after handing back the previous slot.
    bgzf_block &take_block() {
        if (!ring_) {
            bgzf_block &b = *inline_block_;
            bgzf_read_block(src_, b, after_eof_marker_);
            if (!b.end && !b.error) {
                after_eof_marker_ = b.eof_marker;
                bgzf_inflate_block(inline_decompressor_, b);
            }
            return b;
        }
        if (held_) {
            next_ = (next_ + 1) % ring_size_;
            free_slots_.signal();
        }
        slot &s = ring_[next_];
        s.ready.wait();
        held_ = true;
        return s.block;
    }

   public:
    bgzf_stream() = default;
    bgzf_stream(const bgzf_stream &) = delete;
    bgzf_stream &operator=(const bgzf_stream &) = delete;
    ~bgzf_stream() { close(); }

    // Reads `fd` (not owned) from its first block; errors go to `error`.
    bool open(int fd, int threads, std::string &error) {
        src_.open(fd);
        error_ = &error;
        if (threads <= 0) {
            inline_decompressor_ = libdeflate_alloc_decompressor();
            inline_block_ = std::make_unique<bgzf_block>();
            return inline_decompressor_ != nullptr;
        }
        inflater_count_ = threads;
        ring_size_ = kBlocksPerThread * static_cast<size_t>(threads);
        ring_ = std::make_unique<slot[]>(ring_size_);
        free_slots_.signal(static_cast<ssize_t>(ring_size_));
        reader_ = std::thread([this]() { read_ahead(); });
        for (int i = 0; i < threads; ++i) {
            inflaters_.emplace_back([this]() { inflate_loop(); });
        }
        return true;
    }

    // Next non-empty block of decompressed bytes, valid until the next call.
    // False at the end of the stream or on error.
    bool next(const uint8_t *&p, size_t &n) {
        while (!done_) {
            const bgzf_block &b = take_block();
            if (b.end || b.error) {
                if (b.error) *error_ = b.error;
                done_ = true;
                break;
            }
            if (b.data_size == 0) continue;
            p = b.data;
            n = b.data_size;
            return true;
        }
        return false;
    }

    void close() {
        if (reader_.joinable()) {
            stop_.store(true);
            free_slots_.signal(static_cast<ssize_t>(ring_size_));
            reader_.join();
        }
        for (auto &t : inflaters_) t.join();
        inflaters_.clear();
        ring_.reset();
        libdeflate_free_decompressor(inline_decompressor_);
        inline_decompressor_ = nullptr;
    }
};

// Plain gzip, possibly several concatenated members, through zlib.
class gzip_stream {
    static constexpr size_t kBufferSize = 1u << 20;

    int fd_ = -1;
    std::string *error_ = nullptr;
    z_stream zs_{};
    bool zs_open_ = false;
    bool member_open_ = false;
    bool in_eof_ = false;
    std::unique_ptr<uint8_t[]> in_;
    std::unique_ptr<uint8_t[]> out_;

   public:
    gzip_stream() = default;
    gzip_stream(const gzip_stream &) = delete;
    gzip_stream &operator=(const gzip_stream &) = delete;
    ~gzip_stream() { close(); }

    // Reads `fd` (not owned) from its first member; errors go to `error`.
    bool open(int fd, std::string &error) {
        fd_ = fd;
        error_ = &error;
        in_ = std::make_unique<uint8_t[]>(kBufferSize);
        out_ = std::make_unique<uint8_t[]>(kBufferSize);
        // 16 + MAX_WBITS: gzip wrapper only.
        zs_open_ = inflateInit2(&zs_, 16 + MAX_WBITS) == Z_OK;
        return zs_open_;
    }

    // Bytes inflated before an error are still handed out; the error ends
    // the stream at the next call.
    bool next(const uint8_t *&p, size_t &n) {
        if (!error_->empty()) return false;
        zs_.next_out = out_.get();
        zs_.avail_out = static_cast<uInt>(kBufferSize);
        while (zs_.avail_out > 0 && error_->empty()) {
            if (zs_.avail_in == 0 && !in_eof_) {
                const int64_t got = read_full(fd_, in_.get(), kBufferSize);
                if (got < 0) {
                    *error_ = "read error";
                    break;
                }
                in_eof_ = got == 0;
                zs_.next_in = in_.get();
                zs_.avail_in = static_cast<uInt>(got);
            }
            if (zs_.avail_in == 0) {
                if (member_open_) *error_ = "truncated gzip file";
                break;
            }
            if (!member_open_) {
                inflateReset(&zs_);
                member_open_ = true;
            }
            const int rc = inflate(&zs_, Z_NO_FLUSH);
            if (rc == Z_STREAM_END) {
                member_open_ = false;
            } else if (rc != Z_OK) {
                *error_ = "corrupt gzip data";
            }
        }
        n = kBufferSize - zs_.avail_out;
        p = out_.get();
        return n > 0;
    }

    void close() {
        if (zs_open_) inflateEnd(&zs_);
        zs_open_ = false;
    }
};

// One input file as decompressed bytes, block by block.
class input_stream {
    static constexpr size_t kRawChunkSize = 1u << 20;

    enum class compression { none, gzip, bgzf };

    compression compression_ = compression::none;
    int fd_ = -1;
    size_t raw_chunk_size_;
    std::unique_ptr<uint8_t[]> raw_;
    bgzf_stream bgzf_;
    gzip_stream gzip_;
    std::string error_;

   public:
    // `raw_chunk_size` is the block size of uncompressed input.
    explicit input_stream(size_t raw_chunk_size = kRawChunkSize)
        : raw_chunk_size_(raw_chunk_size) {}
    input_stream(const input_stream &) = delete;
    input_stream &operator=(const input_stream &) = delete;
    ~input_stream() { close(); }

    // `io_threads` inflate threads for BGZF; other inputs ignore it.
    bool open(const std::string &path, int io_threads) {
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) {
            error_ = std::string("cannot open: ") + strerror(errno);
            return false;
        }
#if defined(POSIX_FADV_SEQUENTIAL)
        posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
        uint8_t head[kBgzfHeaderSize];
        const ssize_t got = pread(fd_, head, sizeof(head), 0);
        if (got < 0) {
            error_ = std::string("read error: ") + strerror(errno);
            return false;
        }
        const size_t n = static_cast<size_t>(got);
        if (is_bgzf_header(head, n)) {
            compression_ = compression::bgzf;
            if (!bgzf_.open(fd_, io_threads, error_)) {
                error_ = "out of memory";
                return false;
            }
            return true;
        }
        if (is_gzip_magic(head, n)) {
            compression_ = compression::gzip;
            if (!gzip_.open(fd_, error_)) {
                error_ = "cannot initialise zlib";
                return false;
            }
            return true;
        }
        compression_ = compression::none;
        raw_ = std::make_unique<uint8_t[]>(raw_chunk_size_);
        return true;
    }

    // Next non-empty block of bytes, valid until the next call. False at the
    // end of the stream or on error.
    bool next(const uint8_t *&p, size_t &n) {
        if (compression_ == compression::bgzf) return bgzf_.next(p, n);
        if (compression_ == compression::gzip) return gzip_.next(p, n);
        if (!error_.empty()) return false;
        const int64_t got = read_full(fd_, raw_.get(), raw_chunk_size_);
        if (got < 0) {
            error_ = std::string("read error: ") + strerror(errno);
            return false;
        }
        p = raw_.get();
        n = static_cast<size_t>(got);
        return n > 0;
    }

    // Why the stream stopped early; empty after a clean end.
    const std::string &error() const { return error_; }

    void close() {
        bgzf_.close();
        gzip_.close();
        if (fd_ >= 0) {
            ::close(fd_);
            fd_ = -1;
        }
    }
};

enum class input_format { bam, fastx };

// Content sniffing (magic bytes), not file extension. After decompression:
// "BAM\1" is BAM; '>' or '@' (after blank lines) is FASTA/FASTQ, unless '@'
// starts a SAM header record. SAM, CRAM and anything else are rejected. An
// empty file reads as FASTA/FASTQ with no records.
inline std::expected<input_format, std::string> detect_input_format(
    const std::string &path) {
    // A record marker plus the four bytes that tell a SAM header from FASTQ.
    constexpr size_t kSniffAfterMarker = 4;
    constexpr const char *kBlank = " \t\r\n";
    const std::string kConvert =
        "; convert it to BAM or FASTQ first, e.g. samtools view -b -o "
        "reads.bam reads.sam";
    input_stream in;
    if (!in.open(path, 0)) return std::unexpected(in.error());
    std::string head;
    auto enough = [&]() {
        const size_t first = head.find_first_not_of(kBlank);
        return first != std::string::npos &&
               head.size() >= first + kSniffAfterMarker;
    };
    const uint8_t *p = nullptr;
    size_t n = 0;
    while (!enough() && in.next(p, n)) {
        head.append(reinterpret_cast<const char *>(p), n);
    }
    // An error after some bytes is the reader's to report, with the record
    // it hit.
    if (head.empty() && !in.error().empty()) {
        return std::unexpected(in.error());
    }
    if (head.starts_with(std::string_view("BAM\1", 4))) {
        return input_format::bam;
    }
    if (head.starts_with("CRAM")) {
        return std::unexpected("CRAM input is not supported" + kConvert);
    }
    const size_t first = head.find_first_not_of(kBlank);
    if (first == std::string::npos) return input_format::fastx;
    if (head[first] == '>') return input_format::fastx;
    if (head[first] == '@') {
        for (const char *type : {"@HD\t", "@SQ\t", "@RG\t", "@PG\t", "@CO\t"}) {
            if (head.compare(first, 4, type) == 0) {
                return std::unexpected("SAM input is not supported" + kConvert);
            }
        }
        return input_format::fastx;
    }
    return std::unexpected(
        "unrecognised input format (expected BAM, FASTA or FASTQ); SAM or "
        "CRAM input needs converting first, e.g. samtools view -b -o "
        "reads.bam reads.sam");
}

#endif
