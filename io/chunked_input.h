#ifndef _GALLIFREY_IO_CHUNKED_INPUT_H_
#define _GALLIFREY_IO_CHUNKED_INPUT_H_

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "bam_chunks.h"
#include "bam_reader.h"
#include "chunks.h"
#include "fastx.h"
#include "fastx_chunks.h"
#include "input_stream.h"

// An input file cut into record-aligned chunks, each read on its own by a
// chunk_reader.
struct chunked_input {
    std::string path;
    std::unique_ptr<input_stream> stream;  // keeps the file open
    uint64_t size = 0;
    input_format format = input_format::fastx;
    bool bgzf = false;
    fastx_kind kind = fastx_kind::unknown;  // FASTA/FASTQ records
    std::vector<chunk> chunks;
};

// Why cut_into_chunks reads an input whole.
inline constexpr const char *kNotRegularFile = "not a regular file";
inline constexpr const char *kNotBgzf = "not BGZF compressed";
inline constexpr const char *kBamNotCut =
    "not an unaligned BAM file that can be cut (aligned, no end-of-file "
    "marker or too small)";
inline constexpr const char *kFastxNotCut =
    "no record starts to cut at (wrapped FASTQ, one record, no end-of-file "
    "marker or too small)";

// Cuts the input `path` (open as `stream`, status `st`, format `format`)
// into up to `parts` chunks, described in `chunked`. Returns why it is read
// whole instead, or nullptr when it was cut.
inline const char *cut_into_chunks(input_stream &stream, const struct stat &st,
                                   const std::string &path, input_format format,
                                   size_t parts, chunked_input &chunked) {
    if (is_stdin_path(path) || !S_ISREG(st.st_mode)) return kNotRegularFile;
    const bool is_bam = format == input_format::bam;
    if ((is_bam || stream.compressed()) && !stream.bgzf_compressed()) {
        return kNotBgzf;
    }
    chunked.path = path;
    chunked.size = static_cast<uint64_t>(st.st_size);
    chunked.format = format;
    const int fd = stream.descriptor();
    chunked.bgzf = stream.bgzf_compressed();
    if (is_bam) {
        chunked.chunks = cut_bam_file(fd, chunked.size, parts);
        if (chunked.chunks.empty()) return kBamNotCut;
        return nullptr;
    }
    chunked.chunks =
        cut_fastx_file(fd, chunked.size, chunked.bgzf, parts, chunked.kind);
    if (chunked.chunks.empty()) return kFastxNotCut;
    return nullptr;
}

// The file of a chunked input opened again for one chunk, so that the kernel
// and network file systems see one sequential reader per chunk and read
// ahead for each (on a shared descriptor the interleaved offsets look random
// and read-ahead stops). The input's own descriptor when the path no longer
// names the same file.
class chunk_descriptor {
    int fd_ = -1;
    int shared_;

   public:
    explicit chunk_descriptor(const chunked_input &input)
        : shared_(input.stream->descriptor()) {
        fd_ = ::open(input.path.c_str(), O_RDONLY);
        if (fd_ < 0) return;
        struct stat a, b;
        if (fstat(fd_, &a) != 0 || fstat(shared_, &b) != 0 ||
            !same_file(a, b)) {
            ::close(fd_);
            fd_ = -1;
            return;
        }
#if defined(POSIX_FADV_SEQUENTIAL)
        posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    }
    ~chunk_descriptor() {
        if (fd_ >= 0) ::close(fd_);
    }
    chunk_descriptor(const chunk_descriptor &) = delete;
    chunk_descriptor &operator=(const chunk_descriptor &) = delete;
    int fd() const {
        if (fd_ >= 0) return fd_;
        return shared_;
    }
};

// The records of chunk `c` of a chunked input, in file order, like the
// readers of whole files. Every record of the file is returned by exactly one
// of its chunks, or the reading of some chunk fails: a FASTA/FASTQ record of
// a shape the cut did not allow, or a BAM chunk that does not end at a record
// end, means a cut inside a record.
class chunk_reader {
    const chunked_input &input_;
    bool last_;
    chunk_descriptor fd_;
    chunk_source src_;
    std::optional<bam_reader<chunk_source>> bam_;
    std::optional<fastx_reader<chunk_source>> fastx_;
    std::string misfit_;

    // The error of a FASTA/FASTQ record of a shape the cut did not allow.
    std::string misfit_error(seq_record &rec) {
        // The file's last record lacking its quality (or all of it): the file
        // was cut short, which no reading mode should accept.
        const bool truncated =
            last_ &&
            (fastx_->truncated_record() ||
             (input_.kind == fastx_kind::fastq && !fastx_->fastq_record())) &&
            !fastx_->next(rec) && !fastx_->failed();
        if (truncated) return "it ends inside its last record (truncated file)";
        return "it mixes record layouts (wrapped FASTQ, or FASTQ inside "
               "FASTA), so it cannot be read in chunks; read it with "
               "--file-reading-mode sequential";
    }

   public:
    chunk_reader(const chunked_input &input, size_t c)
        : input_(input),
          last_(c + 1 == input.chunks.size()),
          fd_(input),
          src_(fd_.fd(), input.size, input.bgzf, input.chunks[c]) {
        if (input.format == input_format::bam) {
            bam_.emplace(src_);
        } else {
            fastx_.emplace(src_);
        }
    }

    // Fills `rec` with the next record; false at the end of the chunk or on
    // error (see error()).
    bool next(seq_record &rec) {
        if (bam_) return bam_->next(rec);
        if (!misfit_.empty() || !fastx_->next(rec)) return false;
        if (record_fits_chunk(input_.kind, *fastx_)) return true;
        misfit_ = misfit_error(rec);
        return false;
    }

    // Why the chunk could not be read; empty after a clean end.
    std::string error() const {
        if (!src_.error().empty()) return src_.error();
        if (bam_ && bam_->failed()) {
            return "reading it in chunks failed (" + bam_->error() +
                   "); read it with --file-reading-mode sequential";
        }
        if (fastx_ && fastx_->failed()) {
            return "malformed or unreadable FASTA/FASTQ";
        }
        return misfit_;
    }
};

#endif
