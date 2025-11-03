#pragma once
#ifndef _BGZF_H_
#define _BGZF_H_
#include <string>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include "endian.h"
struct BAMReader {
    std::string filename;
    samFile *fp;
    sam_hdr_t *hdr;

    BAMReader(const std::string &filename, int threads [[maybe_unused]] = 2)
        : fp{sam_open(filename.c_str(), "rb")} {
        this->filename = filename;
        // bgzf_mt(fp->fp.bgzf, 2, 0);
        // hts_set_threads(fp,threads);
        if (!(hdr = sam_hdr_read(fp))) {
            throw "Cannot read header from file!";
        }
    }
    ~BAMReader() {
        if (hdr) {
            sam_hdr_destroy(hdr);
        }
        if (fp) {
            sam_close(fp);
        }
    }

    bool read(bam1_t *bamdata) { return sam_read1(fp, hdr, bamdata) >= 0; }
    off_t tell() {
        if (fp->is_bgzf) {
            return bgzf_tell(fp->fp.bgzf);
        }
        return -199;
    }
    int seek(off_t offset) {
        if (fp->is_bgzf) {
            // return hseek(fp->fp.bgzf->fp, offset, SEEK_SET);
            return bgzf_seek(fp->fp.bgzf, offset, SEEK_SET);
        } else
            return -199;
    }
    bam1_t *init_data() {
        bam1_t *data = bam_init1();
        if (data == NULL) {
            throw "Not enough memory";
        }
        return data;
    }
};


struct Chunk {
    int64_t vstart;
    int64_t vend;
};

struct Block {
    off_t start;
    uint32_t len;
    constexpr int64_t vstart() const noexcept {
        return static_cast<int64_t>(start) << 16;
    }
    constexpr int64_t vend() const noexcept {
        return (static_cast<int64_t>(start + len)) << 16;
    }
};


class BlockIterator{
public:
    BlockIterator(const std::string& filename) {
        hfp_ = hopen(filename.c_str(), "r");
        if (!hfp_) {
            std::perror("hopen");
            valid_ = false;
            return;
        }
        if (hseek(hfp_, 0, SEEK_SET) < 0) {
            std::perror("hseek");
            valid_ = false;
        }
    }

    ~BlockIterator() {
        if (hfp_) {
            int result = hclose(hfp_);
            (void)result; // Suppress unused result warning
        }
    }

    bool next(Block& block) {
        if (!valid_) return false;
        off_t start = htell(hfp_);
        unsigned char hdr[18];
        ssize_t got = hread(hfp_, hdr, sizeof(hdr));
        if (got < static_cast<ssize_t>(sizeof(hdr))) return false;
        if (hdr[0]!=0x1f || hdr[1]!=0x8b || hdr[2]!=0x08 || hdr[3]!=0x04) return false;
        uint16_t bsize = le16toh(*reinterpret_cast<uint16_t*>(hdr + 16));
        uint32_t block_len = bsize + 1;
        block.start = start;
        block.len   = block_len;
        if (hseek(hfp_, start + block_len, SEEK_SET) < 0) {
            valid_ = false;
            return false;
        }
        return true;
    }

private:
    hFILE*  hfp_    = nullptr;
    bool    valid_  = true;
};

#endif