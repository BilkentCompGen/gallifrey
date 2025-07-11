#pragma once
#ifndef _BGZF_H_
#define _BGZF_H_
#include <string>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

struct BAMReader {
    std::string filename;
    samFile *fp;
    sam_hdr_t *hdr;

    BAMReader(const std::string &filename, int threads = 2)
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

#endif