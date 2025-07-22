#pragma once
#include <iterator>
#ifndef LCMER_H
#define LCMER_H
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <gtl/vector.hpp>

#include "config.h"
#include "ext/lcptools/lcptools_ho.h"

uint32_t MurmurHash3_32(const void *key, int len, uint32_t seed);

struct lcp_read {
    gtl::vector<uint32_t> cores;
    gtl::vector<std::pair<uint32_t, uint32_t>> pos;
    void reserve(size_t size) {
        cores.reserve(size);
        pos.reserve(size);
    }
};

template <int N>
struct lcmer {
    std::array<uint32_t, N> data;

    lcmer() {}
    lcmer(const lcmer<N> &other) {
        memcpy(this->data.data(), other.data.data(), N * sizeof(uint32_t));
    }
    lcmer(std::forward_iterator auto it) {
        std::copy(it, it + N, std::begin(data));
    }

    bool operator==(const lcmer &other) const {
        return memcmp(data.data(), other.data.data(), 16 * sizeof(uint32_t)) ==
               0;
        // for (int i = 0; i < N; ++i) {
        //     if (data[i] != other.data[i]) return false;
        // }
        // return true;
    }
    bool operator<(const lcmer &other) const {
        for (int i = 0; i < N; ++i) {
            if (data[i] != other.data[i]) {
                if (data[i] < other.data[i]) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }
    void print_name(FILE *out) const {
        fprintf(out, "%0x", data[0]);
        for (int z = 1; z < N; ++z) {
            fprintf(out, "%0x", data[z]);
        }
    }

    void set0(uint32_t new_core) { data[N - 1] = new_core; }
    void setN(uint32_t new_core) { data[0] = new_core; }
    static lcmer<N> left_shift(lcmer<N> l) {
        for (int i = 1; i < N; ++i) {
            l.data[i - 1] = l.data[i];
        }
        return l;
    }

    static lcmer<N> right_shift(lcmer<N> l) {
        for (int i = N - 1; i > 0; --i) {
            l.data[i] = l.data[i - 1];
        }
        return l;
    }

    uint32_t update_lcmer(uint32_t new_core) {
        uint32_t old_core = data[0];
        for (int i = 1; i < N; ++i) {
            data[i - 1] = data[i];
        }
        data[N - 1] = new_core;
        return old_core;
    }
};

template <int N>
struct std::hash<lcmer<N>> {
    std::size_t operator()(const lcmer<N> &l) const noexcept {
        return MurmurHash3_32(&l.data, N * sizeof(uint32_t), 12312342);
    }
};

// https://github.com/ruanjue/wtdbg2/blob/b77c5657c8095412317e4a20fe3668f5bde6b1ac/dna.h#L108
template <class T>
T reverse_complement_kmer(T kmer, int K) {
    uint32_t shift_count = (sizeof(T) * 8) - 2 * K;
    kmer = ~kmer;
    kmer = ((kmer & 0x3333333333333333LLU) << 2) |
           ((kmer & 0xCCCCCCCCCCCCCCCCLLU) >> 2);
    kmer = ((kmer & 0x0F0F0F0F0F0F0F0FLLU) << 4) |
           ((kmer & 0xF0F0F0F0F0F0F0F0LLU) >> 4);
    if constexpr (sizeof(T) == 8) {
        kmer = __builtin_bswap64(kmer);
    } else if constexpr (sizeof(T) == 4) {
        kmer = __builtin_bswap32(kmer);
    } else if constexpr (sizeof(T) == 2) {
        kmer = __builtin_bswap16(kmer);
    }
    kmer = kmer >> shift_count;
    return kmer;
}

template <int _K, typename INT_TYPE = uint32_t, int BIT_WIDTH = 2>
struct KmerViewOnCrumbs {
    using int_type = INT_TYPE;
    static constexpr int k = _K;
    int_type value;
    static constexpr int_type KMASK = ((int_type)1 << (BIT_WIDTH * k)) - 1;
    KmerViewOnCrumbs() : value{0} {}

    int fill(const uint8_t *data_address) {
        int i;
        for (i = 0; i < k; i++) {
            move_get(data_address, i);
        }
        return i;
    }
    int_type move_get(const uint8_t *data_address, int index) {
        value <<= BIT_WIDTH;
        if constexpr (KMASK != -1) value &= KMASK;
        value |= data_address[index];

        value = value & KMER_SPACE_MASK;
        return value;
    }

    constexpr int_type get() const noexcept { return value & KMER_SPACE_MASK; }
};

template <int _K, typename INT_TYPE = uint32_t, int BIT_WIDTH = 2>
struct KmerViewOnNibbles {
    static constexpr int tab[] = {4, 0, 1, 4, 2, 4, 4, 4,
                                  3, 4, 4, 4, 4, 4, 4, 4};

    using int_type = INT_TYPE;
    static constexpr int k = _K;
    int_type value;
    static constexpr int_type KMASK = ((int_type)1 << (BIT_WIDTH * k)) - 1;
    KmerViewOnNibbles() : value{0} {}

    int fill(const uint8_t *data_address) {
        int i;
        for (i = 0; i < k; i++) {
            move_get(data_address, i);
        }
        return i;
    }
    int_type move_get(const uint8_t *data_address, int index) {
        value <<= BIT_WIDTH;
        if constexpr (KMASK != -1) value &= KMASK;
        if constexpr (BIT_WIDTH == 2) {
            value |= tab[bam_seqi(data_address, index)];
        } else {
            value |= bam_seqi(data_address, index);
        }
        value = value & KMER_SPACE_MASK;
        return value;
    }

    constexpr int_type get() const noexcept { return value & KMER_SPACE_MASK; }
};

template <int _K, typename INT_TYPE = uint32_t, int BIT_WIDTH = 2>
struct KmerViewOnBytes {
    static constexpr const int tab[256] = {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4,
        4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
    using int_type = INT_TYPE;
    static constexpr int k = _K;
    int_type value;
    static constexpr int_type KMASK = ((int_type)1 << (BIT_WIDTH * k)) - 1;
    KmerViewOnBytes() : value{0} {}

    int fill(const char *data) {
        int i;
        for (i = 0; i < k; i++) {
            move_get(data, i);
        }
        return i;
    }
    int_type move_get(const char *data, int index) {
        value <<= BIT_WIDTH;
        if constexpr (KMASK != -1) value &= KMASK;
        if constexpr (BIT_WIDTH == 2) {
            value |= tab[static_cast<int>(data[index])];
        } else {
            value |= data[index];
        }
        value = value & KMER_SPACE_MASK;
        return value;
    }

    constexpr int_type get() const noexcept { return value & KMER_SPACE_MASK; }
};

void init_core1bam(struct core *cr, const uint8_t *data, uint64_t distance,
                   uint64_t si, uint64_t ei);

template <typename _T, int WIDTH = 2>
void init_core_kmer(struct core *cr, _T Kval, uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= Kval;
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = WIDTH * SYNCMER_K;
}

template <typename _T>
void init_core_kmer_xor(struct core *cr, _T Kval, uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= (Kval ^ reverse_complement_kmer(Kval, SYNCMER_K));
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * SYNCMER_K;
}

uint8_t *onehotto2bitDNA(const uint8_t *data, size_t len);
uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len);

uint64_t parse1_bam_syncmers(const uint8_t *data, int64_t len,
                             struct core *cores);

uint64_t parse1_bam_lmax(const uint8_t *data, int64_t len, struct core *cores);

/// @brief Initializes an lps struct with a given sequence in the bam compressed
void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len);

lcp_read parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                     int lcp_level);

lcp_read parse_cores_from_kseq(const char *data, size_t len, int lcp_level);
#endif