#pragma once
#ifndef LCMER_H
#define LCMER_H
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <gtl/vector.hpp>
#include <iterator>

#include "ext/lcptools/lcptools_ho.h"

// Makes sure we have enough scratch space
template <class T, double GF = 2>
struct charsi_base {
    T *data;
    size_t capacity;
    charsi_base() : data{NULL}, capacity{0} {}
    // Rule of 3 just delete
    charsi_base(const charsi_base &) = delete;
    charsi_base &operator=(const charsi_base &) = delete;
    // Rule of 5
    charsi_base(charsi_base &&other) noexcept
        : data{other.data}, capacity{other.capacity} {
        other.data = nullptr;
        other.capacity = 0;
    }
    charsi_base &operator=(charsi_base &&other) noexcept {
        if (this != &other) {
            free(data);
            data = other.data;
            capacity = other.capacity;
            other.data = nullptr;
            other.capacity = 0;
        }
        return *this;
    }

    ~charsi_base() { free(data); }
    int require(size_t cap, bool keep_when_realloc = false) {
        if (cap > capacity) {
            capacity = cap * GF;
            if (keep_when_realloc) {
                data = (decltype(data))realloc(data, capacity * sizeof(T));
            } else {
                free(data);
                data = (decltype(data))malloc(capacity * sizeof(T));
            }
        }
        return data == NULL;
    }
};
using charsi = charsi_base<uint8_t, 2.0>;
int onehotto2bitDNA2(const uint8_t *data, char *ret);
struct sv {
    char *from;
    char *to;
    sv() : from(0), to(0) {}
    sv(char *from, char *to) : from(from), to(to) {}

    void expand_until(char c) {
        while (*to != 0 && *to != c) {
            to++;
        }
    }
    void skip_next(char c) {
        expand_until(c);
        if (*to) ++to;
        from = to;
    }
    template <int N, char DEL>
    void get() {
        for (int i = 0; i < N; ++i) {
            skip_next(DEL);
        }
        expand_until(DEL);
    }
    // operator std::string() { return {from, to}; }
};

struct lcp_read {
    gtl::vector<uint32_t> cores;
    gtl::vector<std::pair<uint32_t, uint32_t>> pos;
    const uint8_t *data;
    void reserve(size_t size) {
        cores.reserve(size);
        pos.reserve(size);
    }
    void clear() {
        cores.clear();
        pos.clear();
        data = NULL;
    }
};

template <int N>
struct lcmer {
    // static uint32_t MaskedMurmur3(const std::array<uint32_t, N> &key,
    //                               uint32_t seed) {
    //     //        const uint8_t *data = (const uint8_t *)key;

    //     uint32_t h1 = seed;

    //     const uint32_t c1 = 0xcc9e2d51;
    //     const uint32_t c2 = 0x1b873593;

    //     // Body: Process blocks of 4 bytes at a time
    //     const uint32_t *blocks = key.data() + N;

    //     for (int i = -N; i; i++) {
    //         uint32_t k1 = lcmer<N>::mask[i + N] ? blocks[i] : 0;

    //         k1 *= c1;
    //         k1 = (k1 << 15) | (k1 >> (32 - 15));
    //         k1 *= c2;

    //         h1 ^= k1;
    //         h1 = (h1 << 15) | (h1 >> (32 - 15));
    //         h1 = h1 * 5 + 0xe6546b64;
    //     }
    //     // Finalization: Mix the hash to ensure the last few bits are fully
    //     // mixed
    //     // h1 ^= (N*4);

    //     /* fmix32 */
    //     h1 ^= h1 >> 16;
    //     h1 *= 0x85ebca6b;
    //     h1 ^= h1 >> 13;
    //     h1 *= 0xc2b2ae35;
    //     h1 ^= h1 >> 16;
    //     return h1;
    // }

    std::array<uint32_t, N> data;

    lcmer() {}
    lcmer(const lcmer<N> &other) {
        memcpy(this->data.data(), other.data.data(), N * sizeof(uint32_t));
    }
    lcmer(std::forward_iterator auto it) {
        std::copy(it, it + N, std::begin(data));
    }

    lcmer<N> operator=(const lcmer<N> &other) {
        memcpy(this->data.data(), other.data.data(), N * sizeof(uint32_t));
        return *this;
    }
    bool operator==(const lcmer &other) const {
        return memcmp(data.data(), other.data.data(), N * sizeof(uint32_t)) ==
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

template <typename INT_TYPE = uint32_t, INT_TYPE BIT_WIDTH = 2>
struct KmerViewOnCrumbs {
    using int_type = INT_TYPE;
    int_type k;
    int_type KMASK;
    int_type value;

    static constexpr INT_TYPE BIT_MASK = -1;
    // static constexpr INT_TYPE BIT_MASK = ((INT_TYPE)1 << BIT_WIDTH) - 1;

    KmerViewOnCrumbs(int_type k)
        : k{k}, KMASK{((int_type)1 << (BIT_WIDTH * k)) - 1}, value{0} {}

    int fill(const uint8_t *data_address) {
        int_type i;
        for (i = 0; i < k; i++) {
            move_get(data_address, i);
        }
        return i;
    }

    int fill_rle(const uint8_t *data_address, int_type cap) {
        int_type i;
        int_type it = move_get_rle(-1, data_address, 0, cap);
        for (i = 0; i < k && it < cap; i++) {
            it = move_get_rle(data_address[it - 1], data_address, it, cap);
        }
        return it;
    }

    int_type move_get(const uint8_t *data_address, int_type index) {
        value <<= BIT_WIDTH;
        value &= KMASK;
        value |= (data_address[index] & BIT_MASK);
        return index + 1;
    }

    int_type move_get_rle2(uint8_t prev, const uint8_t *data_address,
                           int_type index, int_type cap) {
        while (index < cap && prev == data_address[index]) {
            index++;
        }
        value <<= BIT_WIDTH;
        value &= KMASK;
        value |= (data_address[index] & BIT_MASK);
        return index + 1;
    }
    int_type move_get_rle(uint8_t prev, const uint8_t *data_address,
                          int_type index, int_type cap) {
        const uint8_t *p = data_address + index;
        const uint8_t *end = data_address + cap;

        while (p < end && *p == prev) {
            ++p;
        }
        if(p==end) { return cap;}
        value = ((value << BIT_WIDTH) & KMASK) | (*p & BIT_MASK);
        return static_cast<int_type>((p - data_address) + 1);
    }
    constexpr int_type get() const noexcept { return value; }
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
        return value;
    }

    constexpr int_type get() const noexcept { return value; }
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
        return value;
    }

    constexpr int_type get() const noexcept { return value; }
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

void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &c);
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &c);

inline uint64_t parse1_bam_syncmers(const uint8_t *data, int64_t len,
                                    struct lps *str,
                                    const auto &syncmer_functor) {
    int core_index = 0;
    KmerViewOnCrumbs<uint64_t, 2> kmer_view{syncmer_functor.k()};
    int64_t it1 = kmer_view.fill_rle(data, len);
    // find lcp cores
    // if(str->capacity < )
    if (str->capacity < len * CONSTANT_FACTOR) {
        str->capacity = 2 * len * CONSTANT_FACTOR;
        str->cores = (struct core *)realloc(
            str->cores, (str->capacity) * sizeof(struct core));
    }
    for (; it1 + SYNCMER_K < len;) {
        // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
        if (syncmer_functor(kmer_view.get())) {
            // if (core_index >= str->capacity){
            //     str->cores = (struct core *)realloc(str->cores,
            //     (str->capacity * 2 + 1) * sizeof(struct core)); str->capacity
            //     = str->capacity * 2 + 1;
            // }
            init_core_kmer<uint32_t>(&(str->cores[core_index]), kmer_view.get(),
                                     it1, it1 + SYNCMER_K);
            core_index++;
        }
        it1 = kmer_view.move_get_rle(data[it1 - 1], data, it1, len);
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}

uint64_t parse1_bam_lmax(const uint8_t *data, int64_t len, struct core *cores);

/// @brief Initializes an lps struct with a given sequence in the bam compressed
// void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len);
inline void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len,
                         const auto &syncmer_functor) {
    lps_ptr->level = 1;
    if (lps_ptr->cores == nullptr) {
        lps_ptr->size = 0;
        lps_ptr->capacity = (len / CONSTANT_FACTOR);
        lps_ptr->cores = (struct core *)malloc((len / CONSTANT_FACTOR) *
                                               sizeof(struct core));
    } else {
        if (lps_ptr->capacity < len / CONSTANT_FACTOR) {
            lps_ptr->capacity = len / CONSTANT_FACTOR;
            lps_ptr->cores = (struct core *)realloc(
                lps_ptr->cores, (len / CONSTANT_FACTOR) * sizeof(struct core));
        }
    }

    //    lps_ptr->size = parse1_bam_syncmers_nthash(data, len, lps_ptr,
    //    syncmer_functor);
    lps_ptr->size = parse1_bam_syncmers(data, len, lps_ptr, syncmer_functor);
}

inline void parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                        int lcp_level, lps *str,
                                        const auto &syncmer_functor,
                                        charsi &chs, lcp_read &lr,
                                        bool reversed = false) {
    // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

    // onehotto2bitDNA_bivec_chs(data, len, chs);  // TODO Make this not
    // allocate work
    //  on a scratch space per thread.
    //
    chs.require(len);
    //    nibble2base((uint8_t *)data, (char *)chs.data, (int)len);
    if (reversed) {
        onehotto2bitDNA_bivec_chs_reverse(data, len, chs);
    } else {
        onehotto2bitDNA_bivec_chs(data, len, chs);
    }

    init_lps_bam(str, chs.data, len, syncmer_functor);
    lps_deepen(str, lcp_level);

    for (int i = 0; i < str->size; ++i) {
        lr.cores.emplace_back(str->cores[i].label);
        lr.pos.emplace_back(str->cores[i].start, str->cores[i].end);
    }
    lr.data = data;
}

lcp_read parse_cores_from_kseq(const char *data, size_t len, int lcp_level);
#endif
