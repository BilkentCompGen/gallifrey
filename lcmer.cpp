#include "lcmer.h"

#include "ext/doctest.h"
#include "gtl/phmap.hpp"
#include "syncmer.h"

#include "ext/lcptools/lcptools_ho.h" 

void init_core1bam(struct core *cr, const uint8_t *data, uint64_t distance,
                   uint64_t si, uint64_t ei) {
    cr->start = si;
    cr->end = ei;
    cr->label = 0;
    cr->label |= ((distance - 2) << 12);
    cr->label |= (bam_seqi(data, si) << 8);
    cr->label |= (bam_seqi(data, si + distance - 2) << 4);
    cr->label |= (bam_seqi(data, si + distance - 1));
    cr->bit_rep = 0x8000000000000000 | cr->label;
    cr->bit_size = 2 * distance;
}

syncmer::syncmer_precomputed<
    syncmer::SYNCMER_RULE<SYNCMER_K, SYNCMER_S, SYNCMER_T>>
    table;
uint64_t parse1_bam_syncmers(const uint8_t *data, int64_t len,
                             struct lps *str) {
    int core_index = 0;
    KmerViewOnCrumbs<SYNCMER_K, uint32_t, 2> kmer_view;
    int64_t it1 = kmer_view.fill(data);
    // find lcp cores
    for (; it1 + SYNCMER_K < len; it1++) {
        // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
#if SYNCMER_RULE == DensitySyncmerRule
        if ( kmer_view.get() < ((double)SYNCMER_S / SYNCMER_T)) 
#else
        if ( table[kmer_view.get()]) 
#endif
        {
            if (core_index >= str->capacity){
                str->cores = (struct core *)realloc(str->cores, (core_index * 2) * sizeof(struct core));
                str->capacity = core_index * 2;
            }
            init_core_kmer<uint32_t, 2>(&(str->cores[core_index]), kmer_view.get(),
                                        it1, it1 + SYNCMER_K);
            core_index++;
        }
        kmer_view.move_get(data, it1);
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}
uint64_t parse1_bam_syncmers_old(const uint8_t *data, int64_t len,
                                 struct core *cores) {
    int core_index = 0;
    KmerViewOnNibbles<SYNCMER_K, uint32_t, 2> kmer_view;
    int64_t it1 = kmer_view.fill(data);
    // find lcp cores
    for (; it1 + SYNCMER_K < len; it1++) {
        // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
        if (table[kmer_view.get()]) {
            init_core_kmer<uint32_t, 2>(&(cores[core_index]), kmer_view.get(),
                                        it1, it1 + SYNCMER_K);
            core_index++;
        }
        kmer_view.move_get(data, it1);
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}
uint64_t parse1_kseq_syncmers(const char *data, int64_t len,
                              struct core *cores) {
    int core_index = 0;
    KmerViewOnBytes<SYNCMER_K, uint32_t, 2> kmer_view;
    int64_t it1 = kmer_view.fill(data);
    // find lcp cores
    for (; it1 + SYNCMER_K < len; it1++) {
        // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
        if (table[kmer_view.get()]) {
            init_core_kmer<uint32_t, 2>(&(cores[core_index]), kmer_view.get(),
                                        it1, it1 + SYNCMER_K);
            core_index++;
        }
        kmer_view.move_get(data, it1);
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}

uint64_t parse1_bam_lmax(const uint8_t *data, int64_t len, struct core *cores) {
    int64_t it1 = 0;
    int64_t it2 = len;

    int core_index = 0;
    // int64_t last_invalid_char_index = -1;

    // find lcp cores
    for (; it1 + 2 < len; it1++) {
        if (bam_seqi(data, it1) == bam_seqi(data, it1 + 1)) {
            continue;
        }
        if (0 == it1) {
            continue;
        }

        // check for LMAX
        if (it1 + 3 < len && bam_seqi(data, it1) < bam_seqi(data, (it1 + 1)) &&
            bam_seqi(data, (it1 + 1)) > bam_seqi(data, (it1 + 2)) &&
            bam_seqi(data, (it1 - 1)) <= bam_seqi(data, (it1)) &&
            bam_seqi(data, (it1 + 2)) >= bam_seqi(data, (it1 + 3))) {
            // // check if there is any SSEQ cores left behind
            // if (it2 < it1 && last_invalid_char_index < it2 - 1) {
            //   init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 -
            //   1,
            //                 it1 + 1);
            //   core_index++;
            // }

            // create LMAX core
            it2 = it1 + 3;
            init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
            core_index++;

            continue;
        }
        // cout << "None\n";
    }
    // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
    return core_index;
}
void init_lps_bam_old(lps *lps_ptr, const uint8_t *data, size_t len) {
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->capacity = (len / CONSTANT_FACTOR);
    lps_ptr->cores =
        (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
    lps_ptr->size = parse1_bam_syncmers_old(data, len, lps_ptr->cores);
}
void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len) {
    lps_ptr->level = 1;
    if(lps_ptr->cores == nullptr){
        lps_ptr->size = 0;
        lps_ptr->capacity = (len / CONSTANT_FACTOR);
        lps_ptr->cores =
            (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
    }
    else{
        if (lps_ptr->capacity < len / CONSTANT_FACTOR){
            lps_ptr->capacity = len / CONSTANT_FACTOR;
            lps_ptr->cores = (struct core *)realloc(lps_ptr->cores, (len / CONSTANT_FACTOR) * sizeof(struct core));
        }
    }
    lps_ptr->size = parse1_bam_syncmers(data, len, lps_ptr);
}

void init_lps_kseq(lps *lps_ptr, const char *data, size_t len) {
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores =
        (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
    lps_ptr->size = parse1_kseq_syncmers(data, len, lps_ptr->cores);
}

uint8_t *onehotto2bitDNA(uint8_t *data, size_t len) {
    uint8_t *ret = (uint8_t *)malloc(len);
    for (size_t i = 0; i < len / 2; ++i) {
        ret[2 * i] = 27 - __builtin_clz(data[i]);
        ret[2 * i + 1] = __builtin_ctz(data[i]);
    }
    return ret;
}
uint8_t *onehotto2bitDNA_bivec(const uint8_t *data, size_t len) {
    uint8_t *ret = (uint8_t *)malloc(len);
#pragma clang loop vectorize(enable) interleave(enable)
    for (size_t i = 0; i < len / 2; ++i) {
        ret[2 * i] = 27 - __builtin_clz(data[i]);
        ret[2 * i + 1] = __builtin_ctz(data[i]);
    }
    if (len % 2 == 1) {
        ret[len - 1] = 27 - __builtin_clz(data[len / 2]);
    }
    return ret;
}

lcp_read parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                     int lcp_level, lps *str) {
    // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

    uint8_t *two_bit_encoded =
        onehotto2bitDNA_bivec(data, len);  // TODO Make this not allocate work
                                           // on a scratch space per thread.
    init_lps_bam(str, two_bit_encoded, len);
    lps_deepen(str, lcp_level);
    lcp_read lr;
    lr.reserve(str->size);
    for (int i = 0; i < str->size; ++i) {
        lr.cores.emplace_back(str->cores[i].label);
        lr.pos.emplace_back(str->cores[i].start, str->cores[i].end);
    }
    // free_lps(&str);
    free(two_bit_encoded);
    return lr;
}
lcp_read parse_cores_from_kseq(const char *s, size_t len, int lcp_level) {
    // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

    lps str;
    init_lps_kseq(&str, s, len);
    lps_deepen(&str, lcp_level);
    lcp_read lr;
    lr.reserve(str.size);
    for (int i = 0; i < str.size; ++i) {
        lr.cores.emplace_back(str.cores[i].label);
        lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
    }
    free_lps(&str);
    return lr;
}

// TESTS //

TEST_CASE("Nibbles to Crumbs") {
    gtl::vector<uint8_t> bytes;
    gtl::flat_hash_map<uint8_t, std::pair<uint8_t, uint8_t>> expected;
    static constexpr int tab[] = {4, 0, 1, 4, 2, 4, 4, 4,
                                  3, 4, 4, 4, 4, 4, 4, 4};

    for (int high : {1, 2, 4, 8}) {
        for (int low : {1, 2, 4, 8}) {
            uint8_t byte = (high << 4) | low;
            bytes.push_back(byte);
            expected[byte] = {tab[high], tab[low]};
        }
    }
    SUBCASE("Even") {
        uint8_t *conversion =
            onehotto2bitDNA_bivec(bytes.data(), 2 * bytes.size());
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(expected[bytes[i]].first == conversion[2 * i]);
            CHECK(expected[bytes[i]].second == conversion[2 * i + 1]);
        }
        free(conversion);
    }
    bytes.push_back((1 << 4));
    expected[bytes.back()] = {tab[1], 0};
    SUBCASE("Odd") {
        uint8_t *conversion =
            onehotto2bitDNA_bivec(bytes.data(), 2 * bytes.size() - 1);
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(expected[bytes[i]].first == conversion[2 * i]);
            CHECK(expected[bytes[i]].second == conversion[2 * i + 1]);
        }
        free(conversion);
    }
}

TEST_CASE("K-mer view on Crumbs") {
    gtl::vector<uint8_t> bytes{{0, 1, 2, 3, 0, 1, 2, 3}};
    SUBCASE("K==4") {
        static constexpr int K = 4;
        KmerViewOnCrumbs<K, uint32_t, 2> kview;
        int i = kview.fill(bytes.data());
        uint32_t expected =
            bytes[0] << 6 | bytes[1] << 4 | bytes[2] << 2 | bytes[3];
        for (; i < bytes.size(); ++i) {
            CHECK(expected == kview.get());
            kview.move_get(bytes.data(), i);
            expected <<= 2;
            expected |= bytes[i];
            expected &= 0xff;
        }
    }
    SUBCASE("K==5") {
        static constexpr int K = 5;
        KmerViewOnCrumbs<K, uint32_t, 2> kview;
        int i = kview.fill(bytes.data());
        uint32_t expected =
            bytes[0] << 8 | bytes[1] << 6 | bytes[2] << 4 | bytes[3] << 2 | bytes[4];
        for (; i < bytes.size(); ++i) {
            CHECK(expected == kview.get());
            kview.move_get(bytes.data(), i);
            expected <<= 2;
            expected |= bytes[i];
            expected &= 0x3ff;
        }
    }
}

TEST_CASE("Character table on kmerview on bytes"){
    using kv = KmerViewOnBytes<16, uint32_t, 2>;

    CHECK_EQ(kv::tab['A'], 0);
    CHECK_EQ(kv::tab['G'], 2);
    CHECK_EQ(kv::tab['C'], 1);
    CHECK_EQ(kv::tab['T'], 3);

    CHECK_EQ(kv::tab['a'], 0);
    CHECK_EQ(kv::tab['g'], 2);
    CHECK_EQ(kv::tab['c'], 1);
    CHECK_EQ(kv::tab['t'], 3);


    CHECK_EQ(kv::tab['u'], 3);
    CHECK_EQ(kv::tab['U'], 3);

    CHECK_EQ(kv::tab['N'], 4);
    CHECK_EQ(kv::tab['Z'], 4);
    
}