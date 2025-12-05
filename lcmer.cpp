#include "lcmer.h"

#include <cstdint>

#include "ext/doctest.h"
#include "ext/lcptools/lcptools_ho.h"
#include "gtl/phmap.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "nibble2base.hpp"

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

// syncmer::syncmer_precomputed<
//     syncmer::SYNCMER_RULE<SYNCMER_K, SYNCMER_S, SYNCMER_T>>
//     table;

// uint64_t parse1_bam_syncmers_old(const uint8_t *data, int64_t len,
//                                  struct core *cores) {
//     int core_index = 0;
//     KmerViewOnNibbles<SYNCMER_K, uint32_t, 2> kmer_view;
//     int64_t it1 = kmer_view.fill(data);
//     // find lcp cores
//     for (; it1 + SYNCMER_K < len; it1++) {
//         // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
//         if (table[kmer_view.get()]) {
//             init_core_kmer<uint32_t, 2>(&(cores[core_index]),
//             kmer_view.get(),
//                                         it1, it1 + SYNCMER_K);
//             core_index++;
//         }
//         kmer_view.move_get(data, it1);
//         // cout << "None\n";
//     }
//     // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
//     return core_index;
// }
// uint64_t parse1_kseq_syncmers(const char *data, int64_t len,
//                               struct core *cores) {
//     int core_index = 0;
//     KmerViewOnBytes<SYNCMER_K, uint32_t, 2> kmer_view;
//     int64_t it1 = kmer_view.fill(data);
//     // find lcp cores
//     for (; it1 + SYNCMER_K < len; it1++) {
//         // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
//         if (table[kmer_view.get()]) {
//             init_core_kmer<uint32_t, 2>(&(cores[core_index]),
//             kmer_view.get(),
//                                         it1, it1 + SYNCMER_K);
//             core_index++;
//         }
//         kmer_view.move_get(data, it1);
//         // cout << "None\n";
//     }
//     // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
//     return core_index;
// }

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
// void init_lps_bam_old(lps *lps_ptr, const uint8_t *data, size_t len) {
//     lps_ptr->level = 1;
//     lps_ptr->size = 0;
//     lps_ptr->capacity = (len / CONSTANT_FACTOR);
//     lps_ptr->cores =
//         (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
//     lps_ptr->size = parse1_bam_syncmers_old(data, len, lps_ptr->cores);
// }

// void init_lps_kseq(lps *lps_ptr, const char *data, size_t len) {
//     lps_ptr->level = 1;
//     lps_ptr->size = 0;
//     lps_ptr->cores =
//         (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
//     lps_ptr->size = parse1_kseq_syncmers(data, len, lps_ptr->cores);
// }
int onehotto2bitDNA2(const uint8_t *data, char *ret) {
    static char table[] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0;; ++i) {
        int clz = __builtin_clz(data[i]);
        int ctz = __builtin_ctz(data[i]);

        if (clz < 0 || ctz > 3) {
            ret[i] = 0;
            return i;
        }
        ret[2 * i] = table[27 - clz];
        ret[2 * i + 1] = table[ctz];
    }
    return -1;
}
uint8_t *onehotto2bitDNA(uint8_t *data, size_t len) {
    uint8_t *ret = (uint8_t *)malloc(len);
    for (size_t i = 0; i < len / 2; ++i) {
        ret[2 * i] = 27 - __builtin_clz(data[i]);
        ret[2 * i + 1] = __builtin_ctz(data[i]);
    }
    return ret;
}
void onehotto2bitDNA_bivec_chs(const uint8_t *data, size_t len, charsi &chs) {
    chs.require(len);

    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * i] = 27 - __builtin_clz(data[i]);
        chs.data[2 * i + 1] = __builtin_ctz(data[i]);
    }
    if (len % 2 == 1) {
        chs.data[len - 1] = 27 - __builtin_clz(data[len / 2]);
    }
}
void onehotto2bitDNA_bivec_chs_reverse(const uint8_t *data, size_t len,
                                       charsi &chs) {
    chs.require(len);
    size_t mid = len / 2 - 1;
    for (size_t i = 0; i < len / 2; ++i) {
        chs.data[2 * (mid - i)] = __builtin_clz(data[i]) - 24;
        chs.data[2 * (mid - i) + 1] = 3 - __builtin_ctz(data[i]);
    }
    if (len % 2 == 1) {
        chs.data[0] = __builtin_clz(data[len / 2]) - 24;
    }
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

// lcp_read parse_cores_from_kseq(const char *s, size_t len, int lcp_level) {
//     // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

//     lps str;
//     init_lps_kseq(&str, s, len);
//     lps_deepen(&str, lcp_level);
//     lcp_read lr;
//     lr.reserve(str.size);
//     for (int i = 0; i < str.size; ++i) {
//         lr.cores.emplace_back(str.cores[i].label);
//         lr.pos.emplace_back(str.cores[i].start, str.cores[i].end);
//     }
//     free_lps(&str);
//     return lr;
// }

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
    SUBCASE("Even reverse") {
        charsi c;
        c.require(2 * bytes.size());
        onehotto2bitDNA_bivec_chs_reverse(bytes.data(), 2 * bytes.size(), c);
        size_t bs = bytes.size() - 1;
        for (unsigned i = 0; i < bytes.size(); ++i) {
            CHECK(3 - expected[bytes[bs - i]].first == c.data[2 * i]);
            CHECK(3 - expected[bytes[bs - i]].second == c.data[2 * i + 1]);
        }
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
        KmerViewOnCrumbs<uint64_t, 2> kview(K);
        size_t i = kview.fill(bytes.data());
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
        KmerViewOnCrumbs<uint64_t, 2> kview(K);
        size_t i = kview.fill(bytes.data());
        uint32_t expected = bytes[0] << 8 | bytes[1] << 6 | bytes[2] << 4 |
                            bytes[3] << 2 | bytes[4];
        for (; i < bytes.size(); ++i) {
            CHECK(expected == kview.get());
            kview.move_get(bytes.data(), i);
            expected <<= 2;
            expected |= bytes[i];
            expected &= 0x3ff;
        }
    }
}

TEST_CASE("Character table on kmerview on bytes") {
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

TEST_CASE("BAM Read Parsing") {
    std::string input = "../tests/test_single.bam";
    samFile *fp = sam_open(input.c_str(), "r");
    if (fp == NULL) {
        FAIL("Failed to open BAM file:", input);
    }

    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    if (sam_read1(fp, hdr, rec) <= 0) {
        FAIL("Failed to read BAM record");
    }

    char *seq_gt = (char *)malloc(rec->core.l_qseq + 1);
    for (int i = 0; i < rec->core.l_qseq; ++i) {
        seq_gt[i] = seq_nt16_str[bam_seqi(bam_get_seq(rec), i)];
    }
    seq_gt[rec->core.l_qseq] = '\0';
    SUBCASE("Simd Nibble2base") {
        char *seq_simd = (char *)malloc(rec->core.l_qseq + 1);
        nibble2base(bam_get_seq(rec), seq_simd, rec->core.l_qseq);

        CHECK_EQ(std::string{seq_gt}, std::string{seq_simd});
        free(seq_simd);
    }

    SUBCASE("Translating intervals from nibbles to bytes") {
        SUBCASE("EVEN") {
            int from = 168;
            int to = 266;
            char *seq = (char *)malloc(to - from + 1);
            nibble2base(bam_get_seq(rec) + from / 2, seq, to - from);

            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{seq});
            free(seq);
        }
        SUBCASE("ODD_EVEN") {
            int from = 167;
            int to = 266;
            char *seq = (char *)malloc(to - from + 2);
            nibble2base(bam_get_seq(rec) + from / 2, seq, to - from + 1);

            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{seq}.substr(1));
            free(seq);
        }
        SUBCASE("ODD_ODD") {
            int from = 167;
            int to = 267;
            char *seq = (char *)malloc(to - from + 2);
            nibble2base(bam_get_seq(rec) + from / 2, seq, to - from + 1);

            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{seq}.substr(1));
            free(seq);
        }
        SUBCASE("EVEN_ODD") {
            int from = 168;
            int to = 267;
            char *seq = (char *)malloc(to - from + 2);
            nibble2base(bam_get_seq(rec) + from / 2, seq, to - from);

            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{seq});
            free(seq);
        }
        SUBCASE("ANY") {
            int seqlen = strlen(seq_gt);
            srand(time(NULL));
            for (int casei = 0; casei < 50; ++casei) {
                int to = rand() % (seqlen - 100) + 100;
                int from = rand() % (to - 100);
                char *seq = (char *)malloc(to - from + 2);

                if (from % 2 == 0) {
                    nibble2base(bam_get_seq(rec) + from / 2, seq, to - from);
                    seq[to - from] = '\0';
                    CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                             std::string{seq});
                } else {
                    nibble2base(bam_get_seq(rec) + from / 2, seq,
                                to - from + 1);
                    seq[to - from + 1] = '\0';
                    CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                             std::string{seq}.substr(1));
                }
                free(seq);
            }
        }
    }
    SUBCASE("Copying Nibble-encoded data and decoding it to bases") {
        int seqlen = strlen(seq_gt);
        srand(time(NULL));
        int N_RANDOM_CASES = 50;
        gtl::vector<bam_nibble_data> copied;

        gtl::vector<size_t> froms;
        gtl::vector<size_t> tos;
        for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
            int to = rand() % (seqlen - 100) + 100;
            int from = rand() % (to - 100);
            copied.emplace_back(bam_get_seq(rec), from, to);
            froms.push_back(from);
            tos.push_back(to);
        }
        char buffer[1024 * 1024];
        for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
            auto &seq = copied[casei];
            size_t from = froms[casei];
            size_t to = tos[casei];
            seq.convert_to_ascii(buffer);
            INFO("from: ", from, " to: ", to, "seq.len: ", seq.len);
            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{buffer}.substr(seq.skip_first ? 1 : 0));
        }
    }
    bam_hdr_destroy(hdr);
    sam_close(fp);
    free(seq_gt);
}

TEST_CASE("BAM Read Parsing With Crumbs") {
    std::string input = "../tests/test_single.bam";
    samFile *fp = sam_open(input.c_str(), "r");
    if (fp == NULL) {
        FAIL("Failed to open BAM file:", input);
    }

    bam_hdr_t *hdr = sam_hdr_read(fp);
    bam1_t *rec = bam_init1();
    if (sam_read1(fp, hdr, rec) <= 0) {
        FAIL("Failed to read BAM record");
    }

    char *seq_gt = (char *)malloc(rec->core.l_qseq + 1);
    for (int i = 0; i < rec->core.l_qseq; ++i) {
        seq_gt[i] = seq_nt16_str[bam_seqi(bam_get_seq(rec), i)];
    }
    seq_gt[rec->core.l_qseq] = '\0';
    SUBCASE(
        "Converting Nibble-encoded data to crumbs and decoding it to bases") {
        int seqlen = strlen(seq_gt);
        srand(time(NULL));
        int N_RANDOM_CASES = 50;

        gtl::vector<CrumbledSubSeq> copied;

        gtl::vector<size_t> froms;
        gtl::vector<size_t> tos;
        for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
            int to = rand() % (seqlen - 100) + 100;
            int from = to - 100 - (rand() % 12);
            copied.emplace_back(bam_get_seq(rec), from, to, CrumbledSubSeq::from_nibbles{});
            froms.push_back(from);
            tos.push_back(to);
        }
        char buffer[1024 * 1024];
        for (int casei = 0; casei < N_RANDOM_CASES; ++casei) {
            auto &seq = copied[casei];
            size_t from = froms[casei];
            size_t to = tos[casei];
            seq.convert_to_ascii(buffer);
            INFO("from: ", from, " to: ", to, "seq.len: ", seq.len);
            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{buffer});
            CrumbledSubSeq fr_chr_cr{seq_gt, from, to,
                                     CrumbledSubSeq::from_ascii{}};
            fr_chr_cr.convert_to_ascii(buffer);
            CHECK_EQ(std::string{seq_gt}.substr(from, to - from),
                     std::string{buffer});
        }
    }
    bam_hdr_destroy(hdr);
    sam_close(fp);
    free(seq_gt);
}
