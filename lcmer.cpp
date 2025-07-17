#include "lcmer.h"

#include "syncmer.h"

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
                             struct core *cores) {
    int core_index = 0;
    KmerViewOnNibbles<SYNCMER_K, uint32_t> kmer_view;
    int64_t it1 = kmer_view.fill(data);
    // find lcp cores
    for (; it1 + SYNCMER_K < len; it1++) {
        // printf("%d\t%d\n", kmer_view.get(), table[kmer_view.get()]);
        if (table[kmer_view.get()]) {
            init_core_kmer<uint32_t>(&(cores[core_index]), kmer_view.get(), it1,
                                     it1 + SYNCMER_K);
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

void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len) {
    lps_ptr->level = 1;
    lps_ptr->size = 0;
    lps_ptr->cores =
        (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
    lps_ptr->size = parse1_bam_syncmers(data, len, lps_ptr->cores);
}

lcp_read parse_cores_from_bam_record(const uint8_t *data, size_t len,
                                     int lcp_level) {
    // uint64_t estimated_core_size = (int)(len / std::pow(1.5, lcp_level));

    lps str;
    init_lps_bam(&str, data, len);
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
