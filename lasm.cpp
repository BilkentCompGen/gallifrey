#include <CLI/CLI.hpp>
#define LCPTOOLS_IMPL
#include "ext/lcptools/lcptools_ho.h"
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <string>

#include "lcmer.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

struct core_view {
  core *from;
  core *to;
};

void init_core1bam(struct core *cr, const uint8_t *data, uint64_t distance,
                   uint64_t si, uint64_t ei) {
  cr->start = si;
  cr->end = ei;
  cr->label = 0;
  cr->label |= ((distance - 2) << 6);
  cr->label |= (bam_seqi(data, si) << 4);
  cr->label |= (bam_seqi(data, si + distance - 2) << 2);
  cr->label |= (bam_seqi(data, si + distance - 1));
  cr->bit_rep = 0x8000000000000000 | cr->label;
  cr->bit_size = 2 * distance;
}

uint64_t parse1_bam(uint8_t *data, int64_t len, struct core *cores) {

  int64_t it1 = 0;
  int64_t it2 = len;

  int core_index = 0;
  int64_t last_invalid_char_index = -1;

  // find lcp cores
  for (; it1 + 2 < len; it1++) {

    if (bam_seqi(data, it1) == bam_seqi(data, it1 + 1)) {
      continue;
    }

    // check for RINT core
    if (bam_seqi(data, (it1 + 1)) == bam_seqi(data, (it1 + 2))) {

      // count middle characters
      uint32_t middle_count = 1;
      int64_t temp = it1 + 2;
      while (temp < len && bam_seqi(data, (temp - 1)) == bam_seqi(data, temp)) {
        temp++;
        middle_count++;
      }
      if (temp != len) {
        // check if there is any SSEQ cores left behind
        if (it2 < it1 && last_invalid_char_index < it2 - 1) {
          init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
                        it1 + 1);
          core_index++;
        }

        // create RINT core
        it2 = it1 + 2 + middle_count;
        init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
        core_index++;

        continue;
      }
    }

    if (bam_seqi(data, it1) > bam_seqi(data, (it1 + 1)) &&
        bam_seqi(data, (it1 + 1)) < bam_seqi(data, (it1 + 2))) {

      // check if there is any SSEQ cores left behind
      if (it2 < it1 && last_invalid_char_index < it2 - 1) {
        init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
                      it1 + 1);
        core_index++;
      }

      // create LMIN core
      it2 = it1 + 3;
      init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
      core_index++;

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

      // check if there is any SSEQ cores left behind
      if (it2 < it1 && last_invalid_char_index < it2 - 1) {
        init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
                      it1 + 1);
        core_index++;
      }

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

/// @brief Initializes an lps struct with a given sequence in the bam compressed
/// format
void init_lps_bam(lps *lps_ptr, uint8_t *data, size_t len) {
  lps_ptr->level = 1;
  lps_ptr->size = 0;
  lps_ptr->cores =
      (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
  lps_ptr->size = parse1_bam(data, len, lps_ptr->cores);
}

struct lcp_read {
  gtl::vector<uint32_t> cores;
  gtl::vector<std::pair<uint32_t, uint32_t>> pos;
  void reserve(size_t size) {
    cores.reserve(size);
    pos.reserve(size);
  }
};

lcp_read parse_cores_from_bam_record(uint8_t *data, size_t len, int lcp_level) {
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
  return lr;
}

struct BAMReader {
  string filename;
  samFile *fp;
  sam_hdr_t *hdr;

  BAMReader(const string &filename) : fp{sam_open(filename.c_str(), "rb")} {
    this->filename = filename;
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
  bam1_t *init_data() {
    bam1_t *data = bam_init1();
    if (data == NULL) {
      throw "Not enough memory";
    }
    return data;
  }
};
#define L 8
int main(int argc, char **argv) {
  CLI::App app("lasm");

  string input;
  string output;
  int lcp_level;
  app.add_option("-i,--input", input, "Input file")->required();
  app.add_option("-o,--output", output, "Output file")->required();
  app.add_option("-l,--level", lcp_level, "LCP level to compute")->required();
  CLI11_PARSE(app, argc, argv);

  cout << "Input: " << input << endl;
  cout << "Output: " << output << endl;

  uint8_t *data = NULL;

  BAMReader reader(input);
  bam1_t *bamdata = reader.init_data();

  while (reader.read(bamdata)) {
    // printf("NAME: %s\n", bam_get_qname(bamdata));

    data = bam_get_seq(bamdata);

    lcp_read lr =
        parse_cores_from_bam_record(data, bamdata->core.l_qseq, lcp_level);
    // cout << bamdata->core.l_qseq << " to " << lr.cores.size() << "\n";
    if(lr.cores.size() < L + 2){
      continue;
    }
    lcmer<L> l{lr.cores.begin()};

    
    for(size_t i = L; i < lr.cores.size(); ++i){
      update_lcmer(l,lr.cores[i]);
    }

    //   for(std::pair<uint32_t, uint32_t> &p : lr.pos){
    //     cout << p.first << " " << p.second << "\n";
    // }
    // printf("SEQ: ");
    // for (int i = 0; i < bamdata->core.l_qseq; ++i) {
    //   printf("%c", seq_nt16_str[bam_seqi(data, i)]);
    // }
  }
  if (bamdata) {
    bam_destroy1(bamdata);
  }
  return 0;
}
