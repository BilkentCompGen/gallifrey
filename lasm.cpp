#include <CLI/CLI.hpp>
#include <cstdint>
#define LCPTOOLS_IMPL
#include "ext/lcptools/lcptools_ho.h"
#include "ext/bloom_filter.hpp"
#include "lcmer.h"
#include <concurrentqueue.h>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

struct core_view {
  core *from;
  core *to;
};
// TODO SYNCMER
// DCT
// CORE
// DEEPEN
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

uint64_t parse1_bam(const uint8_t *data, int64_t len, struct core *cores) {

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
void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len) {
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
struct dbg_value{
  size_t count;
  size_t index;
  gtl::vector<uint32_t> next;
};
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
  int threads = 8;
  app.add_option("-i,--input", input, "Input bam file")->required();
  app.add_option("-o,--output", output, "Output gfa file")->required();
  app.add_option("-l,--level", lcp_level, "LCP level to compute")->required();
  app.add_option("-t,--threads", threads, "Number of threads")
      ->capture_default_str();
  CLI11_PARSE(app, argc, argv);
  cout << "Input: " << input << endl;
  cout << "Output: " << output << endl;


  // bloom_parameters bparam;
  // bparam.projected_element_count = 1'000'000'000;
  // bparam.false_positive_probability = 0.0001; 
  // bparam.random_seed = 0xA5A5A5A5;
  // bparam.compute_optimal_parameters();
  // bloom_filter bfil(bparam);


  uint8_t *data = NULL;

  std::atomic<bool> done{false};
  std::atomic<int> count;
  BAMReader reader(input);
  bam1_t *bamdata = reader.init_data();
  using worktype = std::tuple<uint64_t, uint8_t *>;
  using mutex_type = std::mutex;
  moodycamel::ConcurrentQueue<worktype> work_queue(10, 1, 0);
  using dbg_value_type = dbg_value;
  using dbg_map_t = gtl::parallel_flat_hash_map<
      lcmer<L>, dbg_value_type, std::hash<lcmer<L>>,
      gtl::priv::hash_default_eq<lcmer<L>>,
      std::allocator<std::pair<const lcmer<L>, dbg_value_type>>, 4, mutex_type>;
  std::vector<std::thread> consumers;

  dbg_map_t dbg;
  for (int tid = 0; tid < threads; ++tid) {

    consumers.emplace_back([&]() {
      uint64_t lcmer_index = tid;
      worktype work;
      while (!done.load() || work_queue.size_approx() != 0) {
        if (!work_queue.try_dequeue(work)) {
          continue;
        }
        lcp_read lr = parse_cores_from_bam_record(std::get<1>(work),
                                                  std::get<0>(work), lcp_level);
        // cout << bamdata->core.l_qseq << " to " << lr.cores.size() << "\n";
        if (lr.cores.size() < L + 2) {
          free(std::get<1>(work));
          continue;
        }
        lcmer<L> l{lr.cores.begin()};
        for (size_t i = L; i < lr.cores.size(); ++i) {
          // if(!bfil.contains(l)){
          //   bfil.insert(l);
          //   update_lcmer(l, lr.cores[i]);
          //   continue;
          // }
          dbg.lazy_emplace_l(
              l,
              [&](dbg_map_t::value_type &val) {
                val.second.count++;
                auto &vec = val.second.next;
                if (std::find(vec.begin(), vec.end(), lr.cores[i]) !=
                    vec.end()) { // TODO check performance.
                  vec.push_back(lr.cores[i]);
                }
              },
              [&](const dbg_map_t::constructor &ctor) {
                ctor(l, dbg_value_type{1UL, lcmer_index,
                                        gtl::vector<uint32_t>{lr.cores[i]}});
                lcmer_index += 128;
              });
          update_lcmer(l, lr.cores[i]);
          count.fetch_add(1);
        }
        // if(!bfil.contains(l)){
        //     bfil.insert(l);
        //     continue;
        // }
        dbg.lazy_emplace_l(
            l, [&](dbg_map_t::value_type &val) { 
              val.second.count++; 
            },
            [&](const dbg_map_t::constructor &ctor) {
              ctor(l, dbg_value_type{1UL, lcmer_index, gtl::vector<uint32_t>{}});
              lcmer_index += 128;
            });
        free(std::get<1>(work));
      }
    });
  }

  while (reader.read(bamdata)) {
    // printf("NAME: %s\n", bam_get_qname(bamdata));
    data = bam_get_seq(bamdata);
    uint8_t *data_copy = (uint8_t *)malloc(bamdata->core.l_qseq / 2 +
                                           1); // TODO Fix the extra memory
    memcpy(data_copy, data, bamdata->core.l_qseq / 2 + 1);
    while (!work_queue.try_enqueue(
        std::make_tuple(bamdata->core.l_qseq, data_copy)))
      ;
  }
  done.store(true);
  for (std::thread &t : consumers) {
    t.join();
  }
  if (bamdata) {
    bam_destroy1(bamdata);
  }
  return 0;
}
