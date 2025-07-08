#include <CLI/CLI.hpp>
#include <cstdint>
#define LCPTOOLS_IMPL
#include "ext/bloom_filter.hpp"
#include "ext/lcptools/lcptools_ho.h"
#include "lcmer.h"
#include "syncmer.h"
#include <concurrentqueue.h>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <string>

#include "config.h"

using std::cout;
using std::endl;
using std::string;

struct core_view {
  core *from;
  core *to;
};


//https://github.com/ruanjue/wtdbg2/blob/b77c5657c8095412317e4a20fe3668f5bde6b1ac/dna.h#L108
template<class T>
T reverse_complement_kmer(T kmer, int K){
  uint32_t shift_count = (sizeof(T)*8) - 2 * K;
  kmer = ~kmer;
  kmer = ((kmer & 0x3333333333333333LLU)<< 2) | ((kmer & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
  kmer = ((kmer & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((kmer & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
  if constexpr (sizeof(T)==8){
    kmer= __builtin_bswap64(kmer);
  }
  else if constexpr (sizeof(T)==4){
    kmer= __builtin_bswap32(kmer);
  }
  else if constexpr (sizeof(T)==2){
    kmer= __builtin_bswap16(kmer);
  }
  kmer = kmer >> shift_count;
  return kmer;
}


template <int _K, typename INT_TYPE = uint32_t> struct KmerViewOnNibbles {
  using int_type = INT_TYPE;
  static constexpr int k = _K;
  int_type value;
  static constexpr int_type KMASK = ((int_type)1 << k) - 1;
  uint8_t nibble_mask = 15;
  KmerViewOnNibbles() : value{0} {}

  int fill(const uint8_t *data_address) {
    int i;
    for (i = 0; i < k; i++) {
      move_get(data_address, i);
    }
    return i;
  }
  int_type move_get(const uint8_t *data_address, int index) {
    value <<= 2;
    if constexpr (KMASK != -1)
      value &= KMASK;
    value |= bam_seqi(data_address, index);
    return value;
  }
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

template <typename _T>
void init_core_kmer(struct core *cr, _T Kval, uint64_t si, uint64_t ei) {
  cr->start = si;
  cr->end = ei;
  cr->label = 0;
  cr->label |= Kval;
  cr->bit_rep = 0x8000000000000000 | cr->label;
  cr->bit_size = 2 * SYNCMER_K;
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
syncmer::syncmer_precomputed<syncmer::SYNCMER_RULE<SYNCMER_K, SYNCMER_S, SYNCMER_T>> table;
uint64_t parse1_syncmers(const uint8_t *data, int64_t len, struct core *cores) {
  int core_index = 0;
  KmerViewOnNibbles<SYNCMER_K, uint32_t> kmer_view;
  int64_t it1 = kmer_view.fill(data);
  // find lcp cores
  for (; it1 + SYNCMER_K < len; it1++) {

    if (table[kmer_view.value]) {
      init_core_kmer_xor<uint32_t>(&(cores[core_index]), kmer_view.value, it1,
                                  it1 + SYNCMER_K);
      core_index++;
    }
    kmer_view.move_get(data,it1);
    // cout << "None\n";
  }
  // cout << "CORESSS" << it1 << "\t" << it2 << "\n";
  return core_index;
}

uint64_t parse1_bam(const uint8_t *data, int64_t len, struct core *cores) {

  int64_t it1 = 0;
  int64_t it2 = len;

  int core_index = 0;
  // int64_t last_invalid_char_index = -1;

  // find lcp cores
  for (; it1 + 2 < len; it1++) {

    if (bam_seqi(data, it1) == bam_seqi(data, it1 + 1)) {
      continue;
    }

    // // check for RINT core
    // if (bam_seqi(data, (it1 + 1)) == bam_seqi(data, (it1 + 2))) {

    //   // count middle characters
    //   uint32_t middle_count = 1;
    //   int64_t temp = it1 + 2;
    //   while (temp < len && bam_seqi(data, (temp - 1)) == bam_seqi(data,
    //   temp)) {
    //     temp++;
    //     middle_count++;
    //   }
    //   if (temp != len) {
    //     // check if there is any SSEQ cores left behind
    //     if (it2 < it1 && last_invalid_char_index < it2 - 1) {
    //       init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
    //                     it1 + 1);
    //       core_index++;
    //     }

    //     // create RINT core
    //     it2 = it1 + 2 + middle_count;
    //     init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
    //     core_index++;

    //     continue;
    //   }
    // }

    // if (bam_seqi(data, it1) > bam_seqi(data, (it1 + 1)) &&
    //     bam_seqi(data, (it1 + 1)) < bam_seqi(data, (it1 + 2))) {

    //   // check if there is any SSEQ cores left behind
    //   if (it2 < it1 && last_invalid_char_index < it2 - 1) {
    //     init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
    //                   it1 + 1);
    //     core_index++;
    //   }

    //   // create LMIN core
    //   it2 = it1 + 3;
    //   init_core1bam(&(cores[core_index]), data, it2 - it1, it1, it2);
    //   core_index++;

    //   continue;
    // }

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
      //   init_core1bam(&(cores[core_index]), data, it1 - it2 + 2, it2 - 1,
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

/// @brief Initializes an lps struct with a given sequence in the bam compressed
/// format
void init_lps_bam(lps *lps_ptr, const uint8_t *data, size_t len) {
  lps_ptr->level = 1;
  lps_ptr->size = 0;
  lps_ptr->cores =
      (struct core *)malloc((len / CONSTANT_FACTOR) * sizeof(struct core));
  lps_ptr->size = parse1_syncmers(data, len, lps_ptr->cores);
}

struct lcp_read {
  gtl::vector<uint32_t> cores;
  gtl::vector<std::pair<uint32_t, uint32_t>> pos;
  void reserve(size_t size) {
    cores.reserve(size);
    pos.reserve(size);
  }
};
struct dbg_value {
  size_t count;
  size_t index;
  gtl::parallel_flat_hash_set<uint32_t, gtl::priv::hash_default_hash<uint32_t>,
                              gtl::priv::hash_default_eq<uint32_t>,
                              std::allocator<uint32_t>, 2, std::mutex>
      next;
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
    // fprintf(stderr, "%u\n", str.cores[i].bit_size);
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


int main(int argc, char **argv) {

  // syncmer::syncmer_precomputed<syncmer::OpenSyncmerRule<8, 3, 4, uint16_t>>
  //     table;
  // // syncmer::syncmer_precomputed<syncmer::DensitySyncmerRule<8, 1, 100,
  // // uint16_t>> table;
  // for (size_t i = 0; i < table.size; ++i) {
  //   printf("%lx \t %s\n", i, table[i] ? "syncmer" : "not syncmer");
  // }
  // return 0;

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

  bloom_parameters bparam;
  bparam.projected_element_count = 1'000'000'000;
  bparam.false_positive_probability = 0.0001;
  bparam.random_seed = 0xA5A5A5A5;
  bparam.compute_optimal_parameters();
  bloom_filter bfil(bparam);

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
      lcmer<LCMER_L>, dbg_value_type, std::hash<lcmer<LCMER_L>>,
      gtl::priv::hash_default_eq<lcmer<LCMER_L>>,
      std::allocator<std::pair<const lcmer<LCMER_L>, dbg_value_type>>, 4, mutex_type>;
  using counter = gtl::parallel_flat_hash_map<
      size_t,size_t, std::hash<size_t>,
      gtl::priv::hash_default_eq<size_t>,
      std::allocator<std::pair<size_t, size_t>>, 4, mutex_type>;

    
  std::vector<std::thread> consumers;
  counter core_counter;
  dbg_map_t dbg;
  for (int tid = 0; tid < threads; ++tid) {

    consumers.emplace_back(
        [&](int t) {
          uint64_t lcmer_index = 0;
          worktype work;
          while (!done.load() || work_queue.size_approx() > 0) {
            if (!work_queue.try_dequeue(work)) {
              continue;
            }
            lcp_read lr = parse_cores_from_bam_record(
                std::get<1>(work), std::get<0>(work), lcp_level);
            // cout << bamdata->core.l_qseq << " to " << lr.cores.size() <<
            // "\n";
            core_counter[lr.cores.size()]++;
            if (lr.cores.size() < LCMER_L + 2) {
              free(std::get<1>(work));
              continue;
            }
            lcmer<LCMER_L> l{lr.cores.begin()};
            for (size_t i = LCMER_L; i < lr.cores.size(); ++i) {
              // if(!bfil.contains(l)){
              //   bfil.insert(l);
              //   update_lcmer(l, lr.cores[i]);
              //   continue;
              // }
              dbg.lazy_emplace_l(
                  l,
                  [&](dbg_map_t::value_type &val) {
                    val.second.count++;
                    val.second.next.insert(lr.cores[i]);
                  },
                  [&](const dbg_map_t::constructor &ctor) {
                    ctor(l,
                         dbg_value_type{1UL, lcmer_index | t, {lr.cores[i]}});
                    lcmer_index += 128;
                  });
              l = update_lcmer(l, lr.cores[i]);
              count.fetch_add(1);
            }
            // if(!bfil.contains(l)){
            //     bfil.insert(l);
            //     continue;
            // }
            dbg.lazy_emplace_l(
                l, [](dbg_map_t::value_type &val) { val.second.count++; },
                [&](const dbg_map_t::constructor &ctor) {
                  ctor(l, dbg_value_type{1UL, lcmer_index | t, {}});
                  lcmer_index += 128;
                });
            free(std::get<1>(work));
          }
        },
        tid);
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
  for(auto &p: core_counter){
    fprintf(stderr, "CC\t%lu\t%lu\n", p.first, p.second);
  }
  if (bamdata) {
    bam_destroy1(bamdata);
  }

  FILE *out = fopen(output.c_str(), "w");
  fprintf(stderr, "Something printed\n");
  for (const auto &p : dbg) {
    fprintf(out, "S\t%lu\t*\tLN:i:%d\tKC:i:%lu\n", p.second.index, LCMER_L,
            p.second.count);
    // The code snippet provided seems to be a mix of 
    // Explain: 
  }
  for (const auto &p : dbg) {
    lcmer ln = shift_lcmer(p.first);
    for (uint32_t next : p.second.next) {
      lcmer_set0(ln, next);
      fprintf(out, "L\t%lu\t+\t%lu\t+\t%dM\n", p.second.index, dbg[ln].index,
              LCMER_L - 1);
    }
  }

  

  fclose(out);
  return 0;
}
