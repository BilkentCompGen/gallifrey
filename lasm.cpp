#include "ext/CLI11.hpp"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <string>

using std::cerr;
using std::cout;
using std::endl;
using std::string;

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
  CLI::App app("lasm");

  string input;
  string output;

  app.add_option("-i,--input", input, "Input file")->required();
  app.add_option("-o,--output", output, "Output file")->required();

  CLI11_PARSE(app, argc, argv);

  cout << "Input: " << input << endl;
  cout << "Output: " << output << endl;

  uint8_t *data = NULL;

  BAMReader reader(input);

  // read header
  bam1_t *bamdata = reader.init_data();

  while (reader.read(bamdata)) {
    printf("NAME: %s\n", bam_get_qname(bamdata));

    data = bam_get_seq(bamdata);

    printf("SEQ: ");
    for (int i = 0; i < bamdata->core.l_qseq; ++i) {
      printf("%c", seq_nt16_str[bam_seqi(data, i)]);
    }

    break;
  }
  if (bamdata) {
    bam_destroy1(bamdata);
  }
  return 0;
}
