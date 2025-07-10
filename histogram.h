
#include <stdio.h>
#include <stdlib.h>

#include <cstdio>
#include <vector>


struct histo {
    std::string name;
    std::vector<unsigned> bins;
    unsigned bin_width;
    unsigned scale_start;
    unsigned scale_div;
    unsigned max_c_per_line = 60;
    unsigned max_lines = 16;
    std::vector<char> characters{{'#', '=', '_', '.'}};
    histo(std::string name, unsigned bin_width, unsigned ss, unsigned sd)
        : name{name}, bin_width{bin_width}, scale_start{ss}, scale_div{sd} {}
    int print(FILE *f) {
        int idx = 0;
        fprintf(f, "%s\n", name.c_str());
        for (unsigned c : bins) {
            unsigned scale = scale_start;
            int need_more;
            fprintf(f, "%u%n", bin_width * idx, &need_more);
            for (int i = 0; i < 12 - need_more; ++i) {
                fputc(' ', f);
            };
            fputc(':', f);
            while (max_c_per_line < c / scale_start) {
                scale_start *= 2;
            }
            for (unsigned i = 0; i < c / scale; ++i) fputc(characters[0], f);

            c = c % scale_div;
            unsigned i;
            for (i = 1; i < characters.size(); ++i) {
                scale /= scale_div;
                if (c > scale) {
                    fputc(characters[i], f);
                    break;
                }
            }
            fprintf(f, " %u%n", bins[idx], &need_more);
            for (; i < max_c_per_line - 1 - need_more; ++i) {
                fputc(' ', f);
            }
            fputs("\n", f);
            idx += 1;
        }
        fputs(
            "------                                                            "
            "\n",
            f);
        return idx + 2;
    }
    void resize_bins(unsigned new_width) {
        std::vector<unsigned> new_bins;
        new_bins.resize(bins.size());
        for (unsigned i = 0; i < bins.size(); ++i) {
            unsigned tgt = i * bin_width / new_width;
            new_bins[tgt] += bins[i];
        }
        bins = new_bins;
        bin_width = new_width;
    }
    void reroll(FILE *f, int idx) {
        while (idx > 0) {
            fputs("\033[F", f);
            --idx;
        }
    }
    void push(unsigned value) {
        auto bin_tgt = value / bin_width;
        if (bin_tgt > max_lines) {
            while (bin_tgt > max_lines) {
                resize_bins(bin_width * 1.25);
                bin_tgt = value / bin_width;
            }
            push(value);
            return;
        }
        if (bin_tgt >= bins.size()) {
            bins.resize(1 + bin_tgt);
        }
        bins[bin_tgt]++;
    }
};
// int main(){
//     srand(time(NULL));
//     histo h{10,20,4};
//     for(int i =0; i < 1000; ++i){
//         h.push(rand() % 50);
//         std::this_thread::sleep_for(5ms);
//         h.reroll(stdout, h.print(stdout));
//         fflush(stdout);
//     }
//     h.print(stdout);

// }