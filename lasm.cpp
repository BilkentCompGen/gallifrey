#include <CLI/CLI.hpp>
#include <string>

#include "asm.h"
#include "stats.h"
#include "sthist.h"
#define LCPTOOLS_IMPL
#include "ext/lcptools/lcptools_ho.h" 

using std::string;

struct bloom_parameters;
int asm_program(CLI::App &parent) {
    struct parameters {
        string input;
        string output;
        int lcp_level;
        int threads;
    };
    auto *app_ = parent.add_subcommand("hifi", "run hifi assembly");
    auto &app = *app_;
    parameters *p = new parameters{"", "", 5, 8};
    app.add_option("input", p->input, "Input bam file")->required();
    app.add_option("output", p->output, "Output gfa file")->required();
    app.add_option("-l,--level", p->lcp_level, "LCP level to compute")
        ->required();
    app.add_option("-t,--threads", p->threads, "Number of threads")
        ->capture_default_str();

    bloom_param *bparam = new bloom_param;
    bparam->projected_count = 100'000'000;
    bparam->false_positive_rate = 0.0001;
    bparam->seed = 0xA5A5A5A5;
    auto bfg =
        app.add_option_group("bloomfilter", "Parameters to setup bloom filter");
    bfg->add_option("--projected-element-count", bparam->projected_count,
                    "Projected element count")
        ->capture_default_str();
    bfg->add_option("--false-positive-probability", bparam->false_positive_rate,
                    "False positive probability")
        ->capture_default_str();
    bfg->add_option("--random-seed", bparam->seed, "Random seed")
        ->capture_default_str();

    app.callback([=]() {
        int ret = run_assembly(p->input, p->output, p->threads, p->lcp_level,
                               *bparam);
        delete p;
        delete bparam;
        return ret;
    });
    return 0;
}

void stats_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("bam-stats", "Print stats on input bam");
    auto &app = *app_;
    struct parameters {
        string input_file;
        int update_frequency;
        int lcmer_level;
    };
    parameters *params = new parameters{"", 100, 5};
    app.add_option("input", params->input_file, "Input bam file")->required();
    app.add_option("--update-frequency", params->update_frequency,
                   "Update frequency")
        ->capture_default_str();
    app.add_option("-l,--level", params->lcmer_level,
                   "Depth of lcmers to be computed")
        ->capture_default_str();

    app.callback([=]() {
        int ret = run_stats(params->input_file, params->update_frequency,
                            params->lcmer_level);
        delete params;
        return ret;
    });
}

void stdin_hist_prog(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("mephisto",
                                       "Testing HTSLib Multi File movements");
    auto &app = *app_;
    struct parameters {
        int dims;
        int freq;
    };
    // Allocate parameters as heap and set default values
    parameters *params = new parameters{1, 10};
    app.add_option("--dims", params->dims, "screen dimensions")
        ->capture_default_str();
    app.add_option("--freq", params->freq, "Frequency of histogram updates")
        ->capture_default_str();

    app.callback([=]() {
        int ret = run_hist_on_stdin(params->dims, params->freq);
        delete params;
        return ret;
    });
}
int main(int argc, char **argv) {
    CLI::App app("lasm");
    app.require_subcommand(1);
    app.set_help_all_flag("--help-all", "Expand all help");
    asm_program(app);
    stats_program(app);
    stdin_hist_prog(app);
    CLI11_PARSE(app, argc, argv);
    return 0;
}
