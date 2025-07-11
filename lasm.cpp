#include <CLI/CLI.hpp>

#include <string>
#include "asm.h"
#include "stats.h"
#include <memory>
#define LCPTOOLS_IMPL
#include "ext/lcptools/lcptools_ho.h"

using std::string;

struct bloom_parameters;
int asm_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("hifi", "run hifi assembly");
    auto &app = *app_;
    std::shared_ptr<string> input{new string};
    std::shared_ptr<string> output{new string};
    std::shared_ptr<int> lcp_level{new int};
    std::shared_ptr<int> threads{new int{8}};
    app.add_option("input", *input, "Input bam file")->required();
    app.add_option("output", *output, "Output gfa file")->required();
    app.add_option("-l,--level", *lcp_level, "LCP level to compute")
        ->required();
    app.add_option("-t,--threads", *threads, "Number of threads")
        ->capture_default_str();

    std::shared_ptr<bloom_parameters> bparam{new bloom_parameters};
    bparam->projected_element_count = 1'000'000;
    bparam->false_positive_probability = 0.0001;
    bparam->random_seed = 0xA5A5A5A5;
    auto bfg =
        app.add_option_group("bloomfilter", "Parameters to setup bloom filter");
    bfg->add_option("--projected-element-count",
                    bparam->projected_element_count, "Projected element count")
        ->capture_default_str();
    bfg->add_option("--false-positive-probability",
                    bparam->false_positive_probability,
                    "False positive probability")
        ->capture_default_str();
    bfg->add_option("--random-seed", bparam->random_seed, "Random seed")
        ->capture_default_str();

    app.callback([=]() {
        // cout << "Input: " << *input << endl;
        // cout << "Output: " << *output << endl;
        return run_assembly(*input, *output,  *threads, *lcp_level, *bparam);

    });
    return 0;
}

void stats_program(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("bam-stats", "Print stats on input bam");
    auto &app = *app_;
    std::shared_ptr<string> input{new string};
    std::shared_ptr<uint64_t> update_frequency{new uint64_t{100}};
    std::shared_ptr<int> lcmer_level{new int{5}};
    struct parameters {
        string input_file;
        int update_frequency;
        int lcmer_level;
    };
    auto params = std::make_shared<parameters>("", 100, 5);
    app.add_option("input", params->input_file, "Input bam file")->required();
    app.add_option("--update-frequency", params->update_frequency, "Update frequency")
        ->capture_default_str();
    app.add_option("-l,--level", params->lcmer_level, "Depth of lcmers to be computed")
        ->capture_default_str();

    app.callback([=]() {
        return run_stats(params->input_file, params->update_frequency, params->lcmer_level);
    });
}

void test_htslib_file_movements(CLI::App &parent) {
    auto *app_ = parent.add_subcommand("test-hts",
                                       "Testing HTSLib Multi File movements");
    auto &app = *app_;
    struct parameters {
        string input_file;
        int threads;
        int block_size;
    };
    // Allocate parameters as heap and set default values
    std::shared_ptr<parameters> params{new parameters{"", 8, 4}};
    app.add_option("input", params->input_file, "Input file")->required();
    app.add_option("-t,--threads", params->threads, "Number of threads")
        ->capture_default_str();
    app.add_option("-b,--block-size", params->block_size, "Block size")
        ->capture_default_str();

    app.callback([=]() {

    });
}
int main(int argc, char **argv) {
    CLI::App app("lasm");
    app.require_subcommand(1);
    app.set_help_all_flag("--help-all", "Expand all help");
    asm_program(app);
    stats_program(app);
    test_htslib_file_movements(app);
    CLI11_PARSE(app, argc, argv);
    return 0;
}