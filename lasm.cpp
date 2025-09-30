#include <CLI/CLI.hpp>
#include <string>


#define LCPTOOLS_IMPL
#include "ext/lcptools/lcptools_ho.h" 

#include "asm.h"
#include "stats.h"
#include "sthist.h"
#include "bamdoctor.h"


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
    bparam->false_positive_rate = 0.01;
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

int sam_read_testing(CLI::App &parent) {
    struct parameters {
        string input;
        int lcp_level;
        int threads;
        SamReadTestingMode mode;
        string job;
    };
    auto *app_ = parent.add_subcommand("bamdoctor", "read bam file fast");
    auto &app = *app_;
    parameters *p = new parameters{"", 5, 8, SamReadTestingMode::BASIC, "read_count"};  
    app.add_option("input", p->input, "Input bam file")->required();
    app.add_option("-l,--level", p->lcp_level, "LCP level to compute")
        ->capture_default_str();
    app.add_option("-t,--threads", p->threads, "Number of threads")
        ->capture_default_str();

    std::map<std::string, SamReadTestingMode> SamReadTestingMode_map = {
        {"basic", SamReadTestingMode::BASIC},
        {"pool", SamReadTestingMode::POOL},
    };

    std::map<std::string, std::string> job_map = {
        {"read_count", "read_count"},
        {"core_count", "core_count"},
    };
    app.add_option("--job", p->job, "Job to run")->required()->transform(CLI::CheckedTransformer(job_map, CLI::ignore_case));
    app.add_option("--mode", p->mode, "Mode to run")->required()->transform(CLI::CheckedTransformer(SamReadTestingMode_map, CLI::ignore_case));
    app.callback([=]() {
        int ret = run_sam_read_testing(p->input, p->threads, p->lcp_level, p->mode, p->job);
        delete p;
        return ret;
    });
    return 0;
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


#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

int tests(int argc, char** argv) {
    doctest::Context context;

    // !!! THIS IS JUST AN EXAMPLE SHOWING HOW DEFAULTS/OVERRIDES ARE SET !!!

    // defaults
    context.addFilter("test-case-exclude", "*math*"); // exclude test cases with "math" in their name
    context.setOption("abort-after", 5);              // stop test execution after 5 failed assertions
    context.setOption("order-by", "name");            // sort the test cases by their name

    context.applyCommandLine(argc, argv);

    // overrides
    context.setOption("no-breaks", true);             // don't break in the debugger when assertions fail

    int res = context.run(); // run

    if(context.shouldExit()) // important - query flags (and --exit) rely on the user doing this
        return res;          // propagate the result of the tests
    
    int client_stuff_return_code = 0;
    // your program - if the testing framework is integrated in your production code
    
    return res + client_stuff_return_code; // the result from doctest is propagated here as well
}

int main(int argc, char **argv) {
    if(string{argv[1]} == "tests"){
        return tests(argc-1, argv+1);
    }
    CLI::App app("lasm");
    app.require_subcommand(1);
    app.set_help_all_flag("--help-all", "Expand all help");
    asm_program(app);
    stats_program(app);
    stdin_hist_prog(app);
    sam_read_testing(app);
    CLI11_PARSE(app, argc, argv);
    return 0;
}
