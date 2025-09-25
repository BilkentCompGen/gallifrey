#include <CLI/CLI.hpp>
#include <string>
#include "bamdoctor.h"

using std::string;

int main(int argc, char **argv) {
    CLI::App app("bamdoctor - Fast BAM file reader and analyzer");
    
    struct parameters {
        string input;
        int lcp_level = 5;
        int threads = 8;
        SamReadTestingMode mode = SamReadTestingMode::BASIC;
        string job = "read_count";
    };
    
    parameters params;
    
    app.add_option("input", params.input, "Input BAM file")->required();
    app.add_option("-l,--level", params.lcp_level, "LCP level to compute")
        ->capture_default_str();
    app.add_option("-t,--threads", params.threads, "Number of threads")
        ->capture_default_str();
    app.add_option("--job", params.job, "Job to run")->capture_default_str();

    std::map<std::string, SamReadTestingMode> mode_map = {
        {"basic", SamReadTestingMode::BASIC},
        {"pool", SamReadTestingMode::POOL},
    };

    app.add_option("--mode", params.mode, "Reading mode (basic or pool)")
        ->transform(CLI::CheckedTransformer(mode_map, CLI::ignore_case))
        ->capture_default_str();

    CLI11_PARSE(app, argc, argv);
    
    return run_sam_read_testing(params.input, params.threads, params.lcp_level, params.mode, params.job);
}
