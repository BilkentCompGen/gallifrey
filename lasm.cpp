#include <CLI/CLI.hpp>
#include <map>
#include <string>
// #include <vector>

#include "CLI/CLI.hpp"

#define LCPTOOLS_IMPL
#include "asm.h"
#include "bamdoctor.h"
#include "dump_cores.h"
#include "ext/lcptools/lcptools_ho.h"
using std::string;

// Each subcommand is a command object that owns its parsed parameters and
// registers itself with the parent CLI::App on construction. The objects
// live for all of main (enclosing CLI11_PARSE), so bound options and the
// `this`-capturing callbacks stay valid. Copy/move is deleted.

class AsmProgram {
    asm_config cfg;
    std::string file_reading_mode{"sequential"};
    int olc_lcmer_l_ignored = 0;  // --olc-lcmer-l: accepted, no effect
    std::string multipass_text;

   public:
    explicit AsmProgram(CLI::App &parent) {
        auto &app = *parent.add_subcommand("hifi", "run hifi assembly");
        app.add_option("input", cfg.inputs,
                       "Input read files: BAM, CRAM, SAM, or FASTA/FASTQ "
                       "(plain or gzipped). Format is detected from content; "
                       "formats may be mixed.")
            ->required();
        app.add_option("-o,--output", cfg.output, "Output gfa file")
            ->required();
        auto *level_opt =
            app.add_option("-l,--level", cfg.lcp_level, "LCP level to compute")
                ->check(CLI::Range(1, 64));
        auto *lcmer_opt =
            app.add_option("--lcmer-l,--dbg-lcmer-l", cfg.lcmer_l,
                           "LCMER window length, used for DBG construction "
                           "and overlap detection (--dbg-lcmer-l is the old "
                           "name)")
                ->check(CLI::Range(1, 64))
                ->capture_default_str();
        auto *multipass_opt =
            app.add_option(
                   "--multipass", multipass_text,
                   "Multi-pass assembly: LEVEL:L[,LEVEL:L...] gives the LCP "
                   "level and lcmer window length of each pass, e.g. "
                   "4:8,5:6,5:10 (replaces -l and --lcmer-l). Each pass "
                   "assembles the reads together with the previous pass's "
                   "contigs; only the last pass writes the output. Order "
                   "the passes from small to large windows.")
                ->excludes(level_opt)
                ->excludes(lcmer_opt);
        app.add_option("--olc-lcmer-l", olc_lcmer_l_ignored,
                       "Ignored; overlap detection uses --lcmer-l")
            ->group("");
        app.add_option("-t,--threads", cfg.threads.work,
                       "Number of worker threads")
            ->check(CLI::Range(1, MAXIMUM_WORK_THREADS))
            ->capture_default_str();
        app.add_option("-@,--io-threads", cfg.threads.io,
                       "Number of IO threads (Used by htslib for decompressing)")
            ->check(CLI::Range(0, 64))
            ->capture_default_str();
        app.add_option("--random-seed", cfg.seed,
                       "Seed of the syncmer s-mer hash (32-bit)")
            ->check(CLI::Range(uint64_t{0}, uint64_t{UINT32_MAX}))
            ->capture_default_str();
        app.add_option("--file-reading-mode", file_reading_mode,
                       "Input reader mode: parallel or sequential")
            ->check(CLI::IsMember({"parallel", "sequential"}))
            ->capture_default_str();

        app.add_flag("--rle", cfg.rle, "Enable Run Length Encoding")
            ->capture_default_str();

        app.add_flag("--print-raw-dbg", cfg.print_precleanup_graph,
                     "Print DBG before cleanup (with .raw.gfa suffix)")
            ->capture_default_str();
        app.add_flag("--print-bi-dbg", cfg.print_prelayout_graph,
                     "Print DBG before layout (with .bi.gfa suffix)")
            ->capture_default_str();
        app.add_flag("--silent", cfg.silent, "Do not print progress to stderr")
            ->capture_default_str();
        app.add_flag("--drop-crumbs,!--keep-crumbs", cfg.drop_crumbs,
                     "Drop crumb nodes (fewer than lcmer-l+1 cores, so "
                     "invisible to overlap detection) that have a single "
                     "predecessor and continue into a longer node, instead of "
                     "emitting them as C segments. Halves the contig count at "
                     "a small genome-fraction cost.")
            ->capture_default_str();
        app.add_flag("--prime-reverse,!--no-prime-reverse", cfg.prime_reverse,
                     "Also parse each read's reverse complement and insert "
                     "its lcmer keys into the Bloom filter only (no node, "
                     "count or edge), so reads of the other strand are "
                     "admitted at their first sighting")
            ->capture_default_str();
        app.add_flag("--strict-layout-labels,!--no-strict-layout-labels",
                     cfg.strict_layout_labels,
                     "Layout: accept a core-label match as a seed only if "
                     "both cores have the same homopolymer-compressed "
                     "sequence, so cores that merely share a label are not "
                     "chained")
            ->capture_default_str();
        app.add_flag("--end-joins,!--no-end-joins", cfg.end_joins,
                     "Join segment ends that overlap as reverse complements "
                     "(the two strands of one locus stopping at different "
                     "holes) into one segment")
            ->capture_default_str();
        app.add_option("--syncmer-k", cfg.syncmer.k,
                       "Syncmer k-mer length (in RLE runs with --rle)")
            ->check(CLI::Range(2u, syncmer_params::max_k))
            ->capture_default_str();
        app.add_option(
               "--submer-s", cfg.syncmer.smer.s,
               "s-mer length s for syncmer selection. Without --submer-t "
               "this is a closed syncmer (smallest s-mer at the first or "
               "last position). Default when neither --submer-s nor "
               "--submer-t is given: open syncmer, s=4, t=5.")
            ->check(CLI::Range(1ul, 15ul));
        app.add_option(
               "--submer-t", cfg.syncmer.smer.t,
               "Offset t of the smallest s-mer: selects an open syncmer "
               "(smallest s-mer at position t).")
            ->check(CLI::Range(0ul, 15ul))
            ->needs("--submer-s");
        auto bfg = app.add_option_group("bloomfilter",
                                        "Parameters to setup bloom filter");
        bfg->add_option(
               "--projected-element-count", cfg.bloom.projected_count,
               "Distinct lcmer keys the Bloom filter is sized for at the "
               "configured false-positive probability; construction logs the "
               "fill and estimated false-positive rate and warns when the "
               "filter is undersized")
            ->check(CLI::Range((int64_t)1, (int64_t)10000000000))
            ->capture_default_str();
        bfg->add_option("--false-positive-probability",
                        cfg.bloom.false_positive_rate,
                        "False positive probability")
            ->capture_default_str();
        auto gcg = app.add_option_group("graphcleaning",
                                        "Parameters to setup graph cleaning");
        gcg->add_option("--iterations", cfg.graph_cleaning.iterations,
                        "Number of iterations of graph cleaning")
            ->check(CLI::Range(0, 100))
            ->capture_default_str();
        gcg->add_option("--max-tip-length", cfg.graph_cleaning.max_tip_length,
                        "Maximum tip length to cut in terms of cores")
            ->check(CLI::Range(0, 1000000))
            ->capture_default_str();
        gcg->add_option(
               "--min-bubble-length-difference",
               cfg.graph_cleaning.min_bubble_length_difference,
               "Maximum length ratio of a bubble's removed nodes to its kept "
               "path; more divergent bubbles are kept.")
            ->capture_default_str();
        gcg->add_option("--max-bubble-length",
                        cfg.graph_cleaning.max_bubble_length,
                        "Maximum inner length (in cores) of a bubble to pop; 0 "
                        "disables bubble popping.")
            ->check(CLI::Range(0, 1000000))
            ->capture_default_str();
        gcg->add_option("--ploidy", cfg.graph_cleaning.ploidy,
                        "Copies of a locus (1 or 2): bubbles whose nodes read "
                        "at most this many copies are popped, deeper ones are "
                        "repeats and kept; 2 also cuts allele tips at forks")
            ->check(CLI::Range(1, 2))
            ->capture_default_str();
        app.callback([this, level_opt, multipass_opt]() {
            if (olc_lcmer_l_ignored != 0) {
                std::fprintf(stderr,
                             "[config] --olc-lcmer-l is ignored; overlap "
                             "detection uses --lcmer-l (%d)\n",
                             cfg.lcmer_l);
            }
            if (multipass_opt->count() > 0) {
                auto passes = parse_multipass(multipass_text);
                if (!passes) {
                    throw CLI::ValidationError(passes.error());
                }
                cfg.multipass = std::move(*passes);
            } else if (level_opt->count() == 0) {
                throw CLI::RequiredError(level_opt->get_name());
            }
            const bool low_level = [&] {
                if (cfg.multipass.empty()) return cfg.lcp_level < 3;
                return std::any_of(
                    cfg.multipass.begin(), cfg.multipass.end(),
                    [](const pass_params &p) { return p.lcp_level < 3; });
            }();
            if (low_level) {
                std::fprintf(stderr,
                             "Warning: LCP levels lower than 3 are not "
                             "intended for actual use and may use "
                             "significantly more memory and build "
                             "suboptimal assemblies.\n");
            }
            // Validate the syncmer configuration up front; a CLI error makes
            // CLI11 print the message and exit nonzero.
            cfg.syncmer.apply_default_rule();
            if (auto v = cfg.syncmer.validate(); !v) {
                throw CLI::ValidationError(v.error());
            }
            cfg.parallel_file_reading = (file_reading_mode == "parallel");
            // CLI11 callbacks return void; surface a nonzero rc as a
            // CLI::RuntimeError so the failure reaches the process exit status.
            if (int rc = run_assembly(cfg); rc != 0) {
                throw CLI::RuntimeError("assembly failed", rc);
            }
        });
    }
    AsmProgram(const AsmProgram &) = delete;
    AsmProgram &operator=(const AsmProgram &) = delete;
};

class SamReadTesting {
    string input;
    int lcp_level = 5;
    int threads = 8;
    SamReadTestingMode mode = SamReadTestingMode::BASIC;
    string job = "read_count";
    syncmer_params syncmer;

   public:
    explicit SamReadTesting(CLI::App &parent) {
        auto &app = *parent.add_subcommand("bamdoctor", "read bam file fast");
        app.add_option("input", input, "Input bam file")->required();
        app.add_option("-l,--level", lcp_level, "LCP level to compute")
            ->required();
        app.add_option("-t,--threads", threads, "Number of threads")
            ->check(CLI::Range(1, MAXIMUM_WORK_THREADS))
            ->capture_default_str();
        app.add_option("--syncmer-k", syncmer.k, "Syncmer k-mer length")
            ->check(CLI::Range(2u, syncmer_params::max_k))
            ->capture_default_str();
        app.add_option("--submer-s", syncmer.smer.s,
                       "s-mer length (default with no --submer-t: open "
                       "syncmer, s=4, t=5)")
            ->check(CLI::Range(1ul, 15ul));
        app.add_option("--submer-t", syncmer.smer.t,
                       "Offset t of the smallest s-mer (open syncmer)")
            ->check(CLI::Range(0ul, 15ul))
            ->needs("--submer-s");

        std::map<std::string, SamReadTestingMode> SamReadTestingMode_map = {
            {"basic", SamReadTestingMode::BASIC},
            {"pool", SamReadTestingMode::POOL},
        };
        std::map<std::string, std::string> job_map = {
            {"read_count", "read_count"},
            {"core_count", "core_count"},
        };
        app.add_option("--job", job, "Job to run")
            ->required()
            ->transform(CLI::CheckedTransformer(job_map, CLI::ignore_case));
        app.add_option("--mode", mode, "Mode to run")
            ->required()
            ->transform(CLI::CheckedTransformer(SamReadTestingMode_map,
                                                CLI::ignore_case));
        app.callback([this]() {
            syncmer.apply_default_rule();
            if (auto v = syncmer.validate(); !v) {
                throw CLI::ValidationError(v.error());
            }
            // CLI11 callbacks return void; re-raise a nonzero rc as a CLI error
            // so it reaches the process exit status.
            if (int rc = run_sam_read_testing(input, threads, lcp_level, mode,
                                              job, syncmer);
                rc != 0) {
                throw CLI::RuntimeError("bamdoctor failed", rc);
            }
        });
    }
    SamReadTesting(const SamReadTesting &) = delete;
    SamReadTesting &operator=(const SamReadTesting &) = delete;
};

int main(int argc, char **argv) {
    CLI::App app("lasm");
    app.require_subcommand(1);
    // app.set_help_all_flag("--help-all", "Expand all help");

    AsmProgram hifi{app};
    SamReadTesting bamdoctor{app};
    DumpCoresProgram dumpcores{app};

    CLI11_PARSE(app, argc, argv);
    return 0;
}
