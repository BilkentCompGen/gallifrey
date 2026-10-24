#include <CLI/CLI.hpp>
#include <string>
#include <string_view>
#include <vector>

#define LCPTOOLS_IMPL
#include "asm.h"
#include "config.h"
#include "dump_cores.h"
#include "ext/lcptools/lcptools_ho.h"
#include "syncmer_option.h"
using std::string;

// Lists an option's names without the "{false}" CLI11 appends to negated flag
// names (--no-end-joins{false}), which reads as the flag's default.
class HelpFormatter : public CLI::Formatter {
   public:
    std::string make_option_name(const CLI::Option *opt,
                                 bool is_positional) const override {
        if (is_positional) {
            return CLI::Formatter::make_option_name(opt, true);
        }
        std::vector<std::string> names;
        for (const auto &s : opt->get_snames()) names.push_back("-" + s);
        for (const auto &l : opt->get_lnames()) names.push_back("--" + l);
        return CLI::detail::join(names);
    }
};

// CLI11 shows no default for flags, so negatable flags state it in their help.
static std::string with_default(const std::string &desc, bool on) {
    if (on) return desc + " (default: on)";
    return desc + " (default: off)";
}

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
    // show_advanced: --help-all is on the command line. CLI11 hides options
    // whose group is "", so the advanced groups get headings only then; they
    // parse either way.
    AsmProgram(CLI::App &parent, bool show_advanced) {
        auto &app = *parent.add_subcommand("hifi", "run hifi assembly");
        app.set_help_all_flag("--help-all",
                              "Print help including advanced options and exit");
        auto advanced = [show_advanced](const std::string &heading) {
            if (show_advanced) return heading;
            return std::string{};
        };
        app.add_option("input", cfg.inputs,
                       "Input read files: BAM or FASTA/FASTQ "
                       "(plain or gzipped). Format is detected from content; "
                       "formats may be mixed. - reads standard input.")
            ->required();
        app.add_option("-o,--output", cfg.output, "Output gfa file")
            ->required();
        app.add_option("-t,--threads", cfg.threads.work,
                       "Number of worker threads")
            ->check(CLI::Range(1, MAXIMUM_WORK_THREADS))
            ->capture_default_str();
        app.add_option(
               "-@,--io-threads", cfg.threads.io,
               "Number of decompression threads per BGZF input (BAM, bgzipped "
               "FASTA/FASTQ)")
            ->check(CLI::Range(0, 64))
            ->capture_default_str();
        app.add_option("--file-reading-mode", file_reading_mode,
                       "Input reader mode: parallel or sequential")
            ->check(CLI::IsMember({"parallel", "sequential"}))
            ->capture_default_str();
        auto *level_opt =
            app.add_option("-l,--level", cfg.lcp_level, "LCP level to compute")
                ->check(CLI::Range(1, 64));
        auto *window_opt =
            app.add_option("-w,--window", cfg.lcmer_l,
                           "Cores per lcmer window, used for DBG construction "
                           "and overlap detection")
                ->check(CLI::Range(1, 64))
                ->capture_default_str();
        auto *multipass_opt =
            app.add_option(
                   "--multipass", multipass_text,
                   "Multi-pass assembly: LEVEL:WINDOW[,LEVEL:WINDOW...] gives "
                   "the LCP level and cores per lcmer window of each pass, "
                   "e.g. 4:8,5:6,5:10 (replaces -l and -w). Each pass "
                   "assembles the reads together with the previous pass's "
                   "contigs; only the last pass writes the output. Order "
                   "the passes from small to large windows.")
                ->excludes(level_opt)
                ->excludes(window_opt);
        app.add_option("--ploidy", cfg.graph_cleaning.ploidy,
                       "Ploidy of the sample: 1 for a haploid (or fully "
                       "homozygous) genome, 2 for a diploid one")
            ->check(CLI::Range(1, 2))
            ->capture_default_str();
        app.add_flag("--silent", cfg.silent, "Do not print progress to stderr");
        app.add_option("--olc-lcmer-l", olc_lcmer_l_ignored,
                       "Ignored; overlap detection uses --window")
            ->group("");

        add_syncmer_option(app, cfg.syncmer);
        app.add_option("--random-seed", cfg.seed,
                       "Seed of the syncmer s-mer hash (32-bit)")
            ->check(CLI::Range(uint64_t{0}, uint64_t{UINT32_MAX}))
            ->capture_default_str();

        const std::string bloom_group = advanced("Bloom filter");
        app.add_option(
               "--projected-element-count", cfg.bloom.projected_count,
               "Distinct lcmer keys the Bloom filter is sized for at the "
               "configured false-positive probability; construction logs the "
               "fill and modelled false-positive rate and warns when the "
               "filter is undersized")
            ->check(CLI::Range((int64_t)1, (int64_t)10000000000))
            ->capture_default_str()
            ->group(bloom_group);
        app.add_option("--false-positive-probability",
                       cfg.bloom.false_positive_rate,
                       "False positive probability")
            ->capture_default_str()
            ->group(bloom_group);
        app.add_flag("--prime-reverse,!--no-prime-reverse", cfg.prime_reverse,
                     with_default(
                         "Also parse each read's reverse complement and insert "
                         "its lcmer keys into the Bloom filter only (no node, "
                         "count or edge), so reads of the other strand are "
                         "admitted at their first sighting",
                         cfg.prime_reverse))
            ->group(bloom_group);

        const std::string cleaning_group = advanced("Graph cleaning");
        app.add_option("--iterations", cfg.graph_cleaning.iterations,
                       "Number of iterations of graph cleaning")
            ->check(CLI::Range(0, 100))
            ->capture_default_str()
            ->group(cleaning_group);
        app.add_option("--max-tip-length", cfg.graph_cleaning.max_tip_length,
                       "Maximum tip length to cut in terms of cores")
            ->check(CLI::Range(0, 1000000))
            ->capture_default_str()
            ->group(cleaning_group);
        app.add_option("--min-bubble-length-difference",
                       cfg.graph_cleaning.min_bubble_length_difference,
                       "Maximum length ratio of a bubble's removed nodes to "
                       "its kept path; more divergent bubbles are kept.")
            ->capture_default_str()
            ->group(cleaning_group);
        app.add_option("--max-bubble-length",
                       cfg.graph_cleaning.max_bubble_length,
                       "Maximum inner length (in cores) of a bubble to pop; 0 "
                       "disables bubble popping.")
            ->check(CLI::Range(0, 1000000))
            ->capture_default_str()
            ->group(cleaning_group);

        const std::string layout_group = advanced("Layout");
        app.add_flag("--drop-crumbs,!--keep-crumbs", cfg.drop_crumbs,
                     with_default(
                         "Drop crumb nodes (fewer than window+1 cores, so "
                         "invisible to overlap detection) that have a single "
                         "predecessor and continue into a longer node, instead "
                         "of emitting them as C segments. Halves the contig "
                         "count at a small genome-fraction cost.",
                         cfg.drop_crumbs))
            ->group(layout_group);
        app.add_flag(
               "--strict-layout-labels,!--no-strict-layout-labels",
               cfg.strict_layout_labels,
               with_default("Accept a core-label match as a seed only if both "
                            "cores have the same homopolymer-compressed "
                            "sequence, so cores that merely share a label are "
                            "not chained",
                            cfg.strict_layout_labels))
            ->group(layout_group);
        app.add_flag("--end-joins,!--no-end-joins", cfg.end_joins,
                     with_default(
                         "Join segment ends that overlap as reverse "
                         "complements (the two strands of one locus stopping "
                         "at different holes) into one segment",
                         cfg.end_joins))
            ->group(layout_group);

        const std::string other_group = advanced("Other");
        app.add_flag("--print-raw-dbg", cfg.print_precleanup_graph,
                     "Print DBG before cleanup (with .raw.gfa suffix)")
            ->group(other_group);
        app.add_flag("--print-bi-dbg", cfg.print_prelayout_graph,
                     "Print DBG before layout (with .bi.gfa suffix)")
            ->group(other_group);
        app.add_flag("--retro-admission,!--no-retro-admission",
                     cfg.retro_admission,
                     with_default(
                         "Remember each lcmer window's first sighting on a "
                         "read's own strand until a later sighting on either "
                         "strand admits the window, then add that sighting's "
                         "count and edges as if it had been admitted at once",
                         cfg.retro_admission))
            ->group(other_group);
        app.add_flag("--link-reverse", cfg.link_reverse,
                     with_default(
                         "Let each read's reverse complement add the missing "
                         "links between lcmer windows that are already in the "
                         "graph, so the two strands' graphs fill each other's "
                         "gaps",
                         cfg.link_reverse))
            ->group(other_group);
        app.callback([this, level_opt, multipass_opt]() {
            if (olc_lcmer_l_ignored != 0) {
                std::fprintf(stderr,
                             "[config] --olc-lcmer-l is ignored; overlap "
                             "detection uses --window (%d)\n",
                             cfg.lcmer_l);
            }
            if (cfg.link_reverse && !cfg.prime_reverse) {
                std::fprintf(stderr,
                             "[config] --link-reverse needs --prime-reverse; "
                             "no links are added\n");
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

int main(int argc, char **argv) {
    CLI::App app("lasm");
    app.require_subcommand(1);
    app.formatter(std::make_shared<HelpFormatter>());

    bool show_advanced = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string_view(argv[i]) == "--help-all") show_advanced = true;
    }
    AsmProgram hifi{app, show_advanced};
    DumpCoresProgram dumpcores{app};

    CLI11_PARSE(app, argc, argv);
    return 0;
}
