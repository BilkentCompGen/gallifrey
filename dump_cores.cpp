#include "dump_cores.h"

#include <inttypes.h>
#include <stdint.h>

#include <CLI/CLI.hpp>
#include <cstdio>
#include <string>

#include "asm.h"
#include "io/fastx.h"
#include "lcmer.h"
#include "syncmer.h"
#include "syncmer_option.h"

static int run_dump_cores(const std::string &input_path,
                          const std::string &output_path,
                          const syncmer_params &spara, int lcp_level,
                          int lcmer_l, uint32_t seed) {
    using SYNCMER_HASHER = syncmer::MurMur;
    using SyncmerType = syncmer::GenericSubmerRule<uint32_t, SYNCMER_HASHER>;
    SyncmerType syncmer_functor(spara, seed);

    FILE *out = fopen(output_path.c_str(), "w");
    if (!out) {
        fprintf(stderr, "Cannot open output: %s\n", output_path.c_str());
        return 1;
    }
    if (lcmer_l > 0)
        fprintf(out,
                "read_name\tread_len\twin_idx\tlcmer\tstart\tend\tfrag_start\t"
                "frag_end\n");
    else
        fprintf(
            out,
            "read_name\tread_len\tcore_idx\tcore_label\tstart\tend\tspan\n");

    input_stream in;
    if (!in.open(input_path, 0)) {
        fprintf(stderr, "Cannot open input: %s\n", input_path.c_str());
        fclose(out);
        return 1;
    }
    const auto fmt = detect_input_format(in);
    if (!fmt || *fmt != input_format::fastx) {
        fprintf(stderr, "Not a FASTA/FASTQ input: %s\n", input_path.c_str());
        fclose(out);
        return 1;
    }
    fastx_reader reader(in);

    seq_record rec;
    core_buf str;
    charsi chs;
    pos_scratch pos_buf;
    int read_count = 0;

    while (reader.next(rec)) {
        lcp_read lr;
        parse_cores_from_nibbles(rec.nibbles, rec.len, lcp_level, &str,
                                 syncmer_functor, chs, lr, pos_buf);

        if (lcmer_l > 0) {
            // Window w holds cores w..w+L-1 and is keyed like the DBG node of
            // that window; its node fragment is the bases its last core adds.
            if (lr.cores.size() < static_cast<size_t>(lcmer_l)) {
                ++read_count;
                continue;
            }
            LCmer l{lr.cores.begin(), lcmer_l};
            for (size_t w = 0; w + lcmer_l <= lr.cores.size(); ++w) {
                if (w > 0)
                    l.update_lcmer(lr.cores[w + lcmer_l - 1], lr.cores[w - 1],
                                   lcmer_l);
                const size_t last = w + lcmer_l - 1;
                fprintf(out, "%s\t%u\t%zu\t%" PRIu64 "\t%u\t%u\t%u\t%u\n",
                        reader.name().c_str(), rec.len, w, l.value(),
                        lr.pos[w].first, lr.pos[last].second,
                        last > 0 ? lr.pos[last - 1].second : lr.pos[w].first,
                        lr.pos[last].second);
            }
            ++read_count;
            continue;
        }
        for (size_t i = 0; i < lr.cores.size(); ++i) {
            fprintf(out, "%s\t%u\t%zu\t%" PRIu32 "\t%u\t%u\t%u\n",
                    reader.name().c_str(), rec.len, i, lr.cores[i],
                    lr.pos[i].first, lr.pos[i].second,
                    lr.pos[i].second - lr.pos[i].first);
        }
        ++read_count;
    }
    free(rec.nibbles);
    if (reader.failed()) {
        fprintf(stderr, "Malformed FASTA/FASTQ input: %s\n",
                input_path.c_str());
        fclose(out);
        return 1;
    }

    fclose(out);
    fprintf(stderr, "Processed %d reads, output written to %s\n", read_count,
            output_path.c_str());
    return 0;
}

DumpCoresProgram::DumpCoresProgram(CLI::App &parent) {
    auto &app =
        *parent.add_subcommand("dump-cores", "Read FASTA and dump LCP cores");
    app.add_option("input", input, "Input FASTA file")->required();
    app.add_option("-o,--output", output, "Output TSV file")->required();
    add_syncmer_option(app, syncmer);
    app.add_option("-l,--level", lcp_level, "LCP level to compute")
        ->capture_default_str();
    app.add_option("-w,--window", lcmer_l,
                   "Dump lcmer keys of windows of this many cores (DBG node "
                   "keys) instead of cores; 0 dumps cores")
        ->capture_default_str();
    app.add_option("--seed", seed, "Seed of the syncmer s-mer hash (32-bit)")
        ->check(CLI::Range(uint64_t{0}, uint64_t{UINT32_MAX}))
        ->capture_default_str();
    app.callback([this]() {
        // CLI11 callbacks return void; surface a nonzero rc as a
        // CLI::RuntimeError so CLI11_PARSE prints it and exits with that code.
        if (int rc = run_dump_cores(input, output, syncmer, lcp_level, lcmer_l,
                                    seed);
            rc != 0) {
            throw CLI::RuntimeError("dump-cores failed", rc);
        }
    });
}
