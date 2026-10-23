#pragma once

#ifndef _ASM_H_
#define _ASM_H_

#include <inttypes.h>
#include <stdint.h>

#include <algorithm>
#include <cstdio>
#include <expected>
#include <optional>
#include <string>
#include <vector>
// Syncmer selection parameters (Edgar 2021). A k-mer is selected by the
// position µ of its smallest s-mer: a closed syncmer has µ at the first or
// last s-mer, an open syncmer has µ at offset t.
struct syncmer_params {
    // Selection packs k-mers into 32-bit codes.
    static constexpr unsigned max_k = 16;
    // Rule used when neither s nor t is given: open syncmer, s = 4, t = 6.
    static constexpr unsigned long default_s = 4;
    static constexpr unsigned long default_t = 6;

    enum class typeenum { OpenSyncmer, ClosedSyncmer, Invalid };
    constexpr typeenum type() const {
        if (!smer.s.has_value()) {
            return typeenum::Invalid;
        }
        return smer.t.has_value() ? typeenum::OpenSyncmer
                                  : typeenum::ClosedSyncmer;
    }
    unsigned k = 11;
    // s-mer length s and (open syncmer) offset t.
    struct {
        std::optional<unsigned long> s;
        std::optional<unsigned long> t;
    } smer;

    // With no s and no t, select the default open syncmer. An s without a t
    // stays a closed syncmer.
    void apply_default_rule() {
        if (!smer.s.has_value() && !smer.t.has_value()) {
            smer.s = default_s;
            smer.t = default_t;
        }
    }

    // Validate the whole syncmer configuration once, before any syncmer
    // functor is constructed, so callers can surface a clean error.
    std::expected<void, std::string> validate() const {
        if (type() == typeenum::Invalid) {
            return std::unexpected(std::string("--syncmer: no s-mer length"));
        }
        if (k > max_k) {
            return std::unexpected("--syncmer: k (" + std::to_string(k) +
                                   ") must be <= " + std::to_string(max_k));
        }
        if (*smer.s == 0 || *smer.s >= k) {
            return std::unexpected("--syncmer: s (" + std::to_string(*smer.s) +
                                   ") must be in [1, k = " + std::to_string(k) +
                                   ")");
        }
        if (smer.t.has_value() && *smer.t >= k - *smer.s + 1) {
            return std::unexpected("--syncmer: t (" + std::to_string(*smer.t) +
                                   ") must be <= k - s (" +
                                   std::to_string(k - *smer.s) + ")");
        }
        return {};
    }

    // Parses "k/s[/t]" (e.g. 11/4/6; without /t a closed syncmer) and
    // validates the result.
    static std::expected<syncmer_params, std::string> parse(
        const std::string &text) {
        const std::string malformed =
            "--syncmer: expected k/s or k/s/t (e.g. 11/4/6), got \"" + text +
            "\"";
        unsigned long values[3] = {0, 0, 0};
        int count = 0;
        size_t at = 0;
        while (count < 3) {
            size_t len = 0;
            while (at + len < text.size() && text[at + len] >= '0' &&
                   text[at + len] <= '9') {
                ++len;
            }
            if (len == 0 || len > 9) return std::unexpected(malformed);
            values[count] = std::stoul(text.substr(at, len));
            ++count;
            at += len;
            if (at == text.size()) break;
            if (text[at] != '/') return std::unexpected(malformed);
            ++at;
        }
        if (count < 2 || at != text.size()) return std::unexpected(malformed);
        syncmer_params p;
        p.k = static_cast<unsigned>(values[0]);
        p.smer.s = values[1];
        if (count == 3) p.smer.t = values[2];
        if (auto v = p.validate(); !v) return std::unexpected(v.error());
        return p;
    }

    // "k/s/t", or "k/s" for a closed syncmer; the form parse() accepts.
    std::string text() const {
        std::string out = std::to_string(k) + "/" + std::to_string(*smer.s);
        if (smer.t.has_value()) out += "/" + std::to_string(*smer.t);
        return out;
    }

    void print_cfg(FILE *f) const {
        if (!smer.s.has_value()) {
            fprintf(f, "Syncmer: none\n");
        } else if (smer.t.has_value()) {
            fprintf(f, "Syncmer: %s (open)\n", text().c_str());
        } else {
            fprintf(f, "Syncmer: %s (closed)\n", text().c_str());
        }
    }
};
// Default --random-seed. The seed only drives the syncmer s-mer hash, which
// is 32-bit, so the seed is stored as uint32_t.
inline constexpr uint32_t kDefaultRandomSeed = 1453;
// One pass of --multipass: the LCP level and lcmer window length of one
// assembly pass.
struct pass_params {
    int lcp_level;
    int lcmer_l;
};
// Pass bounds: levels as -l; a window needs at least two cores to have an
// edge between consecutive windows.
inline constexpr int kMinPassLevel = 1, kMaxPassLevel = 64;
inline constexpr int kMinPassLcmer = 2, kMaxPassLcmer = 64;
// Parses "LEVEL:WINDOW[,LEVEL:WINDOW...]" (at least one pass, bounds above).
inline std::expected<std::vector<pass_params>, std::string> parse_multipass(
    const std::string &text) {
    std::vector<pass_params> out;
    size_t at = 0;
    while (true) {
        const size_t end = std::min(text.find(',', at), text.size());
        const std::string item = text.substr(at, end - at);
        const size_t colon = item.find(':');
        int level = 0, lcmer = 0;
        size_t level_len = 0, lcmer_len = 0;
        try {
            if (colon != std::string::npos) {
                level = std::stoi(item.substr(0, colon), &level_len);
                lcmer = std::stoi(item.substr(colon + 1), &lcmer_len);
            }
        } catch (const std::exception &) {
            level_len = 0;
        }
        if (colon == std::string::npos || level_len != colon ||
            lcmer_len == 0 || colon + 1 + lcmer_len != item.size()) {
            return std::unexpected(
                "--multipass: expected LEVEL:WINDOW, got \"" + item + "\"");
        }
        if (level < kMinPassLevel || level > kMaxPassLevel ||
            lcmer < kMinPassLcmer || lcmer > kMaxPassLcmer) {
            return std::unexpected("--multipass: " + item + " needs level " +
                                   std::to_string(kMinPassLevel) + "-" +
                                   std::to_string(kMaxPassLevel) +
                                   " and window " +
                                   std::to_string(kMinPassLcmer) + "-" +
                                   std::to_string(kMaxPassLcmer));
        }
        out.push_back({level, lcmer});
        if (end == text.size()) break;
        at = end + 1;
    }
    return out;
}
struct asm_config {
    std::vector<std::string> inputs;
    std::string output;
    int lcp_level = 3;
    // LCMER window length (cores) of the DBG keys; overlap detection queries
    // the same index, so it uses the same length.
    int lcmer_l = 6;
    // --multipass passes; empty = one pass at (lcp_level, lcmer_l).
    std::vector<pass_params> multipass;
    struct {
        int io = 4;
        int work = 8;
    } threads;
    uint32_t seed = kDefaultRandomSeed;
    bool parallel_file_reading = false;
    syncmer_params syncmer;
    bool print_precleanup_graph = false;
    bool print_prelayout_graph = false;
    // Drop crumb nodes (fewer than lcmer_l+1 cores) that continue into a
    // longer node instead of emitting them as C segments.
    bool drop_crumbs = true;
    // Also parse every read's reverse complement and insert its lcmer keys
    // into the Bloom filter only, so the other strand's reads are admitted at
    // their first sighting.
    bool prime_reverse = true;
    // Remember each window's first own-strand sighting until a later sighting
    // on either strand admits the window, then add that sighting's count and
    // edges, so first sightings are not lost.
    bool retro_admission = true;
    // With prime_reverse, two consecutive windows of a read's reverse
    // complement whose nodes both exist get their missing edge; no node or
    // count changes.
    bool link_reverse = false;
    // Accept a core-label match as an SCD layout seed only if both cores have
    // the same homopolymer-compressed sequence, so cores that merely share a
    // label are not chained.
    bool strict_layout_labels = true;
    // Join segment ends that overlap as reverse complements (the two strands
    // of one locus stopping at different holes) into one segment.
    bool end_joins = true;
    bool silent = false;
    struct {
        uint64_t projected_count = 100'000'000;
        double false_positive_rate = 0.001;
    } bloom;
    struct {
        int iterations = 6;
        int max_tip_length = 8 * 2;
        double min_bubble_length_difference = 2;
        // Inner cores of the largest bubble to pop; 0 disables popping.
        size_t max_bubble_length = 1024;
        // Bubbles within this many copies lie at one locus and are popped.
        int ploidy = 1;
    } graph_cleaning;

    void print_cfg(FILE *f) const {
        if (!inputs.empty()) {
            fprintf(f, "Inputs (%zu):\n", inputs.size());
            for (const auto &in : inputs) {
                fprintf(f, "\t%s\n", in.c_str());
            }
        } else {
            fprintf(f, "Inputs (0):\n");
        }
        fprintf(f, "Output: %s\n", output.c_str());
        if (multipass.empty()) {
            fprintf(f, "LCP Level: %d\n", lcp_level);
            fprintf(f, "Window: %d\n", lcmer_l);
        } else {
            fprintf(f, "Multipass:");
            const char *sep = " ";
            for (const auto &p : multipass) {
                fprintf(f, "%s%d:%d", sep, p.lcp_level, p.lcmer_l);
                sep = ",";
            }
            fprintf(f, "\n");
        }
        fprintf(f, "Threads:\n\tIO: %d\n\tWork: %d\n", threads.io, threads.work);
        fprintf(f, "Seed: %" PRIu32 "\n", seed);
        fprintf(f, "File Reading Mode: %s\n",
                parallel_file_reading ? "parallel" : "sequential");
        syncmer.print_cfg(f);
        fprintf(f,
                "Bloom Filter Parameters:\n\tProjected Count: %" PRIu64
                "\n\tFalse Positive "
                "Rate: %.3f\n",
                bloom.projected_count, bloom.false_positive_rate);
        fprintf(f, "Drop Crumbs: %s\n", drop_crumbs ? "true" : "false");
        fprintf(f, "Prime Reverse: %s\n", prime_reverse ? "true" : "false");
        fprintf(f, "Retro Admission: %s\n", retro_admission ? "true" : "false");
        fprintf(f, "Link Reverse: %s\n", link_reverse ? "true" : "false");
        fprintf(f, "Strict Layout Labels: %s\n",
                strict_layout_labels ? "true" : "false");
        fprintf(f, "End Joins: %s\n", end_joins ? "true" : "false");
        fprintf(f,
                "Graph Cleaning Parameters: \n\tIterations: %d\n\tMax Tip "
                "Length: %d\n\tMin Bubble Length Difference: %.2f\n\tMax "
                "Bubble Length: %zu\n\tPloidy: %d\n",
                graph_cleaning.iterations, graph_cleaning.max_tip_length,
                graph_cleaning.min_bubble_length_difference,
                graph_cleaning.max_bubble_length, graph_cleaning.ploidy);
    }
};

int run_assembly(const asm_config &param);
#endif
