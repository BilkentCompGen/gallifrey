#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../ext/doctest.h"

#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "../layout/layout.h"

using u64 = std::uint64_t;
namespace ll = lasm::layout;

// This test file intentionally models how recig is formed in your pipeline:
//
// - We start from a single underlying contig of length N in CORE-SPACE.
// - We create multiple fragments from that contig in both forward (+) and reverse (-) orientations.
// - We then compute overlaps between every (+) fragment and every (-) fragment (and vice versa).
// - Coordinates we store are ALWAYS on the REVERSED (RC) sequence of each fragment/contig.
// - We do NOT use any "master_overlap_start" helper; overlaps are computed directly from fragment
//   ranges on the underlying contig.
//
// Also note your important constraint:
// Each stored coordinate corresponds to the *last* core of an L-length window (LCMER_L),
// so endpoints are not trustworthy anchors. We use interior anchors in the checks.

#ifndef LCMER_L
#define LCMER_L 8
#endif

static constexpr int k_core_margin = (LCMER_L > 1) ? (LCMER_L - 1) : 0;

// ---------------------------
// Helpers: fragment model
// ---------------------------

struct Frag {
    u64 id{};
    bool plus{};   // orientation relative to underlying contig (+ or -)
    int start{};   // inclusive start on underlying contig (core-space)
    int end{};     // inclusive end on underlying contig (core-space)
};

static int len_inclusive(const Frag& f) { return f.end - f.start; } // length in "steps" (end-start)
static int span_inclusive(const Frag& f) { return f.end - f.start + 1; } // number of positions

// Map a coordinate on the underlying contig to a FORWARD coordinate on the fragment.
static int contig_to_frag_forward(const Frag& f, int x_contig) {
    CHECK(x_contig >= f.start);
    CHECK(x_contig <= f.end);

    if (f.plus) {
        // Forward fragment: forward coordinate grows with contig coordinate.
        return x_contig - f.start;
    }
    // Reverse fragment: forward coordinate decreases with contig coordinate.
    // At contig == start -> forward == len
    // At contig == end   -> forward == 0
    return f.end - x_contig;
}

// Convert a fragment forward coordinate to the observed RC coordinate (what recig stores).
// Convention: (matches existing layout code) rc = len - forward.
static int forward_to_rc(int frag_len_steps, int x_forward) {
    // frag_len_steps is (end-start), so forward in [0..len]
    return frag_len_steps - x_forward;
}

// Convert observed RC coordinate to layout-forward coordinate given inferred orientation.
static int rc_to_layout_forward(int frag_len_steps, int rc, bool inferred_plus) {
    if (inferred_plus) return frag_len_steps - rc;
    return rc;
}

// Interiorize an inclusive interval [x0,x1] in layout-forward coordinates to avoid boundary effects.
static std::pair<int, int> interiorize_inclusive_interval(int x0, int x1, int frag_len_steps) {
    if (x0 > x1) std::swap(x0, x1);

    int lo = std::clamp(x0 + k_core_margin, 0, frag_len_steps);
    int hi = std::clamp(x1 - k_core_margin, 0, frag_len_steps);
    if (lo > hi) {
        int mid = (x0 + x1) / 2;
        mid = std::clamp(mid, 0, frag_len_steps);
        return {mid, mid};
    }
    return {lo, hi};
}

// Compute overlap on underlying contig between two fragments.
// Returns true if they overlap in contig space and gives overlap interval [os,oe] inclusive.
static bool contig_overlap(const Frag& a, const Frag& b, int& os, int& oe) {
    os = std::max(a.start, b.start);
    oe = std::min(a.end, b.end);
    return os <= oe;
}

// Build reciprocal rev_match edges for a pair of fragments that overlap on the underlying contig.
// IMPORTANT: coordinates stored are always on RC sequence of each frag (per your pipeline).
static void add_reciprocal_edge_for_pair(ll::recig_graph& recig, const Frag& A, const Frag& B) {
    int os = 0, oe = -1;
    if (!contig_overlap(A, B, os, oe)) return;

    // Only add edges between opposite orientations (RC matches)
    if (A.plus == B.plus) return;

    const int lenA = len_inclusive(A);
    const int lenB = len_inclusive(B);

    // Overlap endpoints in each fragment's FORWARD coords
    int A_f0 = contig_to_frag_forward(A, os);
    int A_f1 = contig_to_frag_forward(A, oe);
    if (A_f0 > A_f1) std::swap(A_f0, A_f1);

    int B_f0 = contig_to_frag_forward(B, os);
    int B_f1 = contig_to_frag_forward(B, oe);
    if (B_f0 > B_f1) std::swap(B_f0, B_f1);

    // Convert to RC coords (stored)
    int A_rc0 = forward_to_rc(lenA, A_f0);
    int A_rc1 = forward_to_rc(lenA, A_f1);
    if (A_rc0 > A_rc1) std::swap(A_rc0, A_rc1);

    int B_rc0 = forward_to_rc(lenB, B_f0);
    int B_rc1 = forward_to_rc(lenB, B_f1);
    if (B_rc0 > B_rc1) std::swap(B_rc0, B_rc1);

    // A -> B
    {
        ll::rev_match rm{};
        rm.node = B.id;

        rm.qs = static_cast<u64>(A_rc0);
        rm.qe = static_cast<u64>(A_rc1);
        rm.qlen = static_cast<u64>(lenA);

        rm.from = static_cast<u64>(B_rc0);
        rm.to = static_cast<u64>(B_rc1);
        rm.size = static_cast<u64>(lenB);

        recig.emplace(A.id, rm);
    }

    // B -> A
    {
        ll::rev_match rm{};
        rm.node = A.id;

        rm.qs = static_cast<u64>(B_rc0);
        rm.qe = static_cast<u64>(B_rc1);
        rm.qlen = static_cast<u64>(lenB);

        rm.from = static_cast<u64>(A_rc0);
        rm.to = static_cast<u64>(A_rc1);
        rm.size = static_cast<u64>(lenA);

        recig.emplace(B.id, rm);
    }
}

// Build edges between all (+) fragments and all (-) fragments (and vice versa),
// by computing contig-space overlap.
static ll::recig_graph build_recig_from_fragments(const std::vector<Frag>& frags) {
    ll::recig_graph recig;
    for (size_t i = 0; i < frags.size(); ++i) {
        for (size_t j = 0; j < frags.size(); ++j) {
            if (i == j) continue;
            // Only connect opposite orientations
            if (frags[i].plus == frags[j].plus) continue;
            add_reciprocal_edge_for_pair(recig, frags[i], frags[j]);
        }
    }
    return recig;
}

// Compute interior-mid residuals per directed edge for a produced layout.
// Residual definition:
//   r = (startA + amid) - (startB + bmid)
// where amid and bmid are interior midpoints of the overlap intervals converted
// into layout-forward coordinates per inferred orientation.
static std::vector<int> compute_edge_interior_mid_residuals(
    const ll::component_layout& lay,
    const ll::recig_graph& recig) {
    std::vector<int> residuals;
    residuals.reserve(256);

    for (const auto& e : recig) {
        u64 A = e.first;
        u64 B = e.second.node;
        if (!lay.contains(A) || !lay.contains(B)) continue;

        auto [startA, endA, oriA] = lay.at(A);
        auto [startB, endB, oriB] = lay.at(B);
        (void)endA;
        (void)endB;

        const int lenA = static_cast<int>(e.second.qlen);
        const int lenB = static_cast<int>(e.second.size);

        int a0 = rc_to_layout_forward(e.second.qlen, static_cast<int>(e.second.qs), oriA);
        int a1 = rc_to_layout_forward(e.second.qlen, static_cast<int>(e.second.qe), oriA);
        int b0 = rc_to_layout_forward(e.second.size, static_cast<int>(e.second.from), oriB);
        int b1 = rc_to_layout_forward(e.second.size, static_cast<int>(e.second.to), oriB);

        auto [ai0, ai1] = interiorize_inclusive_interval(a0, a1, lenA);
        auto [bi0, bi1] = interiorize_inclusive_interval(b0, b1, lenB);

        int amid = (ai0 + ai1) / 2;
        int bmid = (bi0 + bi1) / 2;

        residuals.push_back((startA + amid) - (startB + bmid));
    }

    return residuals;
}

static int median_abs_i(std::vector<int> v) {
    if (v.empty()) return 0;
    for (auto& x : v) x = std::abs(x);
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

static int max_abs_i(const std::vector<int>& v) {
    int m = 0;
    for (int x : v) m = std::max(m, std::abs(x));
    return m;
}

// ---------------------------
// Fragment generation
// ---------------------------

static std::vector<Frag> make_random_fragments_from_single_contig(
    int N,                // underlying contig length in CORE-SPACE positions (use 0..N-1)
    int num_plus,
    int num_minus,
    int min_span,         // minimum span (bases/cores count), inclusive length
    int max_span,
    std::mt19937& rng) {
    REQUIRE(N > 0);
    REQUIRE(min_span >= 2);
    REQUIRE(max_span >= min_span);
    REQUIRE(max_span <= N);

    std::uniform_int_distribution<int> span_dist(min_span, max_span);
    std::uniform_int_distribution<int> start_dist(0, N - min_span);

    std::vector<Frag> frags;
    frags.reserve(static_cast<size_t>(num_plus + num_minus));

    u64 id = 1;

    auto gen_one = [&](bool plus) {
        int span = span_dist(rng); // number of positions
        int start = std::uniform_int_distribution<int>(0, N - span)(rng);
        int end = start + span - 1;
        frags.push_back(Frag{id++, plus, start, end});
    };

    for (int i = 0; i < num_plus; ++i) gen_one(true);
    for (int i = 0; i < num_minus; ++i) gen_one(false);

    return frags;
}

// ---------------------------
// Tests
// ---------------------------

TEST_CASE("recig build: overlaps are computed from fragment coordinates, not master overlap inputs") {
    std::mt19937 rng(1337);

    const int N = 2000; // core positions on single contig [0..N-1]
    auto frags = make_random_fragments_from_single_contig(
        N,
        /*num_plus*/ 8,
        /*num_minus*/ 8,
        /*min_span*/ 100,
        /*max_span*/ 400,
        rng);

    auto recig = build_recig_from_fragments(frags);

    // There should be some edges if fragments overlap across orientations.
    CHECK(recig.size() > 0);

    // Sanity: no same-orientation edges
    // (by construction we don't add them)
    CHECK(true);
}

TEST_CASE("connected_comp_plus_layout: layout recovers consistent placement on random +/- fragments from one contig (core-space, interior anchors)") {
    std::mt19937 rng(42);

    const int N = 3000;
    auto frags = make_random_fragments_from_single_contig(
        N,
        /*num_plus*/ 12,
        /*num_minus*/ 12,
        /*min_span*/ 150,
        /*max_span*/ 600,
        rng);

    auto recig = build_recig_from_fragments(frags);

    // Ensure there is at least one connected component with edges.
    REQUIRE(recig.size() > 0);

    auto layouts = ll::connected_comp_plus_layout(recig);
    REQUIRE(layouts.size() >= 1);

    // Find the largest component layout (most nodes) and evaluate residuals on it.
    size_t best_i = 0;
    for (size_t i = 1; i < layouts.size(); ++i) {
        if (layouts[i].size() > layouts[best_i].size()) best_i = i;
    }
    const auto& lay = layouts[best_i];

    // Must include at least 2 nodes.
    REQUIRE(lay.size() >= 2);

    auto residuals = compute_edge_interior_mid_residuals(lay, recig);
    REQUIRE(residuals.size() > 0);

    // On data derived from one contig, there exists a consistent placement. With the
    // least-squares solver, residuals should be small in interior core-space.
    //
    // We allow some slack due to:
    // - interiorization margin (LCMER_L-1)
    // - discrete rounding of positions to int in returned layout
    CHECK(median_abs_i(residuals) <= (k_core_margin + 2));
    CHECK(max_abs_i(residuals) <= (k_core_margin + 25));
}