#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../ext/doctest.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include <gtl/btree.hpp>
#include <gtl/phmap.hpp>
#include <gtl/vector.hpp>

#include "../layout/layout.h"

using u64 = std::uint64_t;

// Layout types now come from the extracted implementation.
namespace ll = lasm::layout;
using rev_match = ll::rev_match;

// ---------------------------
// Spec / synthetic generator
// ---------------------------

// A synthetic fragment/copy derived from an underlying (conceptual) long sequence.
struct Frag {
    u64 id;
    int contig_len;      // length of this contig (bases)
    int contig_origin;   // "true" start on the conceptual master coordinate line
    bool orient_plus;    // "true" orientation relative to the master (+ or -)
};

// Map a master coordinate x (on the conceptual line) to the *forward* coordinate on the fragment.
static int master_to_forward(const Frag& f, int master_x) {
    // fragment covers [origin, origin+len) on master
    int rel = master_x - f.contig_origin;
    CHECK(rel >= 0);
    CHECK(rel <= f.contig_len);
    if (f.orient_plus) {
        return rel;
    }
    // If fragment is reverse relative to master, forward coordinate runs opposite.
    // At master_x == origin, forward is len.
    // At master_x == origin+len, forward is 0.
    return f.contig_len - rel;
}

// Convert a forward coordinate to a reverse-complement coordinate (same length).
static int forward_to_revcoord(int len, int forward_x) {
    return len - forward_x;
}

// Helper encoding: for an interval [f0, f1) in forward space, the corresponding
// reverse-complement interval is [len - f1, len - f0) (order preserved).
static std::pair<int, int> forward_interval_to_rev_interval(int len, int f0, int f1) {
    CHECK(f0 <= f1);
    return {len - f1, len - f0};
}

// Build reciprocal rev_match edges from a master-overlap interval.
// We deliberately generate *intervals*, not a single point.
// We also support adding small endpoint jitter to emulate approximate mappings.
static void add_reciprocal_overlap_interval(
    ll::recig_graph& recig,
    const Frag& A,
    const Frag& B,
    int master_overlap_start,
    int master_overlap_end,
    int jitter = 0,
    std::mt19937* rng = nullptr
) {
    REQUIRE(A.id != B.id);
    // Edges represent reverse-complement matches => orientations must differ in truth.
    REQUIRE(A.orient_plus != B.orient_plus);

    REQUIRE(master_overlap_start < master_overlap_end);

    // Convert overlap endpoints to forward coords on each fragment
    int a_f0 = master_to_forward(A, master_overlap_start);
    int a_f1 = master_to_forward(A, master_overlap_end);
    if (a_f0 > a_f1) std::swap(a_f0, a_f1);

    int b_f0 = master_to_forward(B, master_overlap_start);
    int b_f1 = master_to_forward(B, master_overlap_end);
    if (b_f0 > b_f1) std::swap(b_f0, b_f1);

    // Optional jitter (simulate approximate interval boundaries)
    auto apply_jitter = [&](int& x0, int& x1, int len) {
        if (jitter <= 0) return;
        REQUIRE(rng != nullptr);
        std::uniform_int_distribution<int> dist(-jitter, jitter);
        int j0 = dist(*rng);
        int j1 = dist(*rng);
        x0 = std::clamp(x0 + j0, 0, len);
        x1 = std::clamp(x1 + j1, 0, len);
        if (x0 > x1) std::swap(x0, x1);
        if (x0 == x1) {
            // Ensure non-empty when possible
            if (x1 < len) x1++;
            else if (x0 > 0) x0--;
        }
    };

    apply_jitter(a_f0, a_f1, A.contig_len);
    apply_jitter(b_f0, b_f1, B.contig_len);

    // Convert forward intervals to rev-comp intervals (what your pipeline stores)
    auto [a_rc0, a_rc1] = forward_interval_to_rev_interval(A.contig_len, a_f0, a_f1);
    auto [b_rc0, b_rc1] = forward_interval_to_rev_interval(B.contig_len, b_f0, b_f1);

    // Build A->B edge
    {
        rev_match rm{};
        rm.node = B.id;

        // query (A) interval in rev-comp
        rm.qs = static_cast<u64>(a_rc0);
        rm.qe = static_cast<u64>(a_rc1);
        rm.qlen = static_cast<u64>(A.contig_len);

        // target (B) interval in rev-comp
        rm.from = static_cast<u64>(b_rc0);
        rm.to = static_cast<u64>(b_rc1);
        rm.size = static_cast<u64>(B.contig_len);

        recig.emplace(A.id, rm);
    }

    // Build B->A reciprocal
    {
        rev_match rm{};
        rm.node = A.id;

        rm.qs = static_cast<u64>(b_rc0);
        rm.qe = static_cast<u64>(b_rc1);
        rm.qlen = static_cast<u64>(B.contig_len);

        rm.from = static_cast<u64>(a_rc0);
        rm.to = static_cast<u64>(a_rc1);
        rm.size = static_cast<u64>(A.contig_len);

        recig.emplace(B.id, rm);
    }
}

// ---------------------------
// Spec checks (do NOT assume current algorithm is correct)
// ---------------------------

// Convert a rev-comp interval to a representative "layout-forward coordinate" given an inferred orientation.
// We use midpoints to avoid sensitivity to endpoint conventions.
//
// If inferred orient is '+':
//   forward = len - rc_mid
// If inferred orient is '-':
//   forward = rc_mid
static double rc_mid_to_layout_forward(u64 len, u64 rc0, u64 rc1, bool inferred_plus) {
    double mid_rc = (static_cast<double>(rc0) + static_cast<double>(rc1)) / 2.0;
    if (inferred_plus) return static_cast<double>(len) - mid_rc;
    return mid_rc;
}

// Given a recovered layout for a component, compute per-edge residuals:
//   (posA + fA(mid)) - (posB + fB(mid))  should be ~= 0
// using the reciprocal edge's midpoints.
static std::vector<double> compute_edge_residuals(
    const gtl::flat_hash_map<u64, std::pair<int, bool>>& lay,
    const gtl::btree_multimap<u64, rev_match>& recig
) {
    std::vector<double> residuals;
    residuals.reserve(256);

    for (const auto& e : recig) {
        u64 A = e.first;
        u64 B = e.second.node;
        if (!lay.contains(A) || !lay.contains(B)) continue;

        auto [posA, oriA] = lay.at(A);
        auto [posB, oriB] = lay.at(B);

        // Do NOT hard-assert flip here. In inconsistent graphs (e.g., odd cycles or noisy/spurious overlaps),
        // the layout may violate flip constraints on some edges. We still want to compute residuals so we
        // can quantify inconsistency rather than failing early.
        //
        // From e (A->B), we have A.rc interval = [qs,qe], B.rc interval = [from,to]
        double a_f = rc_mid_to_layout_forward(e.second.qlen, e.second.qs, e.second.qe, oriA);
        double b_f = rc_mid_to_layout_forward(e.second.size, e.second.from, e.second.to, oriB);

        double lhs = static_cast<double>(posA) + a_f;
        double rhs = static_cast<double>(posB) + b_f;
        residuals.push_back(lhs - rhs);
    }
    return residuals;
}

static double median_abs(std::vector<double> v) {
    if (v.empty()) return 0.0;
    for (auto& x : v) x = std::abs(x);
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

static double max_abs(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) m = std::max(m, std::abs(x));
    return m;
}

// ---------------------------
// Tests
// ---------------------------

TEST_CASE("spec: empty graph yields no layouts") {
    ll::recig_graph recig;
    auto layouts = ll::connected_comp_plus_layout(recig);
    CHECK(layouts.size() == 0);
}

TEST_CASE("spec: simple chain of +/- fragments produces one component and low interval residuals") {
    // Master axis is conceptual. We create 4 fragments overlapping in a chain with alternating orientation.
    //
    // A(+): [0,  1000)
    // B(-): [400,1400)
    // C(+): [900,1900)
    // D(-): [1300,2300)
    Frag A{1, 1000, 0, true};
    Frag B{2, 1000, 400, false};
    Frag C{3, 1000, 900, true};
    Frag D{4, 1000, 1300, false};

    ll::recig_graph recig;

    add_reciprocal_overlap_interval(recig, A, B, 450, 800);   // overlap within [400,1000)
    add_reciprocal_overlap_interval(recig, B, C, 950, 1200);  // overlap within [900,1400)
    add_reciprocal_overlap_interval(recig, C, D, 1350, 1700); // overlap within [1300,1900)

    auto layouts = ll::connected_comp_plus_layout(recig);
    REQUIRE(layouts.size() == 1);

    const auto& lay = layouts[0];
    REQUIRE(lay.size() == 4);

    // Behavioral check: orientations alternate and satisfy flip edges.
    CHECK(lay.at(A.id).second != lay.at(B.id).second);
    CHECK(lay.at(B.id).second != lay.at(C.id).second);
    CHECK(lay.at(C.id).second != lay.at(D.id).second);

    auto residuals = compute_edge_residuals(lay, recig);
    REQUIRE(!residuals.empty());

    // With perfect synthetic intervals (no jitter), residuals should be very small (0 with an ideal solver).
    // We allow small tolerance because algorithm may pick different endpoints effectively.
    CHECK(median_abs(residuals) <= 1.0);
    CHECK(max_abs(residuals) <= 5.0);
}

TEST_CASE("spec: branched graph remains consistent (uses multiple reciprocal overlaps)") {
    // Build a component with branching constraints.
    //
    // N1(+): [0,  800)
    // N2(-): [300,1100)
    // N3(+): [600,1400)
    // N4(-): [900,1700)
    // N5(+): [400,1200)
    Frag N1{10, 800, 0, true};
    Frag N2{20, 800, 300, false};
    Frag N3{30, 800, 600, true};
    Frag N4{40, 800, 900, false};
    Frag N5{50, 800, 400, true};

    ll::recig_graph recig;

    add_reciprocal_overlap_interval(recig, N1, N2, 350, 650);
    add_reciprocal_overlap_interval(recig, N2, N3, 650, 1000);
    add_reciprocal_overlap_interval(recig, N3, N4, 950, 1300);
    add_reciprocal_overlap_interval(recig, N5, N2, 450, 750);
    add_reciprocal_overlap_interval(recig, N5, N4, 950, 1150);

    auto layouts = ll::connected_comp_plus_layout(recig);
    REQUIRE(layouts.size() == 1);
    const auto& lay = layouts[0];
    REQUIRE(lay.size() == 5);

    auto residuals = compute_edge_residuals(lay, recig);
    REQUIRE(residuals.size() >= 5); // at least one per inserted pair direction

    CHECK(median_abs(residuals) <= 1.0);
    CHECK(max_abs(residuals) <= 10.0);
}

TEST_CASE("spec: multiple components yield multiple layouts; each is internally consistent") {
    Frag A1{1, 500, 0, true};
    Frag A2{2, 500, 200, false};

    Frag B1{3, 500, 1000, true};
    Frag B2{4, 500, 1200, false};

    ll::recig_graph recig;
    add_reciprocal_overlap_interval(recig, A1, A2, 250, 450);
    add_reciprocal_overlap_interval(recig, B1, B2, 1250, 1450);

    auto layouts = ll::connected_comp_plus_layout(recig);
    REQUIRE(layouts.size() == 2);

    for (const auto& lay : layouts) {
        auto residuals = compute_edge_residuals(lay, recig);
        REQUIRE(!residuals.empty());
        CHECK(median_abs(residuals) <= 1.0);
        CHECK(max_abs(residuals) <= 10.0);
    }
}

TEST_CASE("spec: noisy (approximate) intervals still produce low median residuals") {
    std::mt19937 rng(1337);

    // Chain of alternating orientations with jittered endpoints.
    Frag A{1, 1200, 0, true};
    Frag B{2, 1200, 500, false};
    Frag C{3, 1200, 900, true};
    Frag D{4, 1200, 1400, false};
    Frag E{5, 1200, 1750, true};

    ll::recig_graph recig;

    // Add multiple constraints per pair to emulate multiple approximate hits.
    for (int k = 0; k < 4; ++k) {
        add_reciprocal_overlap_interval(recig, A, B, 600 + 10 * k, 850 + 10 * k, /*jitter*/ 8, &rng);
        add_reciprocal_overlap_interval(recig, B, C, 950 + 7 * k, 1150 + 7 * k, /*jitter*/ 8, &rng);
        add_reciprocal_overlap_interval(recig, C, D, 1450 + 11 * k, 1650 + 11 * k, /*jitter*/ 8, &rng);
        add_reciprocal_overlap_interval(recig, D, E, 1800 + 9 * k, 2050 + 9 * k, /*jitter*/ 8, &rng);
    }

    auto layouts = ll::connected_comp_plus_layout(recig);
    REQUIRE(layouts.size() == 1);
    const auto& lay = layouts[0];
    REQUIRE(lay.size() == 5);

    auto residuals = compute_edge_residuals(lay, recig);
    REQUIRE(residuals.size() >= 16);

    // With jitter, some residuals can be non-zero; median should remain small if layout is sensible.
    CHECK(median_abs(residuals) <= 15.0);

    // Max can be larger due to jitter/outliers; keep it bounded.
    CHECK(max_abs(residuals) <= 80.0);
}

TEST_CASE("spec: inconsistent constraints (odd cycle) should be detectable (currently not handled)") {
    // Three nodes in a triangle where every edge is assumed to be an orientation-flipping (RC) match
    // is impossible to satisfy consistently (odd cycle in a bipartite/orientation-flip graph).
    //
    // Since the current layout API returns only (pos, orient) assignments with no explicit "conflict"
    // reporting, this test asserts a weaker but still spec-based condition:
    //   the output must exhibit *some* sign of inconsistency, either
    //   - at least one edge violates the flip-orientation constraint, OR
    //   - the interval alignment residuals are non-trivially large.
    //
    // This avoids assuming in advance which failure mode the implementation will produce.
    Frag A{1, 800, 0, true};
    Frag B{2, 800, 300, false};
    Frag C{3, 800, 600, true};

    ll::recig_graph recig;

    // Create a triangle of reciprocal edges.
    // Even if each pair overlap is "plausible" in isolation, globally this is contradictory under
    // the "every edge flips orientation" model.
    add_reciprocal_overlap_interval(recig, A, B, 350, 650);
    add_reciprocal_overlap_interval(recig, B, C, 650, 750);

    // Force an additional A<->C constraint as if it were also an RC match.
    // We generate it by temporarily flipping C's truth orientation to satisfy the generator's
    // "opposite orientations" requirement; the contradiction is in the graph structure itself.
    Frag C_flip = C;
    C_flip.orient_plus = false;
    add_reciprocal_overlap_interval(recig, A, C_flip, 650, 750);

    auto layouts = ll::connected_comp_plus_layout(recig);
    REQUIRE(layouts.size() == 1);
    const auto& lay = layouts[0];
    REQUIRE(lay.size() == 3);

    // Check whether all edges respect the mandatory flip constraint.
    bool has_flip_violation = false;
    for (const auto& e : recig) {
        u64 u = e.first;
        u64 v = e.second.node;
        if (!lay.contains(u) || !lay.contains(v)) continue;
        if (lay.at(u).second == lay.at(v).second) {
            has_flip_violation = true;
            break;
        }
    }

    auto residuals = compute_edge_residuals(lay, recig);
    REQUIRE(!residuals.empty());

    // If orientations happen to alternate on every checked edge anyway (e.g., due to traversal order),
    // the contradiction should surface as larger residuals.
    double med = median_abs(residuals);
    double mx = max_abs(residuals);

    const bool has_signal = has_flip_violation || (med >= 5.0) || (mx >= 20.0);
    CHECK_MESSAGE(has_signal,
                  "Expected contradiction signal in odd-cycle case; "
                  "has_flip_violation=",
                  has_flip_violation,
                  " median_abs_residual=",
                  med,
                  " max_abs_residual=",
                  mx);
}