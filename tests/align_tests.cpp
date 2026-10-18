// Doctests of the ksw2 aligner wrapper (align.h).
#include <random>
#include <string>

#include "align.h"
#include "ext/doctest.h"

namespace {

std::string random_dna(size_t n, uint32_t seed) {
    std::mt19937 gen(seed);
    std::string s(n, 'A');
    for (char &c : s) c = "ACGT"[gen() & 3];
    return s;
}

gtl::vector<uint8_t> codes_of(const std::string &s) {
    gtl::vector<uint8_t> c(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        switch (s[i]) {
            case 'A':
                c[i] = 0;
                break;
            case 'C':
                c[i] = 1;
                break;
            case 'G':
                c[i] = 2;
                break;
            default:
                c[i] = 3;
                break;
        }
    }
    return c;
}

// Walks a CIGAR over q and r from their starts: consumed lengths, and the
// matches and edits (mismatches + inserted + deleted bases) it implies.
struct replayed {
    uint32_t q = 0, r = 0, matches = 0, edits = 0;
};
replayed replay(const cigar &cg, const std::string &q, const std::string &r) {
    replayed o;
    for (uint32_t x : cg.ops) {
        const uint32_t n = x >> 4;
        switch (x & 0xF) {
            case CIGAR_M:
                for (uint32_t k = 0; k < n; ++k) {
                    if (q[o.q + k] == r[o.r + k]) {
                        ++o.matches;
                    } else {
                        ++o.edits;
                    }
                }
                o.q += n;
                o.r += n;
                break;
            case CIGAR_I:
                o.q += n;
                o.edits += n;
                break;
            case CIGAR_D:
                o.r += n;
                o.edits += n;
                break;
        }
    }
    return o;
}

// Applies `edits` edits spaced `step` apart from `first`: substitutions,
// single-base insertions and deletions, and homopolymer extensions in turn.
std::string mutate(std::string s, uint32_t first, uint32_t step,
                   uint32_t edits) {
    uint32_t pos = first;
    for (uint32_t e = 0; e < edits && pos < s.size(); ++e, pos += step) {
        switch (e % 4) {
            case 0:
                if (s[pos] == 'A')
                    s[pos] = 'C';
                else
                    s[pos] = 'A';
                break;
            case 1:
                if (s[pos] == 'G')
                    s.insert(s.begin() + pos, 'T');
                else
                    s.insert(s.begin() + pos, 'G');
                break;
            case 2:
                s.erase(s.begin() + pos);
                break;
            case 3:
                s.insert(s.begin() + pos, s[pos]);  // homopolymer +1
                break;
        }
    }
    return s;
}

}  // namespace

TEST_CASE("aligner::extend stops where the query turns into junk") {
    aligner dp;
    const std::string r = random_dna(2000, 5);
    const std::string q = r.substr(0, 800) + random_dna(400, 6);
    const auto qc = codes_of(q), rc = codes_of(r);
    cigar cg;
    int q_used = 0, r_used = 0;
    dp.extend(qc.data(), q.size(), rc.data(), 1300, 100, cg, q_used, r_used);
    CHECK(q_used >= 795);
    CHECK(q_used <= 805);
    CHECK(q_used == r_used);
    const replayed rp = replay(cg, q, r);
    CHECK(rp.q == uint32_t(q_used));
    CHECK(rp.r == uint32_t(r_used));
}

TEST_CASE("aligner::score_cigar scores an extension through edits") {
    aligner dp;
    for (uint32_t seed = 1; seed <= 20; ++seed) {
        CAPTURE(seed);
        const std::string r = random_dna(300 + seed * 37, seed);
        const std::string q = mutate(r, 13 + seed, 41, 6);
        const auto qc = codes_of(q), rc = codes_of(r);
        cigar cg;
        int q_used = 0, r_used = 0;
        dp.extend(qc.data(), q.size(), rc.data(), r.size(), 32, cg, q_used,
                  r_used);
        const align_stats st = dp.score_cigar(cg, qc.data(), rc.data());
        const replayed rp = replay(cg, q, r);
        CHECK(rp.q == uint32_t(q_used));
        CHECK(rp.r == uint32_t(r_used));
        CHECK(rp.matches == st.matches);
        CHECK(rp.edits == st.nm());
        CHECK(st.nm() <= 6);
        CHECK(q_used + 10 >= int(q.size()));
        // The reported score is the dual-affine score of the CIGAR.
        int64_t expect = int64_t(st.matches) - 4 * int64_t(st.mismatches);
        for (uint32_t x : cg.ops) {
            const int64_t l = x >> 4;
            if ((x & 0xF) != CIGAR_M) expect -= std::min(6 + 2 * l, 26 + l);
        }
        CHECK(st.score == expect);
    }
    // Empty query: nothing aligned.
    const auto rc = codes_of("ACGTTGCA");
    cigar cg;
    int q_used = 1, r_used = 1;
    dp.extend(rc.data(), 0, rc.data(), 8, 8, cg, q_used, r_used);
    CHECK(cg.ops.empty());
    CHECK(q_used == 0);
    CHECK(r_used == 0);
}
