#pragma once
#ifndef ALIGN_H
#define ALIGN_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <gtl/vector.hpp>

#include "ksw2.h"

// Dual-affine scoring; defaults are minimap2's map-hifi preset (-A1 -B4
// -O6,26 -E2,1 -z400). A gap of length l costs min(o1 + l * e1,
// o2 + l * e2).
struct align_scoring {
    int match = 1;
    int mismatch = 4;
    int ambig = 1;  // either base is not A/C/G/T
    int o1 = 6, e1 = 2;
    int o2 = 26, e2 = 1;
    int zdrop = 400;  // extension stops this far below its best score
};

// BAM-encoded CIGAR (len << 4 | op), adjacent equal ops merged.
enum : uint8_t { CIGAR_M = 0, CIGAR_I = 1, CIGAR_D = 2 };
struct cigar {
    gtl::vector<uint32_t> ops;
    void push(uint8_t op, uint32_t len) {
        if (len == 0) return;
        if (!ops.empty() && (ops.back() & 0xF) == op) {
            ops.back() += len << 4;
        } else {
            ops.push_back(len << 4 | op);
        }
    }
    void clear() { ops.clear(); }
};

// Edit statistics of an alignment, as minimap2 reports them.
struct align_stats {
    int64_t score = 0;
    uint32_t matches = 0;     // aligned identical bases (PAF column 10)
    uint32_t mismatches = 0;  // aligned differing bases
    uint32_t ins = 0, del = 0;
    uint32_t nm() const { return mismatches + ins + del; }  // NM:i
};

// ksw2's SIMD dual-affine kernel (ksw_extd2_sse) with a reusable result
// buffer. Sequences are base codes: 0-3 for A/C/G/T, 4 for anything else.
class aligner {
   public:
    explicit aligner(align_scoring sc = {}) : sc_{sc}, mat_{score_matrix(sc)} {
        std::memset(&ez_, 0, sizeof(ez_));
    }
    aligner(const aligner &) = delete;
    aligner &operator=(const aligner &) = delete;
    ~aligner() { free(ez_.cigar); }

    // Statistics and score of cg aligning q to r from their starts.
    align_stats score_cigar(const cigar &cg, const uint8_t *q,
                            const uint8_t *r) const {
        align_stats st;
        size_t i = 0, j = 0;
        for (uint32_t x : cg.ops) {
            const uint32_t len = x >> 4;
            switch (x & 0xF) {
                case CIGAR_M:
                    for (uint32_t k = 0; k < len;) {
                        const uint32_t same =
                            same_prefix(q + i + k, r + j + k, len - k);
                        st.matches += same;
                        st.score += int64_t(same) * sc_.match;
                        if ((k += same) == len) break;
                        st.score += mat_[q[i + k] * 5 + r[j + k]];
                        ++st.mismatches;
                        ++k;
                    }
                    i += len;
                    j += len;
                    break;
                case CIGAR_I:
                    st.ins += len;
                    st.score -= gap_cost(len);
                    i += len;
                    break;
                case CIGAR_D:
                    st.del += len;
                    st.score -= gap_cost(len);
                    j += len;
                    break;
            }
        }
        return st;
    }

    // Extension from (0, 0): the best-scoring alignment of a prefix of q to
    // a prefix of r (z-drop stops it inside junk). Appends its CIGAR and
    // returns the consumed query and reference lengths.
    void extend(const uint8_t *q, int n, const uint8_t *r, int m, int band,
                cigar &out, int &q_used, int &r_used) {
        q_used = r_used = 0;
        if (n == 0 || m == 0) return;
        ksw_extd2_sse(nullptr, n, q, m, r, 5, mat_.data(), sc_.o1, sc_.e1,
                      sc_.o2, sc_.e2, band, sc_.zdrop, 0, KSW_EZ_EXTZ_ONLY,
                      &ez_);
        if (ez_.max_q < 0) return;  // nothing scores above zero
        q_used = ez_.max_q + 1;
        r_used = ez_.max_t + 1;
        for (int k = 0; k < ez_.n_cigar; ++k) {
            out.push(static_cast<uint8_t>(ez_.cigar[k] & 0xF),
                     ez_.cigar[k] >> 4);
        }
    }

   private:
    // Length of the leading stretch where q and r agree on A/C/G/T bases,
    // compared a word at a time (code 4, N, never counts as agreeing).
    static size_t same_prefix(const uint8_t *q, const uint8_t *r, size_t n) {
        constexpr uint64_t kN = 0x0404040404040404ULL;
        size_t k = 0;
        for (; k + 8 <= n; k += 8) {
            uint64_t a, b;
            std::memcpy(&a, q + k, 8);
            std::memcpy(&b, r + k, 8);
            if (const uint64_t d = (a ^ b) | (a & kN)) {
                return k + (__builtin_ctzll(d) >> 3);
            }
        }
        while (k < n && q[k] == r[k] && q[k] < 4) ++k;
        return k;
    }

    // Score of query base (row) against reference base (column), both
    // A, C, G, T, then anything else.
    static std::array<int8_t, 25> score_matrix(const align_scoring &sc) {
        const int8_t m = sc.match, x = -sc.mismatch, a = -sc.ambig;
        // clang-format off
        return {m, x, x, x, a,
                x, m, x, x, a,
                x, x, m, x, a,
                x, x, x, m, a,
                a, a, a, a, a};
        // clang-format on
    }

    int64_t gap_cost(int l) const {
        return std::min<int64_t>(sc_.o1 + int64_t(l) * sc_.e1,
                                 sc_.o2 + int64_t(l) * sc_.e2);
    }

    align_scoring sc_;
    std::array<int8_t, 25> mat_;
    ksw_extz_t ez_;
};

#endif
