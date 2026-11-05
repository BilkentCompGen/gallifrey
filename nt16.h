#ifndef _GALLIFREY_NT16_H_
#define _GALLIFREY_NT16_H_

#include <cstddef>
#include <cstdint>

// The 4-bit nucleotide code of the BAM SEQ field (SAM spec 4.2.3):
// "=ACMGRSVTWYHKDBN" -> 0..15, two bases per byte, high nibble first. The
// tables carry htslib's names and values, so packed sequences are the same
// whichever reader produced them.

// ASCII -> nt16 code. Either case; U packs as T; '0'..'3' as A, C, G, T;
// anything else as N.
inline constexpr uint8_t seq_nt16_table[256] = {
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    1,  2,  4,  8,  15, 15, 15, 15, 15, 15, 15, 15, 15, 0,  15, 15,  //
    15, 1,  14, 2,  13, 15, 15, 4,  11, 15, 15, 12, 15, 3,  15, 15,  //
    15, 15, 5,  6,  8,  8,  7,  9,  15, 10, 15, 15, 15, 15, 15, 15,  //
    15, 1,  14, 2,  13, 15, 15, 4,  11, 15, 15, 12, 15, 3,  15, 15,  //
    15, 15, 5,  6,  8,  8,  7,  9,  15, 10, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,  //
};

// nt16 code -> letter.
inline constexpr char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

// nt16 code -> 2-bit base (A=0, C=1, G=2, T=3); 4 for anything ambiguous.
inline constexpr int seq_nt16_int[16] = {4, 0, 1, 4, 2, 4, 4, 4,
                                         3, 4, 4, 4, 4, 4, 4, 4};

// nt16 code of base `i` of a packed sequence: the high nibble for even `i`.
constexpr uint8_t bam_seqi(const uint8_t *s, size_t i) {
    const size_t shift = 4 * (~i & 1);
    return static_cast<uint8_t>((s[i >> 1] >> shift) & 0xf);
}

#endif
