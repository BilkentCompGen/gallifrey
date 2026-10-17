# HiFi assembler

## Build

Requirements:

- CMake >= 3.20
- C++23 compiler
- `htslib`
- `zlib`
- POSIX threads

Build (defaults to `Release`):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The main executable is:

```bash
./build/lasm
```

### Build types

Select a build type with `-DCMAKE_BUILD_TYPE=<type>`. Each writes to its own
build directory so you can keep them side by side.

| Type | Flags (summary) | 
|------|-----------------|
| `Release` (default) | `-O3 -g -ftree-vectorize -DNDEBUG` |
| `RelWithDebInfo` | `-O2 -g -pg -DNDEBUG` | 
| `Debug` | `-O0 -g -DDEBUG` + UBSan (`-fsanitize=undefined,null,bounds -fsanitize-trap=all`) | 
| `MinSizeRel` | `-Os -DNDEBUG` |

### Build options
- `LASM_NATIVE` — add `-march=native`
- `LASM_X86_64_V2` (default `ON`) — without `LASM_NATIVE`, x86-64 builds
  target `-march=x86-64-v2` (SSSE3/SSE4.2/POPCNT, any x86-64 CPU from about
  2009 on), so the SIMD paths are used. Turn it off only for older CPUs; the
  scalar fallbacks are several times slower. ARM64 builds always use NEON.

```bash
# Default build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Fastest local build for this machine
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLASM_NATIVE=ON
cmake --build build -j

# Profiling build with native tuning (produces gmon.out after a run)
cmake -S . -B build-prof -DCMAKE_BUILD_TYPE=RelWithDebInfo -DLASM_NATIVE=ON
cmake --build build-prof -j
```

## Input

`lasm hifi` accepts one or more read files. Supported formats:

- BAM, CRAM, SAM (unmapped reads are fine; CRAM with mapped reads needs the reference available to htslib)
- FASTA and FASTQ, plain or gzipped

The format is detected from the file content, so any file name works. Read
names and quality values are ignored; only sequences are used. Different
formats can be mixed in one invocation.

```bash
lasm hifi <reads1> [reads2 ...] -o output.gfa -l <lcp_level> [--syncmer-k K --submer-s S --submer-t T]
```

## Basic Usage

Minimal example (default syncmers: k=11, open syncmer with s=4, t=5):

```bash
build/lasm hifi reads.bam -o assembly.gfa -l 5
```

Example with explicit syncmer and thread settings:

```bash
build/lasm hifi reads.bam \
  -o assembly.gfa \
  -l 5 \
  -t 16 \
  -@ 8 \
  --lcmer-l 6 \
  --syncmer-k 11 \
  --submer-s 4 \
  --submer-t 5 \
  --iterations 2
```

## Main Parameters

Required:

- `input`: one or more BAM/CRAM/SAM or FASTA/FASTQ(.gz) read files
- `-o, --output`: output GFA path
- `-l, --level`: LCP/core depth

Core assembly controls:

- `--lcmer-l`: rolling `lcmer` window length (cores) used for DBG construction
  and overlap detection (default 6). `--dbg-lcmer-l` is accepted as the old
  name; `--olc-lcmer-l` is ignored.
- `--random-seed`: hash/random seed
- `--prime-reverse` / `--no-prime-reverse`: also parse each read's reverse
  complement and insert its lcmer keys into the Bloom filter only (no node,
  count or edge), so reads of the other strand are admitted at their first
  sighting (default on)
- `--strict-layout-labels` / `--no-strict-layout-labels`: during layout,
  accept a core-label match as a seed only if both cores have the same
  homopolymer-compressed sequence, so cores that merely share a label are not
  chained (default on)

Syncmer selection (Edgar 2021). A k-mer is selected by the position of its
smallest s-mer; selection is computed on the fly from s-mer hashes, with no
per-k-mer table. With `--rle` (default), k-mers and s-mers count run-length
runs, not bases.

- `--syncmer-k`: k-mer length, at most 16 (default 11)
- `--submer-s`: s-mer length, 1 <= s < k. Without `--submer-t` this is a
  closed syncmer (smallest s-mer at the first or last position)
- `--submer-t`: offset of the smallest s-mer for an open syncmer,
  0 <= t <= k - s (requires `--submer-s`)
- With neither `--submer-s` nor `--submer-t`, an open syncmer with s=4, t=5
  is used

Parallelism:

- `-t, --threads`: worker thread count
- `-@, --io-threads`: decompression thread count per input file (BAM/CRAM/SAM and gzipped FASTA/FASTQ; no effect on uncompressed FASTA/FASTQ)
- `--file-reading-mode`: `sequential` or `parallel`

Graph cleaning:

- `--iterations`: number of graph-cleaning passes (default 6)
- `--max-tip-length`: maximum tip length in cores
- `--ploidy`: copies of a locus, `1` (default) or `2`; a bubble whose nodes all stay within the depth of one locus (its measured depth plus three times its spread) holds alleles or read errors and is popped, a deeper one is a repeat and is kept. With `2`, short dead-end branches at one-locus forks (allele tips) are cut as well
- `--max-bubble-length`: largest bubble to pop, in inner cores (default 1024); `0` disables bubble popping
- `--min-bubble-length-difference`: maximum length ratio of a bubble's removed nodes to its kept path (default 2)

Layout:

- `--drop-crumbs` / `--keep-crumbs`: drop crumb nodes (too short for overlap detection) that have a single predecessor and continue into a longer node, instead of emitting them as `C` segments (default on; halves the contig count at a small genome-fraction cost)

Output and debugging:

- `--print-raw-dbg`: write the pre-cleaning graph with `.raw.gfa` suffix
- `--print-bi-dbg`: write the graph before layout with `.bi.gfa` suffix
- `--silent`: reduce progress logging

## Output Files

Primary output:

- `<output>`: final merged layout/consensus graph in GFA format
- M[0-9]* nodes indicate layout nodes which are strand resolved
- C[0-9]* nodes are not strand resolved (will have duplicates)

Optional auxiliary outputs:

- `<output>.raw.gfa`: raw DBG before cleaning, if `--print-raw-dbg` is enabled
- `<output>.bi.gfa`: cleaned DBG before layout merging if `--print-bi-dbg` is enabled
