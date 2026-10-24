# HiFi assembler

## Build

Requirements:

- CMake >= 3.20
- C++23 compiler
- `libdeflate`
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
| `Release` (default) | `-O3 -g1 -ftree-vectorize -DNDEBUG` (line-level debug info); link-time optimisation with GCC |
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

- BAM (mapped or unmapped; every record's sequence is used)
- FASTA and FASTQ, plain, gzipped or bgzipped

SAM and CRAM are not read directly; convert them first, for example with
`samtools view -b -o reads.bam reads.sam` (or `reads.cram`), or
`samtools fastq reads.cram > reads.fq`.

The format is detected from the file content, so any file name works. Read
names and quality values are ignored; only sequences are used. Different
formats can be mixed in one invocation.

An input can also be a pipe: `-` reads standard input, and a named pipe or
`/dev/stdin` works like a file, for example
`samtools view -b reads.cram | lasm hifi - -o assembly.gfa -l 5`.
Several pipes can be given, but each must already have its writer running
when lasm starts: lasm opens every input before it reads any, so pipes filled
one after the other would wait forever. `--multipass` reads its inputs once
per pass, so it needs regular files. Each input can be given only once.

```bash
lasm hifi <reads1> [reads2 ...] -o output.gfa -l <lcp_level> [--syncmer K/S[/T]]
```

## Basic Usage

Minimal example (default syncmers: k=11, open syncmer with s=4, t=6):

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
  -w 6 \
  --syncmer 11/4/6 \
  --iterations 2
```

## Main Parameters

`lasm hifi --help` lists these; `lasm hifi --help-all` shows every option,
including the advanced ones below.

Required:

- `input`: one or more BAM or FASTA/FASTQ(.gz) read files; `-` reads
  standard input
- `-o, --output`: output GFA path
- `-l, --level`: LCP/core depth (or `--multipass`)

Parallelism and input:

- `-t, --threads`: worker thread count
- `-@, --io-threads`: decompression thread count per BGZF input file (BAM and bgzipped FASTA/FASTQ; plain gzip is decompressed on one thread, uncompressed input needs none)
- `--file-reading-mode`: `sequential` or `parallel`

Assembly:

- `-w, --window`: cores per lcmer window, used for DBG construction and
  overlap detection (default 6)
- `--multipass LEVEL:WINDOW[,LEVEL:WINDOW...]`: multi-pass assembly, e.g.
  `4:8,5:6,5:10`; each item is the LCP level and window of one
  pass, replacing `-l` and `-w`. Each pass assembles the reads together
  with the previous pass's contigs; only the last pass writes the output (and
  the `--print-raw-dbg`/`--print-bi-dbg` dumps). Small windows bridge
  construction holes and large ones resolve repeats, so order the passes from
  small to large windows. Every pass reads the inputs again, so they must be
  files, not pipes or standard input
- `--ploidy`: ploidy of the sample, `1` for a haploid (or fully homozygous) genome (default) or `2` for a diploid one
- `--silent`: reduce progress logging (with `--multipass`, every pass logs its
  settings, the segments it uses and keeps, and its wall/CPU time)

K-mer selection:

A k-mer is selected by the position of its smallest s-mer; selection is
computed on the fly from s-mer hashes, with no per-k-mer table. K-mers and
s-mers count run-length runs, not bases.

- `--syncmer K/S[/T]`: syncmer k-mer length `K` (at most 16), s-mer length
  `S` (1 <= S < K) and, for an open syncmer, offset `T` of the smallest s-mer
  (0 <= T <= K - S). Without `/T` it is a closed syncmer (smallest s-mer at
  the first or last position). Default `11/4/6`
- `--random-seed`: seed of the s-mer hash

## Advanced Parameters

The defaults suit most runs.

### Bloom filter

- `--projected-element-count`: distinct lcmer keys the Bloom filter is sized
  for (default 100000000); construction warns when the filter is undersized
- `--false-positive-probability`: target false-positive probability
  (default 0.001)
- `--prime-reverse` / `--no-prime-reverse`: also parse each read's reverse
  complement and insert its lcmer keys into the Bloom filter only (no node,
  count or edge), so reads of the other strand are admitted at their first
  sighting (default on)

### Graph cleaning

- `--iterations`: number of graph-cleaning passes (default 6)
- `--max-tip-length`: maximum tip length in cores (default 16)
- `--max-bubble-length`: largest bubble to pop, in inner cores (default 1024); `0` disables bubble popping
- `--min-bubble-length-difference`: maximum length ratio of a bubble's removed nodes to its kept path (default 2)

### Layout

- `--drop-crumbs` / `--keep-crumbs`: drop crumb nodes (too short for overlap detection) that have a single predecessor and continue into a longer node, instead of emitting them as `C` segments (default on; halves the contig count at a small genome-fraction cost)
- `--strict-layout-labels` / `--no-strict-layout-labels`: accept a core-label
  match as a seed only if both cores have the same homopolymer-compressed
  sequence, so cores that merely share a label are not chained (default on)
- `--end-joins` / `--no-end-joins`: join two free segment ends that overlap as reverse complements (the two read strands of one locus stopped at different holes) into one segment; candidates need a unique partner, one-locus depth at both ends, end windows that are not tandem-periodic (collapsed satellite arrays), and an overlap that does not recur near either end (a longer-period array, where the overlap's position is ambiguous) (default on)

### Other

- `--print-raw-dbg`: write the pre-cleaning graph with `.raw.gfa` suffix
- `--print-bi-dbg`: write the graph before layout with `.bi.gfa` suffix
- `--retro-admission` / `--no-retro-admission`: remember each lcmer window's
  first sighting on a read's own strand until a later sighting on either
  strand admits the window, then add that sighting's count and edges as if it
  had been admitted at once. Costs memory for the windows still waiting
  (default on)
- `--link-reverse`: let each read's reverse complement add the missing links
  between lcmer windows that are already in the graph, so the two strands'
  graphs fill each other's gaps (default off)

## Output Files

Primary output:

- `<output>`: final merged layout/consensus graph in GFA format
- M[0-9]* nodes indicate layout nodes which are strand resolved
- C[0-9]* nodes are not strand resolved (will have duplicates)

Optional auxiliary outputs:

- `<output>.raw.gfa`: raw DBG before cleaning, if `--print-raw-dbg` is enabled
- `<output>.bi.gfa`: cleaned DBG before layout merging if `--print-bi-dbg` is enabled
