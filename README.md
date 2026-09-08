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

`lasm hifi` accepts one or more positional HTS (prefer unmapped BAM, FASTX also works) inputs:

```bash
lasm hifi <input1.bam> [input2.bam ...] -o output.gfa -l <lcp_level> --density <kmer_density>
```

## Basic Usage

Minimal example:

```bash
build/lasm hifi reads.bam -o assembly.gfa -l 4 --density 0.05
```

Example that uses syncmer:

```bash
build/lasm hifi reads.bam \
  -o assembly.gfa \
  -l 3 \
  -t 8 \
  -@ 4 \
  --dbg-lcmer-l 8 \
  --syncmer-k 11 \
  --submer-s 4 \
  --submer-t 5 \
  --iterations 2 
```

## Main Parameters

Required:

- `input`: one or more input BAM files
- `-o, --output`: output GFA path
- `-l, --level`: LCP/core depth

Core assembly controls:

- `--dbg-lcmer-l`: rolling `lcmer` window length used for DBG construction
- `--syncmer-k`: syncmer k-mer length
- `--density`: density-based syncmer selection
- `--submer-s`: submer size for mincode/open-syncmer rules
- `--submer-t`: submer position for mincode/open-syncmer rules
- `--spaced-seed-mask`: optional spaced-seed mask
- `--random-seed`: hash/random seed

Parallelism:

- `-t, --threads`: worker thread count
- `-@, --io-threads`: HTSlib decompression thread count per file reader
- `--file-reading-mode`: `sequential` or `parallel`

Graph cleaning:

- `--iterations`: number of graph-cleaning passes
- `--max-tip-length`: maximum tip length in cores

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
