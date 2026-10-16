# Prodigal (the-miint fork)

## What this repo is

A fork of [hyattpd/Prodigal](https://github.com/hyattpd/Prodigal), the
prokaryotic gene finder. It adds a reentrant, embeddable C library
(`prodigal.h` / `prodigal_api.c`) on top of the unchanged upstream
algorithm and keeps the upstream CLI (`main.c`) working.

- **Active branch:** `v2.6.4-miint`. Remotes: `origin` = the-miint/Prodigal,
  `upstream` = hyattpd/Prodigal.
- **Primary consumer:** [GPL-boundary](../GPL-boundary). It vendors this repo
  as the submodule `ext/prodigal`, compiles the library sources listed in its
  `build.rs::build_prodigal()`, and binds the API through hand-written
  `#[repr(C)]` FFI in `src/tools/prodigal.rs`. GPL-boundary in turn serves
  the BSD-licensed DuckDB extension miint (`../duckdb-miint`).
- **Lineage:** this fork is based on upstream `GoogleImport` HEAD `c1e2d36`
  (2023), **not** the `v2.6.3` release tag. It includes upstream's
  post-release fixes `b34b63b` and `7b8614d`, so results differ from the
  `v2.6.3` release on metagenomic contigs shorter than 3 kbp. See
  `docs/upstream-behavior.md`.

## Correctness contract (read this first)

1. **CLI ≡ upstream.** The CLI's output must be byte-identical to an upstream
   `c1e2d36` build for every output type (`-f gbk|gff|sco`, `-a`, `-d`, `-s`,
   `-t`), apart from the version string. Verified 2026-10-16 on
   `anthus_aco.fas` and E. coli K-12 MG1655, in both modes.
2. **API ≡ CLI.** For the same input and options, the library must produce:
   - the same genes, coordinates, strands, partial flags and start types;
   - the same RBS motif and spacer, and the same `gc_cont`;
   - the same scores the CLI prints (`conf`, `score`, `cscore`, `sscore`,
     `rscore`, `uscore`, `tscore`), equal up to the CLI's print rounding.

   This holds in both single-genome and metagenomic mode.
3. **Behaviour changes are deliberate or absent.** Algorithmic output may
   diverge from upstream only by explicit decision. Record each one in
   `docs/upstream-behavior.md` and `CHANGES`, and regenerate the ground truth
   when you make it. Upstream quirks that look like bugs are listed in that
   file; do not "fix" them silently.

As of `dd88690`, (1) holds and (2) does not. The CLI runs its own copy of the
gene-finding pipeline, so CI, which only exercises the CLI, never tests the
API path. See "Two pipelines" in `docs/internals.md`.

## Build and test

```bash
make                    # CLI binary ./prodigal (links zlib; gzip input support)
make libprodigal.a      # static library, built with -DPRODIGAL_NO_MAIN
make test               # builds ./test_api and runs it; run from the repo root
make distclean          # removes objects, binaries and libs
```

- `test_api` opens `anthus_aco.fas` and `testdata/ground_truth/*` by relative
  path, so run it from the repo root. A missing fixture prints `[SKIP]` and
  counts as a pass.
- If a result looks wrong, start with `make distclean && make test`. Stale
  objects in the working tree once made `test_version_constants` fail; a
  clean build of the same commit passes 40/40.
- `training.c` is about 25 MB of generated model data (50 metagenomic bins),
  so a full rebuild takes about a minute.
- The CLI regression checks are in `docs/testing.md`. Run them whenever you
  touch anything outside `prodigal_api.c` / `prodigal.h` / `test_api.c`.

### CI caveats

- `.github/workflows/ci.yml` only triggers on `library-api`, `GoogleImport`
  and `master`. **The active branch `v2.6.4-miint` is not CI-tested.**
- The ground-truth files embed `Prodigal.v2.6.3`, but the CLI has printed
  `2.6.4` since `dd88690`. The CI `diff` steps therefore fail on the current
  tree, even though the outputs otherwise match.

## Layout

| File | Role |
|---|---|
| `prodigal.h` | Public API: config, SOA/AOS output structs, stats, error codes. Heavily commented. |
| `prodigal_api.c` | API implementation: sequence encoding, training, both gene-finding modes, SOA extraction. |
| `prodigal_internal.h` | `struct prodigal_ctx`, `pdg_alloc`/`pdg_free`, the `PDG_LOG`/`PDG_SET_ERROR` macros. |
| `main.c` | Upstream CLI, adapted to hold its state in a `prodigal_ctx_t`. Whole file sits behind `#ifndef PRODIGAL_NO_MAIN`. |
| `node.c/.h` | Start/stop node creation, coding/RBS/upstream scoring, start training, `intergenic_mod`, `-s` start file. |
| `dprog.c/.h` | Dynamic programming over nodes and `eliminate_bad_genes`. |
| `gene.c/.h` | `add_genes`, `tweak_final_starts`, `record_gene_data` (GFF attribute strings), CLI printers, translations, `calculate_confidence`. |
| `sequence.c/.h` | 2-bit encoding helpers, codon tests, `amino()` translation, the CLI's FASTA/GenBank readers, SD motif scanners. |
| `training.c/.h` | `struct _training`, training-file I/O, 50 hard-coded `initialize_metagenome_N` models (generated; never hand-edit). |
| `metagenomic.c/.h` | `initialize_metagenomic_bins` and the bin descriptions. The header declares functions that don't exist (`score_edges` etc.). |
| `bitmap.c/.h` | Bit get/set/clear/toggle (`test`, `set`, …; unprefixed global names). |
| `fptr.h` | zlib or stdio file abstraction, CLI only. |
| `test_api.c` | C test suite for the library (custom assert macros). |
| `testdata/ground_truth/` | Reference outputs from an upstream build on `anthus_aco.fas`. |
| `docs/` | Internal docs: `internals.md`, `testing.md`, `upstream-behavior.md`. |

## Conventions

- **Language:** C99, compiled with `-pedantic -Wall -O3`. Match the
  surrounding style: the upstream core uses 2-space K&R, while
  `prodigal_api.c`, `prodigal.h` and `test_api.c` use 4-space.
- **No I/O or exits in library code.** No `printf` or `fprintf` to stdout or
  stderr, and no `exit` / `abort` / `assert`, on any path reachable from the
  API. Put CLI-only output behind `#ifndef PRODIGAL_NO_MAIN`. Library builds
  (Makefile `.lib.o`, WASM CI, GPL-boundary) define `PRODIGAL_NO_MAIN`.
  GPL-boundary's stdout is its JSON control channel, so a stray `printf`
  corrupts the protocol.
- **No mutable global or static state.** GPL-boundary runs many contexts at
  once on separate threads in one process; TSan was clean with 4 concurrent
  contexts on 2026-10-16. Static tables must be `const`.
- **Errors:** return `PRODIGAL_ERR_*` and set `PDG_SET_ERROR(ctx, …)` with a
  message specific enough to act on. A context must stay reusable after any
  error.
- **Memory:** working buffers go through `pdg_alloc` / `pdg_free`, the
  optional custom allocator. Output structs and the context itself always use
  system `malloc`, so the caller can free outputs without a context. Check
  every allocation.
- **ABI:**
  - `prodigal_config_t` is append-only and keeps `struct_size` as its first
    field.
  - Any change to a public struct (`prodigal_config_t`,
    `prodigal_genes_soa_t`, `prodigal_stats_t`, …) must be mirrored in
    GPL-boundary's `#[repr(C)]` bindings in the same rollout.
  - GPL-boundary only size-checks the config struct, so field order needs a
    manual review.
- **Adding a `.c` file** means updating, together:
  - `CORE_SOURCES` in the `Makefile`;
  - the `CORE=` lists in `ci.yml` (WASM job);
  - GPL-boundary's `build_prodigal()` file list and its `rerun-if-changed`
    lines.
- **Stack:** library code runs on Rust worker threads with 2 MiB stacks.
  `train_starts_nonsd` already uses 1.31 MB of stack, and each
  `initialize_metagenome_N` uses about 558 KB. Don't add large stack arrays;
  put them on the heap.
- **Symbols:** core functions have short global names (`set`, `test`,
  `dprog`, `amino`, …). Today they don't collide with the other libraries
  linked into GPL-boundary. `-fvisibility=hidden` does not protect against
  static link-time collisions, so new non-static functions should carry a
  `prodigal_`/`pdg_` prefix or be `static`.

## Git

- Never `git add -A`, `git add .` or `git commit -a`. Stage explicit paths.
  Build artifacts are gitignored, but the tree also holds `.claude/` and
  other local files.
- Regenerate ground truth only from an upstream `c1e2d36` build (see
  `docs/testing.md`). The one exception is a deliberate behaviour change,
  which must be documented.

## Shipping to GPL-boundary

1. Commit and push to `origin v2.6.4-miint`.
2. In GPL-boundary, bump the `ext/prodigal` submodule to the new commit.
   Never edit `ext/prodigal` in place; all fork work happens here.
3. Update GPL-boundary's FFI structs, adapter, `build.rs` source list, and its
   `--describe`/`schema_version` surfaces to match, in the same GPL-boundary
   change.
