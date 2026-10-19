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
   `-t`), apart from the version string and the deviations listed in
   `docs/upstream-behavior.md` (among them the rejected tables, item 3).
   Verified 2026-10-16 on `anthus_aco.fas` and E. coli K-12 MG1655, in
   both modes.
2. **API ≡ CLI.** For the same input and options, the library must produce:
   - the same genes, coordinates, strands, partial flags and start types;
   - the same RBS motif and spacer, and the same `gc_cont`;
   - the same scores the CLI prints (`conf`, `score`, `cscore`, `sscore`,
     `rscore`, `uscore`, `tscore`), equal up to the CLI's print rounding;
   - with `want_proteins`, the proteins of the CLI's `-a` file.

   This holds in both single-genome and metagenomic mode.
3. **Genetic codes are NCBI's, or rejected.** A translation table is
   accepted only if Prodigal uses it exactly as NCBI's `gc.prt` defines it
   (the rules: `docs/genetic-codes.md`). Upstream's tables that fall short
   are rejected, not fixed (decision of 2026-10-19). `test_api` checks
   every table against `testdata/ncbi/gc.prt`. Never accept a table that
   test does not pass, and change none of the code it checks without the
   decision process in `docs/genetic-codes.md`.
4. **Behaviour changes are deliberate or absent.** Algorithmic output may
   diverge from upstream only by explicit decision. Record each one in
   `docs/upstream-behavior.md` and `CHANGES`, and regenerate the ground truth
   when you make it. Upstream quirks that look like bugs are listed in that
   file; do not "fix" them silently.

The CLI and the API share one gene-finding pipeline ("One pipeline" in
`docs/internals.md`); `make compare` checks (2) on every CI run.

## Build and test

```bash
make                    # CLI binary ./prodigal (links zlib; gzip input support)
make libprodigal.a      # static library, built with -DPRODIGAL_NO_MAIN
make test               # builds ./test_api and runs it; run from the repo root
make compare            # library API vs the CLI on the same build (tests/api_compare.sh)
make regress            # CLI vs testdata/ground_truth, every output type, both modes
make oracle             # CLI vs an upstream c1e2d36 build: flag matrix x corpus
make pedantic           # every source as strict C99, no warnings
make threads            # contexts on four 2 MiB threads (TSan in CI)
make stack              # fails on a stack frame over 1 MiB or unbounded
make fuzz               # libFuzzer on the sequence path, with clang (docs/testing.md)
make distclean          # removes objects, binaries and libs
```

- `test_api` opens `anthus_aco.fas` and `testdata/ground_truth/*` by relative
  path, so run it from the repo root. A missing fixture is a failure.
- If a result looks wrong, start with `make distclean && make test`. Stale
  objects in the working tree once made `test_version_constants` fail (the
  Makefile now rebuilds when the compiler or flags change).
- `training.c` is about 25 MB of generated model data (50 metagenomic bins),
  so a full rebuild takes about a minute.
- Run `make regress` and `make oracle` whenever you touch anything outside
  `prodigal_api.c` / `prodigal.h` / `test_api.c`; set `ECOLI` for the
  genome-scale cases before changing scoring, training or memory sizing.
  See `docs/testing.md`.

### CI

`.github/workflows/ci.yml` runs, for pushes and PRs to `v2.6.4-miint`:
- `make test`, `make regress`: Linux x86-64 and ARM64 (gcc, clang), macOS,
  Windows;
- `make threads`, `make compare`, `make oracle`: Linux and macOS;
- `make pedantic`, `make stack`: Linux;
- `test_api` under WASM on a checked 1 MiB stack;
- the `sanitizers` job: `make test`, `make compare` and the oracle
  comparison under ASan + UBSan, `make threads` under TSan, a minute of
  `make fuzz` (some failures show only there).

## Layout

| File | Role |
|---|---|
| `prodigal.h` | Public API: config, SOA/AOS output structs, stats, error codes. Heavily commented. |
| `prodigal_api.c` | API implementation: sequence encoding, training, both gene-finding modes, SOA extraction. |
| `prodigal_internal.h` | `struct prodigal_ctx`, `pdg_alloc`/`pdg_calloc`/`pdg_free`, the `PDG_LOG`/`PDG_SET_ERROR` macros, and the `pdg_*` pipeline shared by the CLI and the API (defined in `prodigal_api.c`). |
| `main.c` | Upstream CLI: input, progress text and printing around the shared `pdg_*` pipeline, with its state in a `prodigal_ctx_t`. Whole file sits behind `#ifndef PRODIGAL_NO_MAIN`. |
| `node.c/.h` | Start/stop node creation, coding/RBS/upstream scoring, start training, `intergenic_mod`, `-s` start file. |
| `dprog.c/.h` | Dynamic programming over nodes and `eliminate_bad_genes`. |
| `gene.c/.h` | `add_genes`, `tweak_final_starts`, `record_gene_data` (GFF attribute strings), CLI printers, `pdg_protein_letter` (translation, shared by `-a` and the library), `calculate_confidence`. |
| `sequence.c/.h` | 2-bit base tests, codon tests, the supported translation tables (`pdg_trans_table_unsupported`), `rcom_seq`, `amino()` translation, the CLI's FASTA/GenBank/EMBL readers (which encode through `pdg_encode_*` in `prodigal_api.c`), SD motif scanners. |
| `training.c/.h` | `struct _training`, training-file I/O, 50 hard-coded `initialize_metagenome_N` models (generated; never hand-edit). |
| `metagenomic.c/.h` | `initialize_metagenomic_bins` and the bin descriptions. The header declares functions that don't exist (`score_edges` etc.). |
| `bitmap.c/.h` | Bit get/set/clear/toggle (`test`, `set`, …; unprefixed global names). |
| `fptr.h` | zlib or stdio file abstraction, CLI only. |
| `test_api.c` | C test suite for the library (custom assert macros). |
| `testdata/ground_truth/` | Reference outputs from an upstream build on `anthus_aco.fas`. |
| `testdata/ncbi/` | NCBI's genetic code table `gc.prt` (4.6), verbatim: the oracle for translation tables. |
| `docs/` | Internal docs: `internals.md`, `testing.md`, `upstream-behavior.md`, `genetic-codes.md`. |

## Conventions

- **Language:** C99, compiled with `-pedantic -Wall -O3`; `make pedantic`
  checks strict C99 (`-std=c99 -Werror`). Match the
  surrounding style: the upstream core uses 2-space K&R, while
  `prodigal_api.c`, `prodigal.h` and `test_api.c` use 4-space.
- **No I/O or exits in library code.** No `printf` or `fprintf` to stdout or
  stderr, and no `exit` / `abort` / `assert`, on any path reachable from the
  API. Put CLI-only output behind `#ifndef PRODIGAL_NO_MAIN`. Library builds
  (Makefile `.lib.o`, WASM CI, GPL-boundary) define `PRODIGAL_NO_MAIN`.
  GPL-boundary's stdout is its JSON control channel, so a stray `printf`
  corrupts the protocol.
- **No mutable global or static state.** GPL-boundary runs many contexts at
  once on separate threads in one process; `make threads` runs 4 under TSan
  in CI. Static tables must be `const`.
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
  - Such a change, a changed function signature or a changed public
    constant bumps `PRODIGAL_ABI_VERSION` (`prodigal.h`); then pin the new
    names hash in `test_abi_layout` (its failure prints it). A new field
    goes in its struct's list in `prodigal_api.c` (`CONFIG_FIELDS` and the
    like), in order: the build fails when a list misses a field, and
    `test_abi_layout` when the order is wrong or a struct typedef is
    missing.
  - GPL-boundary only size-checks the config struct today (its layout check
    against `prodigal_abi_layout()` is planned), so field order still needs
    a manual review there.
- **Adding a `.c` file** means updating, together:
  - `CORE_SOURCES` in the `Makefile`;
  - the `CORE=` lists in `ci.yml` (WASM job);
  - GPL-boundary's `build_prodigal()` file list and its `rerun-if-changed`
    lines.
- **Floating point:** build with `-ffp-contract=off` (the Makefile does).
  With fused multiply-adds, ARM64 output drifts from x86-64 in the last
  printed digit. GPL-boundary's `build.rs` must use the same flag.
- **Stack:** library code runs on Rust worker threads with 2 MiB stacks.
  Each `initialize_metagenome_N` uses about 558 KB. Don't add large stack
  arrays; put them on the heap. `make stack` fails on a frame over 1 MiB.
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
