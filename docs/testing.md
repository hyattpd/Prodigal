# Testing Prodigal

## Oracles

Correctness is judged against **upstream Prodigal at `c1e2d36`** (upstream's
`GoogleImport` HEAD, which this fork is based on). See
`upstream-behavior.md` for why that commit and not the `v2.6.3` tag.

- **`prodigal_native`** (untracked, in the maintainer's checkout) is an
  upstream build. It reports `V2.6.3` but is byte-identical to a `c1e2d36`
  build: it agrees with the fork's CLI on 2,320 E. coli 2 kbp fragments,
  where the real `v2.6.3` tag differs. To rebuild it from scratch:

  ```bash
  mkdir -p <dir> && git archive c1e2d36 | tar -x -C <dir> && make -C <dir>
  ```

- **`testdata/ground_truth/`** holds upstream outputs on `anthus_aco.fas`.
  That file has 22 contigs of 960 bp each (21,142 bp). Its concatenation
  passes the 20 kbp training minimum but sits below the 100 kbp "ideal",
  hence the warning.

  | File | Command |
  |---|---|
  | `ref_meta.{gff,gbk,sco}` | `-p meta -f gff\|gbk\|sco` |
  | `ref_meta.proteins` / `.nucl` / `.starts` | `-p meta -a` / `-d` / `-s` |
  | `ref_single.{gff,gbk,sco,proteins,nucl,starts}` | `-p single` with the same output flags |
  | `ref_train.bin` | `-t ref_train.bin` (write a training file) |
  | `ref_trained.gff` | `-t ref_train.bin -f gff` (read it back) |

  `ref_single.gbk` and `ref_single.sco` were added on 2026-10-16 from a fresh
  `c1e2d36` build; on that day every other file was also re-checked as
  byte-identical to that build's output.

  The text files embed upstream's version string `Prodigal.v2.6.3`; the fork
  prints its own, so each side is normalised for its own version before
  diffing (`tests/lib.sh`). A script pointed at an upstream binary therefore
  fails instead of passing vacuously.

## CLI regression (CLI ≡ upstream)

Run `make regress` (`tests/cli_regress.sh`) after any change outside the API
files; usage, checks and environment are in the script header. Each failed
check leaves a `.diff` in `tests/out/regress`.

**Training file on ARM64.** `ref_train.bin` was written on x86-64 Linux.
ARM64 builds read it correctly (`struct _training` has the same layout on
every LP64 ABI), but the file they *write* differs at the byte level, first
seen on macOS in `d48c935`. The suspected cause is FMA contraction; CI's
report-only `determinism` job measures it with and without
`-ffp-contract=off`, on macOS and on Linux ARM64 with gcc and clang. Until
that settles, the script skips the byte compare on ARM64 (use `--binary` to
force it) and checks that training is deterministic instead. Prediction from
the committed `ref_train.bin` is still checked on every platform.

**Windows.** Text-mode `fopen` writes CRLF, so only there a line-final CR is
stripped. Mid-line CRs are real output: `anthus_aco.fas` has CRLF line
endings, and upstream keeps the header's CR inside `seqhdr="..."`.

On 2026-10-16 every check passed at `dd88690`.

### Larger inputs (not in the repo)

anthus never exercises genome-scale training, large node arrays, or many
short metagenomic contigs. Before shipping a change to scoring, training or
memory sizing, also compare against `prodigal_native` on:

- **E. coli K-12 MG1655** (`NC_000913.3`, 4.64 Mbp; the maintainer has it at
  `~/dev/G000005845.fasta`). Upstream calls 4,319 genes in both `-p single`
  and `-p meta`. Fork and native were identical for gff/faa/fna in both
  modes on 2026-10-16.
- **E. coli shredded into 2 kbp fragments** (2,320 contigs). Upstream meta
  calls 6,039 genes. This set exercises the `slen < 3000` metagenomic
  penalties and the meta bin selection across many contigs. It is also the
  set on which `c1e2d36` and the `v2.6.3` tag differ.
- Low-complexity and synthetic input (all-N, all-A, long `ATG` repeats,
  random 30–80% GC). Use these for memory-safety checks rather than output
  comparison.

## Library tests (`make test`)

`test_api.c` is a self-contained runner. A failing assert `return`s from
that test and leaves `tests_passed` short; the exit code is non-zero unless
every test passed. A missing fixture prints `[SKIP]` and **counts as a
pass**.

What it does *not* check today:

- **Meta output against the reference:** `test_meta_matches_reference`
  compares gene count, coordinates and strand against `ref_meta.gff`, but
  not partial flags, start type, RBS fields, `gc_cont` or any score.
- **API training against the reference:** `test_train_from_sequences` only
  checks that training succeeds; it never compares against `ref_train.bin`.
  It would have caught the missing `tinf.gc`.
- **API single mode end to end:** no test trains from sequences via the API,
  calls `find_genes` and compares with `ref_single.gff`.
  `test_find_genes_with_training` uses `ref_train.bin`, so it bypasses API
  training, and it checks only value ranges.
- **Proteins or nucleotides:** the API has no such output yet.

The goal is for the library tests to assert **API ≡ CLI** field by field
(CLI print precision for scores), in both modes, on anthus. It is cheap to
add, and the larger sets above can serve as an opt-in slow tier.

## Sanitizers and concurrency

The library has no mutable globals; keep it that way and keep these clean.

```bash
CORE="bitmap.c dprog.c gene.c metagenomic.c node.c sequence.c training.c prodigal_api.c"
# ASan + UBSan
gcc -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
    -DPRODIGAL_NO_MAIN $CORE test_api.c -o test_api_asan -lm && ./test_api_asan
# TSan: drive several contexts from separate pthreads with a small harness
gcc -O1 -g -fsanitize=thread -DPRODIGAL_NO_MAIN $CORE harness.c -o h_tsan -lm -lpthread
```

- On 2026-10-16 ASan/UBSan were clean on lengths 1–199, all-N, all-A,
  all-G, lowercase, IUPAC codes, gapped input, 6,000 N runs, closed ends,
  single mode on short contigs, and repeated runs on one context (no leaks).
- TSan was clean with 4 concurrent contexts, and results were deterministic
  across threads.
- The known ASan failure is the `add_nodes` overflow on long `ATG` repeats
  (see `internals.md`, "Limits and memory").
- **Stack:** use `-fstack-usage` to watch frame sizes. Consumers run the
  library on 2 MiB threads, and `train_starts_nonsd` is already 1.31 MB.

## Downstream (GPL-boundary)

After bumping `ext/prodigal` in `../GPL-boundary`, run `cargo test prodigal`
there; it uses real POSIX shared memory. Its adapter tests use the anthus
sequences. Any public struct change also needs that repo's
`#[repr(C)]` bindings updated, and its config ABI size test passing.
