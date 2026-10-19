# Testing Prodigal

## Oracles

Correctness is judged against **upstream Prodigal at `c1e2d36`** (upstream's
`GoogleImport` HEAD, which this fork is based on). See
`upstream-behavior.md` for why that commit and not the `v2.6.3` tag.

- **An upstream build**, made by `make oracle` (below); it reports `V2.6.3`.
  The untracked `prodigal_native` in the maintainer's checkout is
  byte-identical to an unpatched build. The real `v2.6.3` tag is *not*: it
  differs on the E. coli 2 kbp fragments.
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

- **NCBI's genetic code table**, `testdata/ncbi/gc.prt` (version 4.6,
  verbatim; provenance in its `README.md`), for translation tables. Here
  the oracle is NCBI, not upstream: a table is supported only if Prodigal
  uses it as the file defines it, and upstream's tables that don't match
  are rejected (`genetic-codes.md`).

## CLI regression (CLI ≡ upstream)

Run `make regress` (`tests/cli_regress.sh`) after any change outside the API
files; usage, checks and environment are in the script header. Each failed
check leaves a `.diff` in `tests/out/regress`. Its last checks cover the
CLI's deviations from upstream (`upstream-behavior.md`, "Fork deviations"):
empty records (FASTA and EMBL), a closed standard input, and the rejected
translation tables. Each of these (`expect_exit`) checks an exit code and
a message:
- `-g` with a sample of rejected tables, also in meta mode: exit 15;
- `-g` with a table that is not a number, and a final `-G` or `-d`: exit
  15;
- `-t` with a copy of `ref_train.bin` patched to table 24: exit 4.

A last check runs `-g 32`, which upstream rejects, and expects genes called
with table 32 that differ from table 11's.

`make compare`
runs `manymasks.fa` for the mask cap (its single mode cases call genes, the
meta ones none); no corpus contig has 30,000 genes, so the library tests
check the gene array's growth.

**ARM64 and floating-point contraction.** ARM64 compilers fuse multiply-adds
by default, which changes the last bits of some scores: on 2026-10-16 a
default Linux ARM64 build (gcc or clang) printed one start score in
`meta.starts` as -55.59 instead of -55.58. With `-ffp-contract=off` it
matched the x86-64 ground truth byte for byte in every output, training file
included, so the Makefile builds with that flag (a no-op on x86-64, where the
default build emits no FMA). `make oracle` builds the oracle the same way.

**Training file on ARM64 macOS.** Even so, the training file an ARM64 macOS
build writes differs from `ref_train.bin` at the byte level (first seen in
`d48c935`; also with `-ffp-contract=off`), which points to Apple's libm or
toolchain. Prediction from the committed `ref_train.bin` is unaffected. On
ARM64 macOS the script therefore checks that training is deterministic
instead (`--binary` forces the byte compare). Intel macOS has not been
checked.

**Windows.** Text-mode `fopen` writes CRLF, so only there a line-final CR is
stripped. Mid-line CRs are real output: `anthus_aco.fas` has CRLF line
endings, and upstream keeps the header's CR inside `seqhdr="..."`.

On 2026-10-16 every check passed at `dd88690`.

## CLI ≡ upstream across flags and inputs (`make oracle`)

`tests/oracle_diff.sh` runs the fork's CLI and an upstream `c1e2d36` build,
made on the same machine with the same compiler, on a flag matrix × corpus,
and compares exit status, stderr and every output file. How it builds and
caches the oracle, the matrix, and the pass rules are in the script; the
corpus is in `tests/make_corpus.sh`. Failed cases leave `.diff` files under
`tests/out/oracle/cases/`. CI runs it in every Linux (x86-64 and ARM64, gcc
and clang) and macOS job, and in the `sanitizers` job against a CLI built
with ASan + UBSan. Build `./prodigal` with the same `CC` you pass to
`make oracle`.

- The oracle carries `tests/oracle-c1e2d36.patch`, which defines upstream's
  `nodes[-1]` read so the oracle runs on macOS (`upstream-behavior.md`).
  `ORACLE_PATCH=` builds it unpatched, for checking that rule on Linux.
- Every case must succeed on the oracle, except the deliberate usage errors;
  agreeing on an error is not agreement on output.
- Almost every case runs with `-q`. A few anthus cases run without it, with
  the banner normalised, and their stderr (training progress included) must
  match too.
- **`ECOLI=<fasta>`** adds E. coli K-12 MG1655 (`NC_000913.3`, 4.64 Mbp; the
  maintainer has it at `~/dev/G000005845.fasta`) and its 2,320 2 kbp
  fragments: genome-scale training, and the `slen < 3000` metagenomic
  penalties across many contigs. Upstream calls 4,319 genes in both modes on
  the genome and 6,039 (meta) on the fragments.
- The inputs upstream mishandles are left out, and `make compare` runs them
  (`fork_corpus` in `tests/lib.sh`): `lowcomplex.fa` (`CAT` and `ATG`
  repeats with more nodes than the node array is first sized for; upstream
  overflows it) and `manymasks.fa` (5,100 masked regions with `-m`;
  upstream exits).

On 2026-10-16: 150/150 cases identical with `ECOLI` set (about 50 s at
`JOBS=24`, including the oracle build), 98/98 without, against both the
patched and the unpatched oracle.

## Library tests (`make test`)

`test_api.c` is a self-contained runner; run it from the repo root. A
failing assert `return`s from that test, and the exit code is non-zero
unless every test passed. A missing fixture is a failure.

**Golden tests (phase 14): API ≡ CLI ground truth.** The library's output,
in the form of `tests/api_gff.h`, is compared line by line with the anthus
references, every gene with its protein (the tests' contexts set
`want_proteins`; the `.proteins` files are the CLI's `-a` output): meta
mode with `ref_meta.gff` and `.proteins`; training through the API with
`ref_train.bin` byte for byte and its predictions with `ref_single.gff`
and `.proteins`, also when the context was trained before or had another
model loaded; prediction from `ref_train.bin` with `ref_trained.gff` (equal
to `ref_single.gff`) and `ref_single.proteins`, and so must a model trained
through the API, exported and loaded into a fresh context. Two more tests
check that a contig's output doesn't depend on what the context did before
(a second pass, each record on a fresh context, and `prodigal_find_genes`
twice per sequence). Anthus contigs carry about one gene each, so the meta
golden test sees little of metagenomic scoring; `make compare` covers more.
The low-complexity tests run an `ATG` repeat whose nodes don't fit the node
array: growing it must give what the grown array gives (Phase 15). Phase 16
checks the other limits: the gene array grows, more than 5,000 masks are
kept, and the length limits sit at the CLI's boundary. Every comparison
also checks each protein's layout (`api_proteins_ok` in `tests/api_gff.h`:
offsets from 0, a letter per codon, `M` and `*`). Phase 17 checks the
letters of genes made by hand (no `M` at an Edge start, `X`, the table) and
that there are no proteins without `want_proteins`.

**Genetic codes (phase 18).** `test_genetic_codes_ncbi` parses
`testdata/ncbi/gc.prt`, taking the codon order from its `Base1`-`Base3`
rows. For each table number from -1 to 99, the library must support the
table if and only if Prodigal's own code meets the rules of
`genetic-codes.md`:
- `is_stop`, `is_start` and `amino` on all 64 codons. Of the start codons,
  only ATG, GTG and TTG count, the three Prodigal models.
- `pdg_no_stop`, the ORF-length score's stop probability, at three GC
  contents.
- `is_stop` on the spacer that `prodigal_set_training_sequences` writes
  between contigs, in every frame.
- `is_stop` and `is_start` on the 61 codons that hold an N.

A codon that is a stop or sense by context rules a table out. The outcome
is pinned too, as lists of supported and rejected NCBI tables, and every
metagenomic bin's table must be supported. `test_genetic_code_entry_points`
checks that `config_validate`, `create`, `set_translation_table`,
`training_info` and `load_training` reject the same tables with the same
whole message, and change nothing when they do. `genetic-codes.md` has the
rules and the tables.

**ABI.** A struct field missing from the ABI table is a compile error
(`abi_complete_*` in `prodigal_api.c`). `test_abi_layout` checks the field
order, that the table has every struct typedef in `prodigal.h`, and that
the table's names hash to the value pinned for `PRODIGAL_ABI_VERSION`, so
changing them without a new version fails.

On ARM64 macOS, whose training writes different bytes, the training tests
compare only the model fields training copies or counts (`gc`,
`trans_table`, `st_wt`, `uses_sd`), which every platform checks, and the
predictions. `make compare` still checks the bytes there, against the CLI
of the same build.

Under WASM (CI) the tests need `-sNODERAWFS=1` to see the fixtures. CI
gives them a 1 MiB stack with `-sSTACK_OVERFLOW_CHECK=2` (emsdk pinned),
so every call chain the tests run must fit in 1 MiB; they need about
600 KB, mostly a metagenomic bin's 558 KB model, against WASM's 64 KiB
default. Without the check, an overflow shows up as an unrelated error
("null function or function signature mismatch").

## API ≡ CLI on larger inputs (`make compare`)

`tests/api_compare.sh` runs the CLI and `tests/api_dump` (the public API
only) with the same options over the corpus and a flag matrix, in meta mode,
single mode, and single mode from a CLI-written training file, and diffs the
two in the form of `tests/api_gff.h`, every gene with its protein (the
CLI's `-a` records, `api_dump --proteins`); every case must be identical.
In the single mode cases the model the API trains (`api_dump
--write-model`) must also equal the CLI's `-t` file byte for byte. `ECOLI`
adds E. coli and its fragments. For any other input,
`tests/api_dump -p MODE [FLAGS] --proteins FASTA` and
`tests/api_dump --from-cli -p MODE --faa cli.faa cli.gff` print the two
sides.

## Performance baselines

Wall time and peak RSS (MiB = 1024 KiB), CLI vs `tests/api_dump` (the
public API, printing the same GFF lines), median of 3 runs on 2026-10-17:
AMD Ryzen 9 9950X3D, Debian, gcc 14.2, `make` defaults (`-O3`), one thread,
output to a file. Measured at `6617375`, so the library no longer clears its
gene array per call. "From model" predicts from a training file the CLI
wrote for E. coli. At `6617375` the API's own single-mode training still
built a different model (more genes; fixed since), so the first row is only
roughly comparable.

| Input | Mode | CLI | API |
|---|---|---|---|
| E. coli K-12 (4.64 Mbp) | single, trained on input | 2.50 s, 209 MiB | 2.21 s, 212 MiB |
| E. coli K-12 | single, from model | 1.19 s, 81 MiB | 1.07 s, 141 MiB |
| E. coli K-12 | meta | 14.32 s, 136 MiB | 14.14 s, 181 MiB |
| E. coli 2 kbp fragments (2,320) | meta | 6.01 s, 88 MiB | 6.03 s, 92 MiB |
| E. coli 2 kbp fragments | single, from model | 0.68 s, 49 MiB | 0.65 s, 53 MiB |

The API's extra peak RSS came from zeroing the whole grown node array with
`memset`, which commits every page. Since `c572679` the array is grown with
`calloc` (`reserve_nodes`), and the API is within 3–5 MiB of the CLI (E.
coli single from model 84 MiB, meta 138 MiB, single trained 152 MiB).
Compare later measurements with this table on the same machine.

## Sanitizers, threads, fuzzing and the stack

The library has no mutable globals; keep it that way and keep these clean.
Each change of compiler or flags rebuilds every object (about a minute,
mostly training.c), and the next plain `make` rebuilds them again.

```bash
SAN="-fsanitize=address,undefined -fno-sanitize-recover=all -g -O1"
make test compare CFLAGS="$SAN"                  # ASan + UBSan
make threads CFLAGS="-fsanitize=thread -g -O1"   # TSan
make fuzz FUZZ_TIME=600                          # libFuzzer, 10 minutes (clang)
make stack                                       # frames over 1 MiB fail
```

- CI's `sanitizers` job (clang, aborting on any report) runs the library
  tests, the API vs CLI comparison and the oracle comparison under ASan +
  UBSan, `make threads` under TSan, and `make fuzz` for 60 s (steps in
  `ci.yml`). Some failures show only there, e.g. the AOS strings outliving
  the SOA.
- **Threads:** `tests/threads.c` (its header says what it runs). A static
  counter planted in `prodigal_find_genes` makes TSan fail it.
- **Fuzzing:** `tests/fuzz_find_genes.c` (input format and checks in its
  header; it checks the proteins too); `make fuzz` builds everything with
  `FUZZ_CFLAGS` (Makefile). New
  inputs collect in `tests/out/fuzz/corpus` (CI caches it, so it grows
  across runs) and failing ones in `tests/out/fuzz/` (CI uploads them).
  About 300–500 inputs/s under ASan, mostly metagenomic scoring across the
  50 bins (comparison tracing is off: it cost two thirds of the time).
- **Stack:** consumers run the library on 2 MiB threads. The library
  objects are built with `-fstack-usage`; `make stack` lists their largest
  frames and fails on any over 1 MiB or of unbounded size (a VLA or
  `alloca`). The largest are `initialize_metagenome_N`'s (training.c builds
  each bin's model on the stack; `static const` tables would remove those
  frames and most of training.c's compile time but add 16 MB to the
  binary, so they stay). Frame sizes don't add up call chains; the WASM job
  (1 MiB, checked) and `make threads` (2 MiB threads) check those on the
  paths the tests take.
- On 2026-10-16 ASan/UBSan were clean on lengths 1–199, all-N, all-A,
  all-G, lowercase, IUPAC codes, gapped input, 6,000 N runs, closed ends,
  single mode on short contigs, and repeated runs on one context (no leaks).

## Downstream (GPL-boundary)

After bumping `ext/prodigal` in `../GPL-boundary`, run `cargo test prodigal`
there; it uses real POSIX shared memory. Its adapter tests use the anthus
sequences. Any public struct change also needs that repo's
`#[repr(C)]` bindings updated, and its config ABI size test passing.
