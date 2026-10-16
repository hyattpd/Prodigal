# Upstream lineage and preserved behaviour

## Lineage

```
v2.6.3 tag (004218f, 2016-02-11)        ← the "Prodigal 2.6.3" most users run (bioconda etc.)
  │  upstream GoogleImport branch:
  │   6e54227 / e83fb71  gzip input support (CLI only; zlib)
  │   f729050 …          Travis CI, README tweaks
  │   b34b63b (2021)     score_nodes: abs(x < 120) → abs(x) < 120   ← behaviour change
  │   7b8614d (2023)     SD scanners: initialise all of match[6]     ← removes UB
  │   7cb28c9 (2023)     SD scanners: refactor comparisons (no behaviour change intended)
c1e2d36 (2023-01-27)                     ← fork base = correctness oracle
  │  the-miint fork (v2.6.4-miint):
  │   6b81440  reentrant library API + test suite
  │   fea252a  CLI moved onto the library context (output byte-identical)
  │   a6a8206  prodigal.h API docs
  │   b0b2d78 / d48c935 / 0a9a9ef  GitHub Actions CI (Linux/macOS/Windows/WASM)
  │   8558797  PRODIGAL_NO_MAIN guards on library stderr writes
  │   4da852a / fe531b8 / dd88690  version → 2.6.4
```

Upstream has not released anything since `v2.6.3`, and its version string
still says 2.6.3 at `c1e2d36`. A "2.6.3" binary can therefore be either
lineage.

**Measured difference.** On E. coli K-12 shredded into 2,320 fragments of
2 kbp, run with `-p meta`:

| Build | Genes |
|---|---|
| `v2.6.3` tag | 6,025 |
| `c1e2d36` / this fork | 6,039 |

2,136 GFF lines differ between them. On `anthus_aco.fas` (960 bp contigs) the
two agree. Most of the difference comes from `b34b63b`. Before that fix,
`abs(ndx - stop_val < 120)` evaluated `abs()` of a boolean, so in meta mode
every internal gene on a contig under 3 kbp took the short-gene penalty, not
only genes under 120 bp. `7b8614d` removes reads of uninitialised stack in
`shine_dalgarno_exact`/`_mm`. Results with that UB present are
compiler-dependent.

**Decision on record:** the fork follows `c1e2d36`. Both upstream changes are
bug fixes, and pyrodigal carries the same ones. Results should be described
as matching Prodigal at `c1e2d36`, not "Prodigal 2.6.3", whenever that
distinction matters, e.g. comparisons with published 2.6.3 output on short
metagenomic contigs.

## Quirks preserved on purpose

All of these are upstream behaviour, the CLI reproduces them exactly, and
"fixing" any of them changes gene calls or reported values relative to the
oracle. Change one only by explicit decision: record it here and in
`CHANGES`, and regenerate the ground truth.

1. **`score_nodes` mutates node structure.** Starts within 3 bp of a contig
   end become edge nodes (`node.c:459-466`: `edge = 1`, scores reset).
   `reset_node_scores` does not restore `edge`, and in meta mode the nodes
   are rebuilt only when the translation table changes between consecutive
   bins. Later bins therefore score those nodes as pre-existing edge nodes,
   which takes a different branch of the edge-bonus logic.
2. **Meta GFF column 6 differs from the `score=` attribute.** Column 6 is
   printed from the recovery-pass nodes, which lack `intergenic_mod`. The
   attributes (`conf`, `score`, `sscore`, …) come from the in-loop values,
   which include it. In single mode the two agree. *For the library, the
   attribute values are the ones that count:* they are what drove the gene
   choice and what `conf=` reports.
3. **Ambiguous bases are encoded as C.** The encoder doesn't count them
   toward the contig GC fraction (they lower it, like A/T). But `is_gc()`
   sees a C, so they *do* count as GC in the per-ORF `gc_cont`, in the
   GC-frame plot, and in codon tests. Long N runs and all-N contigs can
   yield genes (an all-N 5 kbp contig gives one gene at 99.99% confidence).
   Only `-m` masking prevents genes across N runs of 50 bp or more.
4. **Non-letters are dropped, not counted.** Characters outside `'A'..'z'`
   don't advance the position, so coordinates shift on gapped or otherwise
   decorated input. The characters `[ \ ] ^ _ \`` count as ambiguous.
5. **A trailing N run is never masked.** A mask is recorded only when a
   non-N closes the run.
6. **`calc_most_gc_frame` leaves the last `slen % 3` positions at `-1`.**
   This slightly biases `record_gc_bias` during training.
7. **Confidence is clamped to `[50, 99.99]`.** Scores with
   `score/st_wt ≥ 41` report 99.99; values just under that can print as
   `100.00` at 2 decimals.
8. **Circular genomes are not handled.** Genes cannot span the origin.
9. **Meta mode ignores `-g`, `-n` and the start weight.** Each of the 50
   bins brings its own translation table (11 or 4), SD usage and start
   weight. The best bin's model is the one that applies to the contig's
   output and translations.
10. **Per-contig caps:**
    - more than 30,000 genes in one contig: truncated (CLI warns);
    - training input over `MAX_SEQ`: truncated (CLI warns);
    - more than 5,000 masks: the CLI exits.

## Upstream bugs present in the CLI too (fix, but preserve output)

These are defects in upstream code paths that both the CLI and the library
reach. Fixes must not change output on any input that does not trigger the
defect.

- **`nodes[-1]` read in the metagenomic loop (fixed).** When `dprog` finds
  no gene path it returns -1, and `main.c` then compared
  `nodes[-1].score` with the best score so far. That reads before the node
  array: ASan reports it on ordinary input (anthus), and it segfaults on
  current macOS. **Rule:** a bin without a gene path scores 0.0, the score of
  the empty path (`dprog` starts every node at 0 and only raises a score
  together with its traceback). A no-path bin therefore "wins" only when no
  earlier in-window bin has a path, and it adds no genes. That is also what
  upstream does on glibc, where the read lands in the allocator's header for
  the mmap'd node array and returns 0.0. `tests/oracle-c1e2d36.patch` applies
  the rule to the oracle. On 2026-10-16 on Linux the fixed CLI matched the
  patched and the unpatched oracle (`ORACLE_PATCH=`) on all 150 `make oracle`
  cases with E. coli, and every zero-gene meta contig in the synthetic corpus
  reports the first in-window bin (23/23 default, 43/43 `-c`, 24/24 `-m`,
  23/23 `-g 4`).
  The library's own meta loop skips no-path bins instead, so it can report a
  different bin for a contig where no bin has a path (to be unified in the
  shared pipeline).
- **Unbounded `add_nodes`** (`node.c:28`) overflows the node array on
  low-complexity input whose node density exceeds the capacity, e.g.
  `"ATG" × 110000`. The result is heap corruption. The fix is bounded
  growth or a clean error, not truncation that alters results.
- **Unchecked allocations:** `calc_mer_bg` does not check its `malloc`
  (`sequence.c:645`), and `calc_most_gc_frame` leaks on partial allocation
  failure (`sequence.c:576-580`).

## Fork deviations (CLI)

- **`-ffp-contract=off`.** The fork builds without fused multiply-adds, so
  ARM64 output equals x86-64 output. A default ARM64 build of upstream can
  differ in the last printed digit of some scores (`docs/testing.md`).
  `make oracle` builds the oracle with the same flag.

- **Training progress on stderr** (since `fea252a`). Without `-q`, upstream
  prints progress lines while training ("Looking for GC bias in different
  frames…frame bias scores: …", "Building initial set of genes to train
  from…" and so on). The fork's CLI trains through the library's
  `prodigal_train`, which prints nothing, so in single mode and with `-t`
  the fork prints only "Locating all potential starts and stops...done!".
  Output files are unaffected. To be restored through the log callback.
