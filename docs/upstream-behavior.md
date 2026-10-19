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
   end become edge nodes (in `score_nodes`: `edge = 1`, scores reset).
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
10. **`dprog`'s giant-ORF exception scans from the first node.** `dprog`
    scores each node's connections to the 500 nodes before it
    (`MAX_NODE_DIST`), except for a reverse-strand start, or a
    forward-strand stop, whose own stop lies further back. For those the
    loop meant to find that stop, `while(min >= 0 && nod[i].ndx !=
    nod[i].stop_val) min--`, tests `nod[i]` instead of `nod[min]`, so it
    runs to -1 and the node is scored against every node before it. One
    long reverse ORF with many starts therefore takes quadratic time: a
    330 kbp `CAT` repeat (110,000 starts in one ORF) takes 67 s in meta mode
    and 17 s in single mode, against about 2 s for the `ATG` repeat, whose
    one long ORF is forward and has one stop (2026-10-18, Ryzen 9 9950X3D).
    It needs an ORF spanning more than 500 nodes, which typical genomes
    rarely have. Kept: the wider scan can change the chosen path, so a fix
    changes output.
11. **Low-complexity training input gives a NaN.** Training on
    low-complexity input (30 kbp of one base, N included, or of a short
    repeat such as ACGT) leaves `no_mot` NaN in the model. The CLI writes
    and loads such a model, and so does the library, which therefore
    doesn't reject non-finite values when it loads one.
12. **Training input over the length limit is truncated** (CLI). When the
    training input reaches 31,990,000 bases, the CLI warns and trains on
    what it has read (`read_seq_training`). The library returns an error
    instead. (Upstream's other caps, 30,000 genes per contig and 5,000
    masks, are gone: see "Fork deviations".)
13. **Three translation tables differ from NCBI's** (`amino`, `is_stop`
    in `sequence.c`). Table 10: TGA translates as W (NCBI: C). Table 16:
    TAG, a sense codon there, translates as X (NCBI: L). Table 24: TGA is
    a stop (NCBI: W), which also changes gene calls, and AGA and AGG
    translate as R (NCBI: S and K). The other tables match NCBI. Separately,
    an ambiguous base reads as C (quirk 3), so in table 22, where TCA is a
    stop, a TNA codon ends a gene, and its translation is X, not `*`.
    These affect the CLI's `-a` output and the library's proteins alike.

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
  The library's meta loop applies the same rule (until 2.6.5 it skipped
  no-path bins, so a contig on which no bin had a path reported bin 0).
- **The last training spacer could overflow (fixed in 2.6.5).** When the
  CLI cut multi-record training input at the limit (quirk 12) on a line
  long enough to end within 12 bases of `MAX_SEQ` (unwrapped sequence), it
  then added the trailing 12-base spacer past the end of its 2-bit buffers,
  up to 3 bytes, and `rcom_seq` read and wrote as far (ASan:
  heap-buffer-overflow in `read_seq_training`, checked on 2026-10-18 with a
  13,185-base record and a 32 Mbp unwrapped one). The buffers now have 16
  bases of slack (`PDG_SEQ_CAP`), so the model is what upstream computed
  when the bytes past its buffers happened to be free.
- **Unbounded `add_nodes` (fixed in 2.6.5).** `add_nodes` (`node.c`) wrote
  past the node array on low-complexity input with more nodes than the
  `slen/8` (at least `STT_NOD`) it is sized for, e.g. `"ATG" × 110000`: heap
  corruption. The array now grows (`internals.md`, "Limits and memory"), so
  output is unchanged wherever upstream didn't overflow.
- **Unchecked allocations (fixed in 2.6.5).** `calc_mer_bg` did not check
  its `malloc`; it now counts in its output array and allocates nothing.
  `calc_most_gc_frame` leaked on a partial allocation failure; it now makes
  one allocation, through the context's allocator.

## Fork deviations (CLI)

Each changes output only on input upstream mishandled (tests:
`docs/testing.md`).

- **Empty records are skipped.** `next_seq_multi` returned -1, its
  end-of-input value, for a record with no sequence, so a header followed by
  another header silently dropped every record after it (and a file
  starting with one failed with "no input sequences to analyze"). Now it
  warns and reads on; the genes are those of the file without the record.
  An empty last record is ignored without a warning, as before (so with
  several empty records at the end, all but the last are warned about).
  Training on the input still adds upstream's 12-base spacer for it, so
  there the model can differ slightly. The library takes one record at a
  time and rejects an empty one.
- **Closed standard input** is reported as `/dev/stdin`; upstream printed
  `input_file`, which is NULL there ("(null)" on glibc).
- **No gene cap.** Upstream kept at most 30,000 genes per contig (CLI:
  "max # of genes exceeded, truncating..."); the gene array now grows.
- **No mask cap.** With `-m`, upstream exited at the 5,001st masked region
  (exit 52 in training, 55 otherwise); the mask list now grows.
- **`-ffp-contract=off`.** The fork builds without fused multiply-adds, so
  ARM64 output equals x86-64 output. A default ARM64 build of upstream can
  differ in the last printed digit of some scores (`docs/testing.md`).
  `make oracle` builds the oracle with the same flag.
