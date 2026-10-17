# Prodigal internals

How the algorithm and the library are put together, written for someone
changing the code. Line references are to `v2.6.4-miint` at `dd88690`.

## Sequence encoding

Each sequence lives in three bitmaps inside the context (`seq`, `rseq`,
`useq`, each sized for `MAX_SEQ`):

- **`seq`**: 2 bits per base. Base `i` uses bits `2i` and `2i+1`:
  `A=00`, `G=01` (bit `2i` set), `C=10` (bit `2i+1` set), `T=11`.
  Complementing is bit inversion (A↔T, G↔C). `mer_ndx(len, seq, pos)` packs
  `len` bases into an integer, base `k` contributing `val << 2k`.
  `mer_text` inverts that using `{A, G, C, T}`.
- **`rseq`**: the reverse complement, built by `rcom_seq`
  (`sequence.c:243`). An index `i` in `rseq` is forward position
  `slen-1-i`.
- **`useq`**: 1 bit per base, set for any character other than A/C/G/T
  (upper or lower case). Such bases are **encoded as C** in `seq`.
  `rcom_seq` toggles them back so they are C on both strands. Code that cares
  checks `is_n(useq, forward_pos)`; on the reverse strand that is
  `is_n(useq, slen-1-j)`.

Input filtering is the same in the CLI readers (`sequence.c:70`, `:181`) and
in the API (`prodigal_set_sequence`, `prodigal_set_training_sequences`):

- Characters outside `'A'..'z'` are skipped and do not advance the position.
  That includes `-`, `.`, `*`, digits and whitespace, so coordinates refer to
  the *filtered* sequence.
- The six characters between `Z` and `a` (`[ \ ] ^ _ \``) and all IUPAC
  codes count as ambiguous.

**Masking** (`-m` / `mask_regions`): a run of at least `MASK_SIZE` (50) `N`/`n`
becomes a mask region, recorded when the run is closed by a non-N. A run
reaching the end of the sequence is never recorded. `cross_mask()` stops
nodes from spanning a mask.

**Training input** (single mode) is all contigs concatenated, with
`TTAATTAATTAA` between them and after the last one when there is more than
one. That spacer puts a stop in all six frames. The training GC is computed
over the concatenation, spacers included.

## Nodes

A **node** (`struct _node`, `node.h`, 176 bytes) is a potential start or stop
codon on one strand.

`add_nodes` (`node.c:28`) scans each strand from the 3' end. In every frame
it records each stop, and each start (ATG/GTG/TTG, as allowed by
`is_start(trans_table)`) that is at least `MIN_GENE` (90 bp) from its stop,
or `MIN_EDGE_GENE` (60 bp) for genes that run off an edge. Fields:

| Field | Meaning |
|---|---|
| `type` | `ATG=0`, `GTG=1`, `TTG=2`, `STOP=3` |
| `strand` | `+1` / `-1` |
| `ndx` | 0-based forward coordinate of the codon. For reverse-strand nodes this is the codon's rightmost forward base (`slen-1-i`). |
| `stop_val` | For a start: the `ndx` of its stop. For a stop: the previous stop in frame. |
| `edge` | 1 if the node runs off the sequence end (no real start or stop). |
| `cscore`, `sscore`, `rscore`, `uscore`, `tscore` | Coding, start, RBS, upstream and type scores. |
| `rbs[2]`, `mot` | Best SD motif (exact and mismatch) and best non-SD upstream motif. |
| `gc_cont`, `gc_score[3]`, `gc_bias` | Per-ORF GC and frame GC statistics. |
| `traceb`, `tracef`, `score`, `ov_mark`, `star_ptr[3]`, `elim` | Dynamic programming state. |

Nodes are sorted with `compare_nodes` (`ndx` ascending, `+` strand first).
`add_nodes` has **no capacity bound**: the caller must size `nodes` for the
worst case. See "Limits and memory".

## Scoring

`score_nodes` (`node.c:375`) runs `calc_orf_gc`, then `raw_coding_score`
(6-mer log-odds plus a length factor), then RBS scoring. That is either
`rbs_score` (SD) or `find_best_upstream_motif` (non-SD). Finally it combines
everything per start:

- `tscore` = `type_wt[type] * st_wt`. `rscore` is the SD or motif score.
  `uscore` comes from `score_upstream_composition`.
- Edge handling: `EDGE_BONUS` / `EDGE_UPS`. Starts within 3 bp of an end are
  turned into edge nodes, which **mutates `edge`**; see
  `upstream-behavior.md`.
- Penalties: short genes (< 250 bp), and in meta mode, fragments under 3 kbp
  (`META_PEN`).
- `sscore = tscore + rscore + uscore`, adjusted for negative coding.

After the dynamic programming, `eliminate_bad_genes` (`dprog.c:306`) adds
`intergenic_mod(...)` to the `sscore` of the starts on the chosen path. That
is the operon/overlap bonus or the strand-switch penalty. It then marks genes
with `cscore + sscore < 0` for elimination. **The scores the CLI reports
include this adjustment.**

`calculate_confidence(score, st_wt)` (`gene.c:524`) is
`logistic(score / st_wt) * 100`, clamped to `[50, 99.99]`, where
`score = cscore + sscore`.

## Training (single-genome mode)

`pdg_train` (`prodigal_api.c`; `prodigal_train` and the CLI call it) runs
on the concatenated training sequence. It requires at least 20,000 bp; the
CLI warns below 100,000 bp.

1. `add_nodes`, sort.
2. `calc_most_gc_frame`, then `record_gc_bias`, which sets `tinf->bias[3]`.
3. `record_overlapping_starts(flag=0)` and `dprog(flag=0)`: an initial gene
   set chosen by GC frame bias alone.
4. `calc_dicodon_gene`, which sets `tinf->gene_dc[4096]`, then
   `raw_coding_score`.
5. `rbs_score`, then `train_starts_sd` (`type_wt`, `rbs_wt`, `ups_comp`),
   then `determine_sd_usage`.
6. If `force_nonsd`, set `uses_sd = 0`. If `!uses_sd`, run
   `train_starts_nonsd` (`mot_wt`, `no_mot`, `type_wt`, `ups_comp`). Its
   stack frame is 1.31 MB.

`tinf->gc` **must be set before step 1**, because `raw_coding_score` and the
`ups_comp` log conversion read it. The CLI sets it in `read_seq_training`
(called from `main.c`). The model (`struct _training`, 558,392 bytes) is
pointer-free, so `-t` / `prodigal_export_training` simply `fwrite`/`memcpy`
it. The binary layout depends on the architecture's padding.

## Gene finding

### Single-genome mode (per contig)

`pdg_score_nodes`, then `pdg_single_genes` (`prodigal_api.c`):

```
add_nodes → qsort → score_nodes(is_meta=0) → record_overlapping_starts(1)
→ dprog(1) → eliminate_bad_genes → add_genes → tweak_final_starts
→ record_gene_data
```

- `add_genes` (`gene.c:25`) walks the path into `struct _gene`, up to
  `MAX_GENES` (30,000) per contig. Genes beyond that are silently dropped by
  the library; the CLI warns.
- `tweak_final_starts` may move a gene's start to a nearby alternative.
- `record_gene_data` captures what is reported for each gene **from the
  node scores at that moment** into `struct _gene`: strand, partial flags,
  start type, the scores, confidence, `gc_cont`, and the RBS motif and
  spacer text. The library's output (`extract_soa`) is copied from those
  fields; the CLI prints the GFF attribute strings (`gene_data` and
  `score_data`, 500 bytes each) that are formatted from them.

### Metagenomic mode (per contig)

`pdg_meta` (`prodigal_api.c`):

1. GC window from the contig's GC:
   - `low = 0.88495*gc - 0.0102337` (at most 0.65)
   - `high = 0.86596*gc + 0.1131991` (at least 0.35)
2. For each of the 50 bins (`NUM_META`):
   - Rebuild the nodes, but only when the bin's `trans_table` differs from
     the previous bin's. Bins use table 11 (47 bins) or 4 (3 bins).
   - Skip the bin if its GC is outside `[low, high]`.
   - Otherwise run `reset_node_scores`, `score_nodes(is_meta=1)`,
     `record_overlapping_starts`, `dprog`.
   - If the path score is the best so far: `eliminate_bad_genes`,
     `add_genes`, `tweak_final_starts`, `record_gene_data`. The genes and
     their attribute strings for the best bin are captured **inside the
     loop**.
3. **Recovery pass** (the CLI only, `pdg_score_nodes` after `pdg_meta`):
   rebuild and re-score the nodes with the best bin's model, *without*
   `dprog` / `eliminate_bad_genes`. The CLI uses these nodes for the `-s`
   start file, the GFF column-6 score, and the strand and edge flags in
   `-a`/`-d`. Their `sscore` **lacks** the `intergenic_mod` adjustment. The
   library reads everything from `struct _gene` and skips this pass.

Meta mode ignores `trans_table`, `start_weight` and `force_nonsd`; every bin
carries its own model. That is true of both the CLI and the API.

### Output coordinates

`add_genes` produces 1-based, inclusive, forward-strand coordinates with
`begin < end` on both strands. The stop codon is included.

- Forward: `begin = start.ndx+1`, `end = stop.ndx+3`.
- Reverse: `begin = stop.ndx-1`, `end = start.ndx+1`.

Partial flags (`gene.c:271`):

- `partial_left` is set when a forward gene's start is an edge node, or a
  reverse gene's stop is.
- `partial_right` is the mirror case.
- `start_type` is 3 ("Edge") when the start node is an edge node.

### Translation and gene sequences (CLI `-a` / `-d`)

`write_translations` (`gene.c:450`) walks the codons from `begin` to `end`,
in reverse complement for `-` genes:

- **First codon:** `M` when it is a valid start for the table and the gene is
  not edge-started (`is_init`).
- **Middle codons:** `amino(seq, pos, tinf, 0)`, which applies the per-table
  differences.
- **Any ambiguous base in a codon:** `X`.
- **Stop codon:** included and printed as `*`. A gene that runs off the right
  edge has no `*`.
- **Table:** the model's `trans_table`. In meta mode that is the best bin's.
- Lines wrap at 60 residues.

`write_nucleotide_seqs` prints A/C/G/T, with `N` for any ambiguous base
(IUPAC codes are not preserved). Lines wrap at 70.

## One pipeline

The CLI and the library run the same gene-finding code, declared in
`prodigal_internal.h` and defined in `prodigal_api.c` (so GPL-boundary's
`build.rs` needs no new file):

- `pdg_train`: training, in stages with a report hook (`prodigal_train`
  logs from it; the CLI prints upstream's progress text);
- `pdg_score_nodes`, then `pdg_single_genes`: single genome mode, split so
  that the CLI can write its start file (`-s`) between scoring and DP;
- `pdg_meta`: metagenomic mode, every bin whose GC window holds the contig,
  keeping the best bin's genes.

`main.c` adds only what is the CLI's own: reading input, the progress text,
the start file, and printing from `ctx->genes` and `ctx->nodes`. In
metagenomic mode it calls `pdg_score_nodes` after `pdg_meta` with the best
bin's model (the recovery pass); GFF column 6 and the `-s` table come from
those nodes (`upstream-behavior.md`, quirk 2). The library's out-of-memory
messages are upstream's, so the CLI prints them as they are.

Node array invariant: every node past `ctx->nn` is zero. `build_nodes`
clears the previous list before `add_nodes` (which sets `edge` only when it
is 1, and a stale `elim` would drop a gene), and `reserve_nodes` grows the
array zeroed. Nothing else writes nodes outside `[0, ctx->nn)`.

Known divergences at `dd88690`, still open:

- **Training GC.** `prodigal_set_training_sequences` computes `ctx->gc` but
  never sets `ctx->tinf.gc`, so API training runs with `gc = 0`. On E. coli
  the API calls 4,871 genes; the CLI calls 4,319.
- **Missing from the API:**
  - protein translations (`-a`) and gene nucleotide sequences (`-d`);
  - the all-starts table (`-s`);
  - the gene ID (`<seqnum>_<n>`) and the total `score`;
  - per-contig model data (GFF `# Model Data`) beyond `prodigal_stats_t`.
- **Input limits.** Masks beyond `MAX_MASKS` and training input beyond
  `MAX_SEQ` are dropped silently by the API; the CLI exits or warns.
  `prodigal_set_sequence` checks the raw character count against `MAX_SEQ`,
  while the CLI checks the filtered base count plus `MAX_LINE`.

## Limits and memory

| Constant | Value | Where it bites |
|---|---|---|
| `MAX_SEQ` | 32,000,000 | Per-contig maximum. Also sizes `seq` and `rseq` (8 MB each) and `useq` (4 MB) at create time. |
| `STT_NOD` | 100,000 | Initial node capacity (17.6 MB). For `slen > 800,000` the array is reallocated to `slen/8` nodes (22 B/bp) and never shrinks. |
| `MAX_GENES` | 30,000 | `genes` array (33.4 MB; `struct _gene` is 1,112 bytes). Per-contig gene cap; overflow is silently truncated. |
| `MAX_MASKS` | 5,000 | Mask list (inline in the context). |
| `MASK_SIZE` | 50 | Minimum N run to mask. |
| `MAX_LINE` | 10,000 | Header length, CLI line length. |
| `MIN_GENE` / `MIN_EDGE_GENE` | 90 / 60 | Minimum ORF length. |
| `NUM_META` | 50 | Metagenomic bins: 50 × 558,392 B = 27.9 MB, lazily allocated on first meta `find_genes`. |

- **Per context:** about 71 MB at `prodigal_create`, plus 27.9 MB in meta
  mode, plus node growth for long contigs.
- **Node capacity is not enforced.** Real genomes produce about 0.04–0.06
  nodes/bp, well under the 0.125/bp budget. Low-complexity input can exceed
  it: `"ATG" × 110000` overflows the heap in `add_nodes` in both the CLI and
  the API.

## Build integration

- **Makefile:** builds two object sets. `%.o` for the CLI keeps stderr
  messages; `%.lib.o` / `%.pic.o` define `PRODIGAL_NO_MAIN` for
  `libprodigal.a` and `libprodigal.so`. Only the CLI links zlib.
- **GPL-boundary** compiles the eight library sources itself
  (`build_prodigal()`), with `-DPRODIGAL_NO_MAIN -O3 -fvisibility=hidden`,
  into a static archive linked into a Rust binary alongside FastTree,
  SortMeRNA, Bowtie2 and RocksDB.
