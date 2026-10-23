# Genetic codes (translation tables)

**This fork accepts a translation table only if Prodigal uses it exactly
as NCBI defines it.**

Upstream Prodigal 2.6.3 accepts 19 tables. In 11 of them, Prodigal's
handling disagrees with NCBI's definition, and upstream uses them without a
warning: **1, 2, 3, 6, 10, 12, 14, 16, 22, 23 and 24**. With these tables,
upstream:
- starts genes at the wrong codons (1, 3, 12, 24);
- ends genes at a codon that is not a stop (24), or at an ambiguous base
  (22);
- scores open reading frames against the wrong stop codons (1, 2, 6, 12,
  14, 22, 23);
- in training, lets genes run across the joins between contigs (6, 14);
- prints wrong amino acids in the proteins (10, 16, 24).

This fork rejects those tables with an error that says why:
- in the CLI, for `-g` and for a `-t` training file that uses one;
- in the library.

Every table it still accepts from upstream gives upstream's output, byte
for byte. It also supports tables 26 and 32, which NCBI added after
Prodigal 2.6 and which Prodigal can use faithfully with a few clauses of
the fork's.

Decided on 2026-10-19:
- Support a genetic code if and only if Prodigal uses it faithfully, and
  reject the others with an error, rather than produce wrong results
  without a warning.
- Do not change Prodigal's algorithm for now.
- So support table 32, which needs only two codon clauses and nothing else
  changed, and table 26, which needs three clauses. The clauses apply only
  to the new tables, so no table upstream accepted changes.

Decided on 2026-10-22:
- Metagenomic mode refuses any table but 11, rather than ignore it as
  upstream did (below).
- Prodigal's three start codons are a documented limit of its model, in
  table 11 too ("What every table lacks", below).

## Supported tables

| `-g` | NCBI name | Stops | Starts Prodigal uses |
|---|---|---|---|
| 4 | Mold, Protozoan and Coelenterate Mitochondrial; Mycoplasma; Spiroplasma | TAA, TAG | ATG, GTG, TTG |
| 5 | Invertebrate Mitochondrial | TAA, TAG | ATG, GTG, TTG |
| 9 | Echinoderm and Flatworm Mitochondrial | TAA, TAG | ATG, GTG |
| 11 | Bacterial, Archaeal and Plant Plastid (**the default**) | TAA, TAG, TGA | ATG, GTG, TTG |
| 13 | Ascidian Mitochondrial | TAA, TAG | ATG, GTG, TTG |
| 15 | Blepharisma Macronuclear | TAA, TGA | ATG |
| 21 | Trematode Mitochondrial | TAA, TAG | ATG, GTG |
| 25 | Candidate Division SR1 and Gracilibacteria | TAA, TAG | ATG, GTG, TTG |
| 26 | Pachysolen tannophilus Nuclear (new in the fork; upstream rejects it) | TAA, TAG, TGA | ATG |
| 32 | Balanophoraceae Plastid (new in the fork; upstream rejects it) | TAA, TGA | ATG, GTG, TTG |

Tables 15 and 26 are nuclear codes of eukaryotes, a ciliate and a yeast.
"Supported" means Prodigal uses their genetic code as NCBI defines it.
Prodigal still finds genes as in a prokaryote: it does not model introns,
so a gene with one is called wrongly or missed.

> **Metagenomic mode takes table 11 only.** Metagenomic mode (`-p meta`,
> `meta_mode`) has no model of the genome. It scores each sequence with 50
> pretrained bins and keeps the best one, and each bin brings the table it
> was trained under: 11 (47 bins) or 4 (3 *Mycoplasma* bins). So it
> chooses table 11 or 4 for each sequence, and reports the one it used
> (`transl_table=` in the GFF header; `translation_table` in the library's
> `prodigal_stats_t`).
>
> It cannot use another table faithfully. The bins are models of genomes
> that use tables 11 and 4: their coding statistics, start codons and
> stops are those genomes'. The clearest case: all 50 give every hexamer
> that starts with an in-frame TAG Prodigal's lowest coding score, and the
> 47 table-11 bins do the same for TGA. Run under a table where TAG or TGA
> is an amino acid, a bin would read through the codon but score it as
> non-coding.
>
> Upstream ignored `-g` in metagenomic mode without a warning (quirk 9 in
> `docs/upstream-behavior.md`). The fork refuses any table but 11 there,
> and a model loaded into a metagenomic context ("Where it is enforced",
> below). `-g 11` is accepted and changes nothing.
>
> For a genome that uses another table, use single mode, which needs at
> least 20,000 bases of training input. Many mitochondrial genomes are
> smaller than that, and so are the plastid genomes of *Balanophora*
> (about 16 kb), the genus table 32 is for. For those, train on several
> genomes together, or load a model trained on such a set with `-t`.

## The rules

The reference is NCBI's genetic code table `gc.prt`, version 4.6, committed
verbatim as `testdata/ncbi/gc.prt`. A table is supported only if it meets
all seven rules:

1. **NCBI defines it.**
2. **No context-dependent codons.** No codon may be a stop in some contexts
   and an amino acid in others. Prodigal ends every gene at its first
   in-frame stop codon, so it can't model such a codon.
3. **Stops and amino acids.** For all 64 codons, Prodigal agrees with NCBI
   on which are stops (`is_stop`, used for gene calls) and on the amino
   acid (`amino`, used for `-a` and the library's proteins through
   `pdg_codon_letters`, a table of its answers).
4. **Starts.** Prodigal starts genes only at ATG, GTG and TTG, its three
   start types. Of those three, it must start genes at exactly the ones
   NCBI lists as start codons for the table (`is_start`), and translate a
   gene's start codon as M.
5. **The ORF-length score.** `raw_coding_score` rewards a long open reading
   frame by how unlikely it would be by chance. For that it uses the
   probability that a codon is a stop, computed from the GC content
   (`pdg_no_stop`). That probability must count exactly the table's stops.
   Upstream counts TAA, TAG and TGA in table 11, and TAA and one of TAG and
   TGA in every other table. For other stop sets the estimate is badly
   off. At GC 0.3, for example, it assumes 6.1% of codons are stops, while
   the real share is:

   | Tables | Real share of stops |
   |---|---|
   | 6, 14 | 1.8% |
   | 22 | 8.0% |
   | 2 | 8.75% |
   | 23 | 12.25% |
6. **The training spacer.** In single-genome mode, Prodigal trains on all
   contigs joined together, with `TTAATTAATTAA` between them, so that every
   frame has a stop and no gene runs from one contig into the next. The
   spacer must hold a stop in all three frames on both strands. Its only
   stops are TAA (and TTA, in table 23). In a table where TAA is not a
   stop, genes can run across the joins.
7. **Ambiguous bases.** Prodigal reads an ambiguous base (N or another
   IUPAC code) as a C (`docs/internals.md`, "Sequence encoding"). No codon
   holding one may be a stop or a start. This holds when no stop or
   modeled start codon of the table has a C in it. Such a codon translates
   as X (`docs/upstream-behavior.md`, quirk 13). In table 22, where TCA is
   a stop, a TNA codon would end a gene.

### What every table lacks: other start codons

NCBI lists other start codons for some tables: CTG, ATT, ATC, ATA and TTA.
Prodigal has no start type for them, so it never starts a gene at one. A
gene that starts at one of them is called from another start codon, or
missed.

This is a limit of Prodigal's model, not a table error. The fork accepts
tables despite it; rejecting them would reject the default, table 11. The
start codons NCBI lists that Prodigal never uses, per supported table:

| `-g` | NCBI start codons Prodigal never uses |
|---|---|
| 4 | TTA, CTG, ATT, ATC, ATA |
| 5 | ATT, ATC, ATA |
| 11 | CTG, ATT, ATC, ATA |
| 13 | ATA |
| 26 | CTG |
| 32 | CTG, ATT, ATC, ATA |

Tables 9, 15, 21 and 25 have none.

**Table 11, the default, is one of them.** NCBI lists seven start codons
for it (ATG, GTG, TTG, CTG, ATT, ATC, ATA), and Prodigal uses the first
three. In E. coli K-12 MG1655 (NCBI's annotation of `NC_000913.3`, checked
on 2026-10-19), 6 of the 4,300 annotated genes start at one of the others:
4 at ATT, among them infC (IF-3) and pcnB, and 2 at CTG. Prodigal calls 4
of those 6 from another start codon (3 shorter, 1 longer) and misses 2.
For scale, it agrees with the annotation's start on 3,889 of the 4,300
genes and picks a different start on 239. This is a limit of Prodigal, as
published and as upstream ships it; the fork documents it rather than
change Prodigal's start model, which would mean new start types, a new
training file format and retraining the 50 metagenomic bins.

## Rejected tables

"Effect upstream" says what upstream Prodigal gets wrong with the table.
"Gene calls" means genes start, end or score wrongly; "proteins" means only
the `-a` file and the library's proteins are wrong. "None" means upstream
rejects the table too.

| `-g` | NCBI name | Why it is rejected | Effect upstream |
|---|---|---|---|
| 1 | Standard | Prodigal starts genes only at ATG, but NCBI also lists TTG (rule 4). The ORF-length score counts two of its three stops (rule 5). | gene calls |
| 2 | Vertebrate Mitochondrial | The ORF-length score leaves out the AGA and AGG stops (rule 5). | gene calls |
| 3 | Yeast Mitochondrial | Prodigal starts genes only at ATG, but NCBI also lists GTG, added in `gc.prt` 4.4 (rule 4). | gene calls |
| 6 | Ciliate, Dasycladacean and Hexamita Nuclear | The only stop is TGA. The ORF-length score and the training spacer take TAA, which is Q here, for a stop (rules 5, 6). | gene calls |
| 10 | Euplotid Nuclear | Prodigal translates TGA as W; NCBI has C (rule 3). | proteins |
| 12 | Alternative Yeast Nuclear | Prodigal starts genes at TTG, which NCBI does not list (rule 4). The ORF-length score counts two of its three stops (rule 5). | gene calls |
| 14 | Alternative Flatworm Mitochondrial | The only stop is TAG. The ORF-length score and the training spacer take TAA, which is Y here, for a stop (rules 5, 6). | gene calls |
| 16 | Chlorophycean Mitochondrial | Prodigal translates TAG as X; NCBI has L (rule 3). | proteins |
| 22 | Scenedesmus obliquus Mitochondrial | The ORF-length score leaves out the TCA stop (rule 5), and an ambiguous base in TNA makes a stop (rule 7). | gene calls |
| 23 | Thraustochytrium Mitochondrial | The ORF-length score counts two of its four stops: TAA, TAG, TGA, TTA (rule 5). | gene calls |
| 24 | Rhabdopleuridae Mitochondrial | Prodigal reads TGA as a stop (NCBI: W) and translates AGA and AGG as R (NCBI: S and K) (rule 3). It never starts genes at TTG, which NCBI lists (rule 4). | gene calls (genes end at TGA) and proteins |
| 27 | Karyorelict Nuclear | TGA is a stop or W depending on context (rule 2). | none |
| 28 | Condylostoma Nuclear | TAA and TAG are stops or Q, and TGA a stop or W, depending on context (rule 2). | none |
| 29 | Mesodinium Nuclear | Not implemented: NCBI added it after Prodigal 2.6. Even with the codons implemented, the ORF-length score and the training spacer would take TAA, which is Y here, for a stop (rules 5, 6). | none |
| 30 | Peritrich Nuclear | Not implemented, as 29; TAA is E here. | none |
| 31 | Blastocrithidia Nuclear | TAA and TAG are stops or E, depending on context (rule 2). | none |
| 33 | Cephalodiscidae Mitochondrial | Not implemented, as 29; TAA is Y here. | none |

Numbers NCBI does not define are rejected as well, as upstream rejected
them: 0 and below, 7, 8, 17 to 20, and 34 and above.

**Supporting a rejected table again** means changing Prodigal: its codon
logic, its ORF-length score (as the 3.0 rewrite did, below) or its
training spacer. A table NCBI added after Prodigal 2.6 may need less:
table 32 took two codon clauses, and table 26 three clauses (Provenance,
below).

## Where it is enforced

`pdg_trans_table_unsupported` (`sequence.c`, next to the codon logic) is the
one list. It returns NULL for a supported table, and otherwise the reason.
`prodigal_config_validate` and `prodigal_set_translation_table` add
metagenomic mode's rule, any table but 11 refused, and the CLI asks
`prodigal_config_validate` once it has read all its options. Everything
that sets a table goes through them:

| Entry point | Behaviour on a rejected table |
|---|---|
| CLI `-g` | prints the reason and the usage; exit 15, the code upstream used for an invalid table. In metagenomic mode any table but 11 is refused, whether `-g` comes before or after `-p`. `-g` is also checked with `-t`, where upstream warned and used the file's table. |
| CLI `-t`, reading a training file | prints the reason; exit 4. The file is checked as `prodigal_training_info` checks a model, which also covers GC, start weight and SD flag. Training files that upstream wrote with a rejected table are refused. |
| `prodigal_config_validate` / `prodigal_create` | `PRODIGAL_ERR_INVALID_CONFIG` / NULL; with `meta_mode`, also any table but 11 |
| `prodigal_set_translation_table` | `PRODIGAL_ERR_INVALID_INPUT`, the table unchanged; in a metagenomic context, also any table but 11 |
| `prodigal_training_info` / `prodigal_load_training` | `PRODIGAL_ERR_INVALID_INPUT`; the context is unchanged. A metagenomic context refuses any model, as the CLI refuses `-p meta -t` (exit 2). |

The messages give the table, the reason and this page, and fit in 256 bytes:

```
Translation table 10 is not supported: Prodigal translates TGA as W (NCBI: C) (the-miint/Prodigal docs/genetic-codes.md).
Translation table 4 is not supported: metagenomic mode takes table 11 only, as its bins choose table 11 or 4 for each sequence; use single mode (the-miint/Prodigal docs/genetic-codes.md).
Training data has translation table 24, which is not supported: Prodigal reads TGA as a stop (NCBI: W) ...
```

## How it is tested

- **The rules, against NCBI's file.** `test_genetic_codes_ncbi`
  (`test_api.c`) checks every table number from -1 to 99 against
  `testdata/ncbi/gc.prt`, using Prodigal's own code, so a table is
  supported if and only if it meets the seven rules. It also checks that
  every metagenomic bin's table is supported.
- **Pinned lists.** The test also holds the lists of supported and rejected
  tables. A change to Prodigal's code, to the list or to NCBI's file fails
  the test until those lists are updated.
- **The entry points.** `test_genetic_code_entry_points` checks that every
  entry point rejects the same tables with the same message.
  `test_genetic_code_meta_mode` checks that in metagenomic mode the config
  and the setter take table 11 only, with the reason, and
  `test_meta_mode_refuses_model_settings` that a metagenomic context
  refuses a model.
- **The CLI.** `tests/cli_regress.sh` checks the CLI's rejections (in
  metagenomic mode too, and that `-p meta -g 11` changes nothing), that
  `-g 26` and `-g 32` call genes with those tables, and that a training
  file written with either carries its table. `make compare` checks that
  the library and the CLI give the same genes and proteins with `-g 26`
  and `-g 32`.
- **Not checked by any test:** the copies of the lists on this page, in
  `README.md`, in `prodigal.h` and in `CHANGES`. Update them together.

Details are in `docs/testing.md`, "Genetic codes".

## Changing this

- **A new `gc.prt` from NCBI:**
  1. Replace the file and update `testdata/ncbi/README.md`.
  2. Set `GC_PRT_VERSION` in `test_api.c`.
  3. Run `make test`.
  4. Update this page and the other copies of the lists.
- **Making Prodigal faithful for a rejected table** changes upstream's
  output for it, so it is an explicit decision: record it here, in
  `CHANGES` and in `docs/upstream-behavior.md`. The test then requires
  moving the table to the supported list.

## Provenance

- **This fork's code is upstream's, unchanged, except for tables 26 and
  32.** That covers `is_stop`, `is_start` and `amino`, the ORF-length stop
  probability and the training spacer. All of it is upstream's 2.6.x code,
  identical in `v2.6.3` and `c1e2d36`. The stop probability moved,
  verbatim, into `pdg_no_stop`, which then gained table 26's clause.
- **Earlier releases accepted fewer tables.** `v2.6.0` and `v2.6.1`
  accepted tables up to 23 only. `v2.6.2` added 24 and 25.
- **The unreleased 3.0 rewrite fixed some of this.** Upstream's
  unreleased Prodigal 3.0 rewrite is on upstream `master` (`82b7aff`, 2014)
  and `develop` (`9c6483d`, 2016). Its commit `ff660d3` (2013-11-22)
  corrected the amino acids of tables 10, 16 and 24 and the start codons of
  12 and 24. Its `prob_stop` computes the stop probability from the table's
  own stops. The rewrite was never released, and none of this was ported to
  2.6.x. Tables 1 and 3 still differ from NCBI there.
- **NCBI's later changes:** NCBI added table 3's GTG start in `gc.prt` 4.4,
  and tables 26 to 33 in versions 4.1 to 4.5.
- **Table 32's codon logic is the fork's:** TAG is not a stop and
  translates as W (two clauses in `is_stop` and `amino`, marked "the
  fork's"). It needs nothing else, since its stops, TAA and TGA, are a
  pair that `pdg_no_stop` counts exactly.
- **Table 26's logic is the fork's:** genes start only at ATG (`is_start`),
  CTG translates as A (`amino`), and the ORF-length score counts TAA, TAG
  and TGA, as for table 11 (`pdg_no_stop`). Each of the three clauses is
  marked "the fork's". Its stops are table 11's, so `is_stop` needs
  nothing.
- **Training files for tables 26 and 32 work only with this fork.**
  Upstream Prodigal does not check the table of a `-t` file, so it would
  load one and use the table wrongly. With table 26 it would start genes at
  GTG and TTG, translate CTG as L and count two of the three stops in the
  ORF-length score; with table 32 it would end genes at TAG.
- **Last check:** verified on 2026-10-19 against `gc.prt` version 4.6.
