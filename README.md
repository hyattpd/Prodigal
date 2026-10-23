# Prodigal

  Fast, reliable protein-coding gene prediction for prokaryotic genomes.

```bash
prodigal -i my.genome.fna -o my.genes -a my.proteins.faa
prodigal -i my.metagenome.fna -o my.genes -a my.proteins.faa -p meta
prodigal -h
```

### Genetic codes: this fork rejects the tables Prodigal gets wrong

  Upstream Prodigal 2.6.3 accepts 19 translation tables (`-g`). In 11 of them, Prodigal's handling disagrees with NCBI's genetic code table, and it uses them without a warning: **1, 2, 3, 6, 10, 12, 14, 16, 22, 23 and 24**. With them it starts or ends genes at the wrong codons, scores open reading frames against the wrong stop codons, or prints wrong amino acids. **This fork refuses those tables with an error.** It accepts a table only if Prodigal's stops, amino acids, ATG/GTG/TTG starts and stop statistics all match NCBI's:

  * The accepted tables are 4, 5, 9, 11 (the default), 13, 15, 21 and 25, with results unchanged from upstream, and 26 (Pachysolen tannophilus nuclear) and 32 (Balanophoraceae plastid), which NCBI added after Prodigal 2.6.
  * Metagenomic mode (`-p meta`) chooses table 11 or 4 for each sequence from its pretrained models, so its results are unaffected. Upstream silently ignored `-g` there; this fork refuses any `-g` but 11 in metagenomic mode, and `-n`, which it would ignore too. For a genome that uses another table (most mitochondria, say), use single mode.
  * Prodigal never starts a gene at a codon other than ATG, GTG or TTG, whatever the table, including the default, table 11. NCBI lists other start codons, such as CTG, ATT and ATA, for some tables, and Prodigal never uses them. In E. coli K-12, 6 of 4,300 annotated genes start at one (infC among them); Prodigal calls 4 of them from another start and misses 2. This is a limit of Prodigal itself, which the fork documents rather than changes.

  The reasons for each table, the start codons Prodigal never uses, and how this is tested against NCBI's own file are in [docs/genetic-codes.md](docs/genetic-codes.md).

### New in 2.6.4 (April 2026)
  * Reentrant API created

### New in 2.6.3 (February 2016)
  * Fixed a bug in protein translation output of partial genes where TTG/GTG
codons were being incorrectly translated to methionine.

### Getting Started

Prodigal consists of a single binary, which is provided for Linux, Mac OS X, and Windows with each official release.  You can also install from source (you will need Cygwin or MinGW on Windows) as follows:

```bash
$ make install
```

  For more detail, see [Installing Prodigal](https://www.github.com/hyattpd/Prodigal/wiki/installation).

  To see a complete list of options:

```bash
$ prodigal -h
```

### Features

  * **Predicts protein-coding genes**: Prodigal provides fast, accurate protein-coding gene predictions in GFF3, Genbank, or Sequin table format.
  * **Handles draft genomes and metagenomes**: Prodigal runs smoothly on finished genomes, draft genomes, and metagenomes.
  * **Runs quickly**: Prodigal analyzes the *E. coli K-12* genome in 10 seconds on a modern MacBook Pro.
  * **Runs unsupervised**: Prodigal is an unsupervised machine learning algorithm.  It does not need to be provided with any training data, and instead automatically learns the properties of the genome from the sequence itself, including RBS motif usage, start codon usage, and coding statistics.
  * **Handles gaps and partial genes**: The user can specify if Prodigal should build genes across runs of N's as well as how to handle genes at the edges of contigs.
  * **Identifies translation initiation sites**: Prodigal predicts the correct translation initiation site for most genes, and can output information about every potential start site in the genome, including confidence score, RBS motif, and much more.

### More Information

  * [Website](http://prodigal.ornl.gov/)
  * [Wiki Documentation](https://github.com/hyattpd/prodigal/wiki)
  * [Options Cheat Sheet](https://github.com/hyattpd/prodigal/wiki#cheat-sheet)
  * [Google Discussion Group](https://groups.google.com/group/prodigal-discuss)

#### Contributors

 * Author: [Doug Hyatt](https://github.com/hyattpd/)

#### License

  [GPL](LICENSE)
