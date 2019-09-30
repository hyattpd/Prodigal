# Prodigal

  Fast, reliable protein-coding gene prediction for prokaryotic genomes.

```bash
prodigal -i my.genome.fna -o my.genes -a my.proteins.faa
prodigal -i my.metagenome.fna -o my.genes -a my.proteins.faa -p meta
prodigal -h
```

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
