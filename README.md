[![Prodigal Logo](http://i57.tinypic.com/n3rygn.png)](http://prodigal.ornl.gov/)

  Fast, reliable protein-coding gene prediction for prokaryotic genomes.

```bash
prodigal -i my.genome.fna -o my.genes -a my.proteins.faa
prodigal -i my.metagenome.fna -o my.genes -a my.proteins.faa -p meta
```

### New in 3.0.0-devel.1.0

  * Automatic detection of genetic code 11 vs. 4 (new default for *-g* option).
  * Comprehensive gap/scaffold handling options (new *-e* option, deprecated the *-m* option).
  * Addition of summary statistics for each run (new *-w* option).
  * Stop codon type added to gene data.
  * Sequin table output format added (via *-f sqn*).
  * Support for long and short options added.
  * Training file creation made explicit (via *-p train*).
  * Deprecated *-p meta* and *-p single* but will continue to support them for backwards compatibility.
  * New warnings if Prodigal suspects gene decay (average gene length is low).

### Getting Started

Prodigal consists of a single binary, which is provided for Linux, Mac OS X, and Windows with each official release.  You can also install from source (you will need Cygwin or MinGW on Windows) as follows:

```bash
$ make install
```

  For more detail, see [Installing Prodigal](https://www.github.com/hyattpd/Prodigal/wiki/installation).

  To see a complete list of options:

```bash
$ prodigal --help
```

### Features

  * **Predicts protein-coding genes**: Prodigal provides fast, accurate protein-coding gene predictions in GFF3, Genbank, or Sequin table format.
  * **Handles draft genomes and metagenomes**: Prodigal runs smoothly on finished genomes, draft genomes, and metagenomes.
  * **Runs quickly**: Prodigal analyzes the *E. coli K-12* genome in 10 seconds on a modern MacBook Pro.
  * **Runs unsupervised**: Prodigal is an unsupervised machine learning algorithm.  It does not need to be provided with any training data, and instead automatically learns the properties of the genome from the sequence itself, including genetic code, RBS motif usage, start codon usage, and coding statistics.
  * **Handles gaps, scaffolds, and partial genes**: The user can specify how Prodigal should deal with gaps and has numerous options for allowing or forbidding genes to run into or span gaps.
  * **Identifies translation initiation sites**: Prodigal predicts the correct translation initiation site for most genes, and can output information about every potential start site in the genome, including confidence score, RBS motif, and much more.
  * **Outputs detailed summary statistics for each genome**: Prodigal makes available many statistics for each genome, including contig length, gene length, GC content, GC skew, RBS motifs used, and start and stop codon usage.

### More Information

  * [Website](http://prodigal.ornl.gov/)
  * [Wiki Documentation](https://github.com/hyattpd/prodigal/wiki)
  * [Options Cheat Sheet](https://github.com/hyattpd/prodigal/wiki#cheat-sheet)
  * [Google Discussion Group](https://groups.google.com/group/prodigal-discuss)

#### Contributors

 * Author: [Doug Hyatt](https://github.com/hyattpd/)

#### License

  [GPL](LICENSE)
