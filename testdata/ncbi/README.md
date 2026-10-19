# NCBI genetic codes

`gc.prt` is NCBI's genetic code table, **version 4.6**, unmodified:

- Source: <https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt>
- Retrieved: 2026-10-19
- SHA-256: `2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04`

It is produced by NCBI, a US government agency, which places no
restrictions on the use or distribution of its data.

`test_api` checks Prodigal's codon logic against it, codon by codon: a
translation table is supported only if Prodigal implements it as this file
defines it (`docs/genetic-codes.md`). To update the file, see "Changing
this" there.
