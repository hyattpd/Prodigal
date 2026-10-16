#!/usr/bin/env bash
# Generate the deterministic test corpus used by tests/oracle_diff.sh.
#
# Usage: tests/make_corpus.sh OUTDIR [ECOLI_FASTA]
#
# Writes into OUTDIR:
#   synthetic.fa   edge cases built from anthus_aco.fas slices and seeded
#                  pseudo-random sequence (see the case list below)
#   formats.gbk    GenBank records, including an "Expand ... gap(N)" line
#   formats.embl   EMBL records (SQ ... //)
#   anthus.fa.gz   anthus_aco.fas, gzipped (the CLI's zlib input path)
#   big.fa         one 1.2 Mbp contig alternating anthus slices and random
#                  sequence: above 800 kbp the CLI grows its node array
#   ecoli_2k.fa    ECOLI_FASTA cut into 2,000 bp fragments, remainder dropped
#                  (only when ECOLI_FASTA is given)
#
# The generator is plain awk with a Park-Miller PRNG; every value stays an
# exact integer in a double, so output is identical on every awk.

set -eu
. "$(dirname "$0")/lib.sh"
[ $# -ge 1 ] || { echo "usage: $0 OUTDIR [ECOLI_FASTA]" >&2; exit 2; }
out=$(abspath "$1")
ecoli=${2:+$(abspath "$2")}
cd "$ROOT"
mkdir -p "$out"

# Shared by the generators below: a Park-Miller PRNG (seed it in END), R(n,
# gc) random sequence, A(off, n) a slice of the anthus concatenation, which
# the main rule collects.
lib='
function rnd() { seed = (seed * 16807) % 2147483647; return seed / 2147483647 }
function R(n, gc,   s, i, r) {
  s = ""
  for (i = 0; i < n; i++) {
    r = rnd()
    if (r < gc) s = s (rnd() < 0.5 ? "G" : "C")
    else        s = s (rnd() < 0.5 ? "A" : "T")
  }
  return s
}
function A(off, n) { return substr(cat, off % (length(cat) - n) + 1, n) }
{ sub(/\r$/, "") }
/^>/ { next }
{ cat = cat $0 }
'

# Synthetic FASTA of edge cases. A() slices carry real gene signal; R()
# gives training mass and GC extremes. Order matters: meta mode carries
# state between contigs, so no-node contigs come first and after
# gene-bearing ones.
awk "$lib"'
function rec(hdr, s, width,   i) {
  print ">" hdr
  if (width == 0) { print s; return }
  for (i = 1; i <= length(s); i += width) print substr(s, i, width)
}
function rep(unit, k,   s, i) { s = ""; for (i = 0; i < k; i++) s = s unit; return s }
# every(s, n, chars, replace): after (or, if replace, instead of) every n-th
# character of s, insert the next character of chars, cycling.
function every(s, n, chars, replace,   o, i, k) {
  o = ""; k = 0
  for (i = 1; i <= length(s); i++) {
    if (!replace || i % n) o = o substr(s, i, 1)
    if (i % n == 0) o = o substr(chars, k++ % length(chars) + 1, 1)
  }
  return o
}
function mixcase(s,   o, i) {
  o = ""
  for (i = 1; i <= length(s); i += 7)
    o = o (int(i / 7) % 2 ? tolower(substr(s, i, 7)) : substr(s, i, 7))
  return o
}
END {
  seed = 20261016

  rec("first_nonode_50", A(0, 50), 60)
  rec("cat_2990", A(1000, 2990), 60)
  rec("after_gene_1", A(1, 1), 60)
  rec("after_gene_2", A(2, 2), 60)
  rec("after_gene_59", A(3, 59), 60)
  rec("cat_2999", A(4000, 2999), 80)
  rec("after_gene_30", A(5, 30), 60)
  rec("cat_3000", A(7000, 3000), 60)
  rec("cat_3001", A(10000, 3001), 60)
  rec("cat_3010 with a description", A(13000, 3010), 60)

  rec("rand_gc30_8000", R(8000, 0.30), 60)
  rec("rand_gc50_12000", R(12000, 0.50), 80)
  rec("rand_gc70_6000", R(6000, 0.70), 7)
  # One line longer than the 10,000-char line buffer of the CLI reader.
  rec("rand_gc50_unwrapped_12500", R(12500, 0.50), 0)
  # The only contigs past the clamps of the meta GC window: high < 0.35
  # needs GC < 0.274, low > 0.65 needs GC > 0.746.
  rec("rand_gc20_6000", R(6000, 0.20), 60)
  rec("rand_gc80_6000", R(6000, 0.80), 60)
  rec("after_gene_45", A(6, 45), 60)

  split("1 2 3 4 5 10 20 30 45 59 60 61 75 89 90 91 100 120 150 180 199 200", L, " ")
  for (i = 1; i in L; i++) rec("ladder_" L[i], A(i * 977, L[i]), 60)

  rec("allN_100", rep("N", 100), 60)
  rec("cat_after_allN", A(16000, 2500), 60)
  rec("allN_1000", rep("N", 1000), 60)
  rec("allN_5000", rep("N", 5000), 60)

  # N runs around the -m threshold of 50; a trailing run is never masked.
  s = rep("N", 60) A(2000, 600) rep("N", 10) A(2600, 600) rep("N", 49) \
      A(3200, 600) rep("N", 50) A(3800, 600) rep("N", 51) A(4400, 600) \
      rep("N", 200) A(5000, 600) rep("N", 70)
  rec("nruns", s, 60)

  rec("iupac", every(A(6000, 3000), 97, "RYKMSWBDHV", 1), 60)
  rec("lower", tolower(A(9000, 3000)), 60)
  rec("mixed_case", mixcase(A(12000, 3000)), 60)
  # The reader drops - . * and digits, and codes [ and _ as ambiguous.
  rec("gapped", every(every(every(A(15000, 3000), 50, "-"), 71, ".*"), 131, "0123456789[_"), 60)

  rec("atg_repeat_3000", rep("ATG", 1000), 60)
  rec("taa_repeat_900", rep("TAA", 300), 60)
  rec("allA_2000", rep("A", 2000), 60)
  rec("", A(18000, 2500), 60)
  rec("last_nonode_20", A(7, 20), 60)
}' anthus_aco.fas >"$out/synthetic.fa"

# GenBank and EMBL readers.
awk -v gbk="$out/formats.gbk" -v embl="$out/formats.embl" "$lib"'
function seqlines(s, pre, f,   i, j, line) {
  for (i = 1; i <= length(s); i += 60) {
    line = (pre ? sprintf("%9d", i) : "    ")
    for (j = i; j < i + 60 && j <= length(s); j += 10) line = line " " tolower(substr(s, j, 10))
    print line > f
  }
}
function gb(name, s, gap) {
  printf "LOCUS       %s %d bp    DNA     linear   BCT 01-JAN-2000\n", name, length(s) > gbk
  printf "DEFINITION  %s test record.\nORIGIN\n", name > gbk
  seqlines(s, 1, gbk)
  if (gap) print "     Expand Ns gap(100)" > gbk
  print "//" > gbk
}
function em(name, s) {
  printf "ID   %s; SV 1; linear; genomic DNA; STD; PRO; %d BP.\n", name, length(s) > embl
  printf "SQ   Sequence %d BP;\n", length(s) > embl
  seqlines(s, 0, embl)
  print "//" > embl
}
END {
  gb("gb_one", substr(cat, 1, 7000), 0)
  gb("gb_gap", substr(cat, 7001, 3000), 1)
  gb("gb_three", substr(cat, 10001, 11000), 0)
  em("embl_one", substr(cat, 1, 9000))
  em("embl_two", substr(cat, 9001, 12000))
}' anthus_aco.fas

gzip -9 -n -c anthus_aco.fas >"$out/anthus.fa.gz"

# One large contig, printed as it is generated.
awk "$lib"'
END {
  seed = 4242; total = 1200000; len = 0; buf = ""
  print ">big_mosaic_1200k"
  for (k = 0; len < total; k++) {
    piece = (k % 2 == 0) ? A(k * 7919, 900) : R(300, 0.5)
    if (len + length(piece) > total) piece = substr(piece, 1, total - len)
    len += length(piece); buf = buf piece
    while (length(buf) >= 80) { print substr(buf, 1, 80); buf = substr(buf, 81) }
  }
  if (buf != "") print buf
}' anthus_aco.fas >"$out/big.fa"

if [ -n "$ecoli" ]; then
  awk -v w=2000 '
    /^>/ { s = ""; next }   # fragments never span records
    { sub(/\r$/, ""); s = s $0
      for (i = 1; i + w - 1 <= length(s); i += w) printf ">ecoli_2k_%d\n%s\n", ++n, substr(s, i, w)
      s = substr(s, i) }
  ' "$ecoli" >"$out/ecoli_2k.fa"
fi
