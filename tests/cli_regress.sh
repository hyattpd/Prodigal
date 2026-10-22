#!/usr/bin/env bash
# CLI regression against the committed upstream ground truth (docs/testing.md).
#
# Runs the CLI on anthus_aco.fas for every output type in both modes and
# compares each result with testdata/ground_truth/. Also predicts from the
# committed ref_train.bin, and writes a training file, predicts from it and
# byte-compares it with ref_train.bin. Last, the fork's deviations from
# upstream (docs/upstream-behavior.md) are checked against its own output.
#
# Usage: tests/cli_regress.sh [--binary | --no-binary]
#   The training-file byte compare is on by default except on ARM64 macOS,
#   where it is replaced by a determinism check (see docs/testing.md); the
#   flags force either behaviour.
# Environment:
#   PRODIGAL     binary under test (default ./prodigal in the repo root)
#   REGRESS_OUT  scratch directory for outputs (default tests/out/regress)
# Exit status: 0 all checks passed, 1 a check failed, 2 setup error.

set -u
. "$(dirname "$0")/lib.sh"

case $(uname -sm) in "Darwin arm64") binary=0 ;; *) binary=1 ;; esac
case "${1:-}" in
  --binary) binary=1 ;;
  --no-binary) binary=0 ;;
  "") ;;
  *) echo "usage: $0 [--binary | --no-binary]" >&2; exit 2 ;;
esac

P=$(prodigal_bin) || exit 2
O=$(abspath "${REGRESS_OUT:-$ROOT/tests/out/regress}")
cd "$ROOT" || exit 2
G=testdata/ground_truth
IN=anthus_aco.fas

# The CLI writes a training file when -t names a missing one; never let it
# create a reference.
[ -f $G/ref_train.bin ] || { echo "error: missing $G/ref_train.bin" >&2; exit 2; }
mkdir -p "$O" || exit 2
# Remove only what this script writes, so a stale result can never pass.
rm -f "$O"/meta.* "$O"/single.* "$O"/trained.* "$O"/fromref.* \
      "$O"/train*.bin "$O"/stderr "$O"/empty* "$O"/noempty* "$O"/g26.* "$O"/g32.*

fail=0
run() {  # run the CLI; any non-zero exit is a failure
  if ! "$P" "$@" 2>"$O/stderr"; then
    echo "FAIL prodigal $*"; sed 's/^/     /' "$O/stderr"; fail=1
  fi
}
expect_exit() {  # expect_exit CODE PATTERN OPTIONS...: the CLI on $IN fails so
  local code=$1 pattern=$2; shift 2
  "$P" -i $IN -q -o /dev/null "$@" 2>"$O/stderr"
  local rc=$?
  if [ $rc = "$code" ] && grep -q "$pattern" "$O/stderr"; then
    echo "OK   $* rejected"
  else
    echo "FAIL $* (exit $rc, want $code and \"$pattern\")"; sed 's/^/     /' "$O/stderr"; fail=1
  fi
}
check() {  # check <output> <reference>
  if diff_norm "$G/$2" "$O/$1" "$O/$1.diff"; then
    echo "OK   $1"
  else
    echo "DIFF $1 (vs $2)"; head -20 "$O/$1.diff" | sed 's/^/     /'; fail=1
  fi
}

# Outputs are named like their references: <mode>.<type> vs ref_<mode>.<type>.
for mode in meta single; do
  run -i $IN -p $mode -f gff -q -o "$O/$mode.gff" \
      -a "$O/$mode.proteins" -d "$O/$mode.nucl" -s "$O/$mode.starts"
  run -i $IN -p $mode -f gbk -q -o "$O/$mode.gbk"
  run -i $IN -p $mode -f sco -q -o "$O/$mode.sco"
  for t in gff gbk sco proteins nucl starts; do check $mode.$t ref_$mode.$t; done
done

# Predict from the committed training file, then write one and predict from it.
run -i $IN -t $G/ref_train.bin -f gff -q -o "$O/fromref.gff"
check fromref.gff ref_trained.gff
run -i $IN -t "$O/train.bin" -q
run -i $IN -t "$O/train.bin" -f gff -q -o "$O/trained.gff"
check trained.gff ref_trained.gff
if [ $binary = 1 ]; then
  want=$G/ref_train.bin
else
  run -i $IN -t "$O/train2.bin" -q
  want=$O/train2.bin
fi
if cmp "$want" "$O/train.bin"; then
  echo "OK   train.bin (vs ${want##*/})"
else
  echo "DIFF train.bin (vs ${want##*/})"; fail=1
fi

# Deviations. An empty record (header, no sequence), first or later, is
# skipped with a warning: the genes are those of the file without it, in
# FASTA and in EMBL (default names) alike. Upstream stopped reading there.
# Training on the input still adds a spacer for it, as upstream did, so
# there only the record count is checked.
fa2embl='/^>/ { if (n++) print "//"; print "ID   " substr($0, 2) ";"; print "SQ   Sequence;"; next }
         { print } END { print "//" }'
# Both sides of each comparison come from awk, so they have the same line
# endings (awk may drop the CRs of anthus_aco.fas, and headers keep them).
awk '{ print }' $IN >"$O/noempty.fa"
awk '/^>/ { n++; if (n == 1 || n == 4) print ">empty_record_" n } { print }' \
    $IN >"$O/empty.fa"
awk "$fa2embl" $IN >"$O/noempty.embl"
awk "$fa2embl" "$O/empty.fa" >"$O/empty.embl"
empty_case() {  # empty_case WITHOUT WITH OPTIONS...: same genes, 2 warnings
  local without=$1 with=$2; shift 2
  run -i "$without" "$@" -f gff -q -o "$O/noempty.gff"
  run -i "$with" "$@" -f gff -q -o "$O/empty.gff"
  if [ "$(grep -c 'no sequence' "$O/stderr")" != 2 ]; then
    echo "FAIL empty records (${with##*/} $*): expected 2 warnings"; sed 's/^/     /' "$O/stderr"; fail=1
  elif ! cmp -s "$O/noempty.gff" "$O/empty.gff"; then
    echo "DIFF empty records (${with##*/} $*): genes differ from the file without them"; fail=1
  else
    echo "OK   empty records skipped (${with##*/} $*)"
  fi
}
empty_case "$O/noempty.fa" "$O/empty.fa" -p meta
empty_case "$O/noempty.fa" "$O/empty.fa" -t $G/ref_train.bin
empty_case "$O/noempty.embl" "$O/empty.embl" -p meta
run -i "$O/empty.fa" -f gff -q -o "$O/empty.gff"
n=$(grep -c '^# Sequence Data' "$O/empty.gff")
if [ "$n" = "$(grep -c '^>' $IN)" ]; then
  echo "OK   empty records skipped (training on the input)"
else
  echo "FAIL empty records (training on the input): $n records in the output"; fail=1
fi

# Input from a closed standard input names /dev/stdin (upstream: "(null)").
"$P" -p meta -q <&- >/dev/null 2>"$O/stderr"
rc=$?
if [ $rc = 5 ] && grep -q "can't open input file /dev/stdin" "$O/stderr"; then
  echo "OK   closed stdin"
else
  echo "FAIL closed stdin (exit $rc)"; sed 's/^/     /' "$O/stderr"; fail=1
fi

# Translation tables Prodigal does not use as NCBI defines them are
# rejected (docs/genetic-codes.md), with the reason; upstream used 1, 2, 3,
# 6, 10, 12, 14, 16, 22, 23 and 24.
for t in 1 2 3 6 10 12 14 16 22 23 24 27 29 7; do
  expect_exit 15 "^Translation table $t is not supported: " -g $t
done
# Meta mode chooses table 11 or 4 for each sequence, so it refuses any -g
# but 11, before or after -p (upstream ignored -g there); -g 11 changes
# nothing
expect_exit 15 "^Translation table 1 is not supported: Prodigal" -p meta -g 1
expect_exit 15 "^Translation table 4 is not supported: metagenomic mode" -p meta -g 4
expect_exit 15 "^Translation table 25 is not supported: metagenomic mode" -g 25 -p meta
# A table upstream rejects still fails at -g, before a later -h or -g
expect_exit 15 "^Translation table 7 is not supported: NCBI" -g 7 -h
expect_exit 15 "^Translation table 7 is not supported: NCBI" -g 7 -g 11
run -i $IN -p meta -g 11 -f gff -q -o "$O/meta.g11.gff"
check meta.g11.gff ref_meta.gff
# A table that is not a number, and a missing parameter at the end
# (upstream read argv[argc] after a final -G or -d)
expect_exit 15 "^Translation table 'abc' is not a number" -g abc
expect_exit 15 "^Translation table '11x' is not a number" -g 11x
expect_exit 15 "options require parameters" -G
expect_exit 15 "options require parameters" -d
# So is a training file that uses one (upstream's -g 24 -t writes one). The
# table is the int after the model's first double (struct _training).
cp $G/ref_train.bin "$O/train24.bin"
if [ "$(od -An -tu4 -j8 -N4 "$O/train24.bin" | tr -d ' ')" != 11 ]; then
  echo "FAIL ref_train.bin: no table 11 at byte 8"; fail=1
else
  printf '\030\000\000\000' | dd of="$O/train24.bin" bs=1 seek=8 conv=notrunc 2>/dev/null
  expect_exit 4 "translation table 24, which is not supported: " -t "$O/train24.bin"
fi

# Table 32 (Balanophoraceae plastid), which upstream rejects, is supported:
# genes read through TAG, so they differ from table 11's (single.gff)
run -i $IN -g 32 -f gff -q -o "$O/g32.gff"
if grep -q "transl_table=32;" "$O/g32.gff" && grep -q "CDS" "$O/g32.gff" &&
   ! cmp -s <(grep CDS "$O/single.gff" | cut -f4,5,7) \
            <(grep CDS "$O/g32.gff" | cut -f4,5,7); then
  echo "OK   -g 32 calls genes with table 32"
else
  echo "FAIL -g 32: no genes, not with table 32, or table 11's genes"; fail=1
fi

# So is table 26 (Pachysolen nuclear). Of ATG, GTG and TTG it starts genes
# only at ATG, so no gene starts at GTG or TTG, as many of table 11's do
run -i $IN -g 26 -f gff -q -o "$O/g26.gff"
if grep -q "transl_table=26;" "$O/g26.gff" && grep -q "CDS" "$O/g26.gff" &&
   grep -q "start_type=[GT]TG" "$O/single.gff" &&
   ! grep -q "start_type=[GT]TG" "$O/g26.gff"; then
  echo "OK   -g 26 calls genes with table 26, starting at ATG only"
else
  echo "FAIL -g 26: no genes, not with table 26, or a GTG or TTG start"; fail=1
fi

# A training file written with either table carries it: predicting from the
# file calls the same genes as -g does
for t in 26 32; do
  run -i $IN -g $t -t "$O/g$t.bin" -q
  run -i $IN -t "$O/g$t.bin" -f gff -q -o "$O/g$t.t.gff"
  if grep -q "transl_table=$t;" "$O/g$t.t.gff" &&
     cmp -s <(grep CDS "$O/g$t.gff") <(grep CDS "$O/g$t.t.gff"); then
    echo "OK   -g $t -t writes a model that carries table $t"
  else
    echo "FAIL -g $t -t: predicting from the model is not -g $t"; fail=1
  fi
done

[ $fail = 0 ] && echo "all CLI regression checks passed"
exit $fail
