#!/usr/bin/env bash
# CLI regression against the committed upstream ground truth (docs/testing.md).
#
# Runs the CLI on anthus_aco.fas for every output type in both modes and
# compares each result with testdata/ground_truth/. Also predicts from the
# committed ref_train.bin, and writes a training file, predicts from it and
# byte-compares it with ref_train.bin.
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
      "$O"/train*.bin "$O"/stderr

fail=0
run() {  # run the CLI; any non-zero exit is a failure
  if ! "$P" "$@" 2>"$O/stderr"; then
    echo "FAIL prodigal $*"; sed 's/^/     /' "$O/stderr"; fail=1
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

[ $fail = 0 ] && echo "all CLI regression checks passed"
exit $fail
