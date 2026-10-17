#!/usr/bin/env bash
# Library API ≡ CLI on the same build.
#
# For each case, runs the CLI (GFF output) and tests/api_dump (the public
# API) with the same options and compares their output in the form of
# tests/api_gff.h. Single mode is checked both training on the input, as the
# CLI does, and (the "model" cases) predicting from a CLI-written training
# file; the "-iso" cases also run with api_dump --isolated, which neutralises
# the known defects A and C.
#
# Usage: tests/api_compare.sh
# Environment:
#   PRODIGAL     CLI under test (default ./prodigal in the repo root)
#   ECOLI        optional FASTA of E. coli K-12 MG1655; adds the genome and
#                its 2 kbp fragments (tests/make_corpus.sh)
#   COMPARE_OUT  scratch directory (default tests/out/compare)
#   JOBS         parallel jobs (default: number of CPUs)
# Exit status: 0 every case as expected, 1 a case differs (or no longer
# differs, see expect_diff), 2 setup error.

set -u
. "$(dirname "$0")/lib.sh"
D=$ROOT/tests/api_dump

# Worker mode: one case. Writes RESULT (one line) into the case directory.
# Arguments: CLI, case directory, input, mode (single, meta, or model: single
# mode from a training file the CLI writes first; "-iso" appended adds
# --isolated), then options for both.
if [ "${1:-}" = --case ]; then
  P=$2 dir=$3 in=$4 mode=$5 iso=; shift 5
  mkdir -p "$dir" && cd "$dir" || exit 1
  name=${dir##*/}
  case $mode in *-iso) mode=${mode%-iso} iso=--isolated ;; esac
  if [ "$mode" = model ]; then
    mode=single
    "$P" -i "$in" -t train.bin -q 2>cli.err ||
      { echo "FAIL $name: CLI training failed" >RESULT; exit 0; }
    set -- "$@" -t "$dir/train.bin"
  fi
  if ! "$P" -i "$in" -p "$mode" "$@" -f gff -q -o cli.gff 2>cli.err; then
    echo "FAIL $name: CLI failed: $(head -1 cli.err)" >RESULT
  elif ! "$D" $iso -p "$mode" "$@" "$in" >api.txt 2>api.err; then
    echo "FAIL $name: API failed: $(head -1 api.err)" >RESULT
  elif ! "$D" --from-cli $iso -p "$mode" cli.gff >cli.txt 2>>api.err; then
    echo "FAIL $name: cannot translate the CLI's GFF: $(tail -1 api.err)" >RESULT
  elif ! grep -q '^#' cli.txt; then
    echo "FAIL $name: no contigs in the CLI's GFF" >RESULT
  elif diff cli.txt api.txt >diff; then
    echo "OK   $name ($(grep -vc '^#' cli.txt) genes)" >RESULT
  else
    echo "DIFF $name: $(grep -c '^<' diff) of $(wc -l <cli.txt | tr -d ' ') lines differ" >RESULT
  fi
  exit 0
fi

SELF=$(abspath "$0")
P=$(prodigal_bin) || exit 2
O=$(abspath "${COMPARE_OUT:-$ROOT/tests/out/compare}")
ECOLI=${ECOLI:+$(abspath "$ECOLI")}
cd "$ROOT" || exit 2
[ -x "$D" ] || { echo "error: $D not found; run make tests/api_dump" >&2; exit 2; }
xargs_safe "$ROOT" "$P" "$O" "$ECOLI"

mkdir -p "$O" || exit 2
rm -rf "$O/cases" "$O/jobs" "$O/cases.list"
corpus_setup "$O/corpus" || exit 2

# Cases that differ today because of the library defects A-C
# (docs/testing.md); P1 fixes them and empties the list. A listed case must
# keep differing, so the list stays current; every other case must match.
# Single mode trains wrongly (A); meta reports recovery-pass scores (B,
# anthus happens to agree); prediction from a model misreports non-SD motifs
# (C, so E. coli agrees). The -iso cases are never listed.
expect_diff=("*.single.*" "synthetic.meta.*" "big.meta.*" "ecoli.meta.*" "ecoli2k.meta.*"
             "anthus.model" "synthetic.model" "big.model")

add() {  # add <case> <input> <mode> <options...>
  local name=$1; shift
  echo "$name" >>"$O/cases.list"
  echo "$P $O/cases/$name $*" >>"$O/jobs"
}
for entry in "${corpus[@]}"; do
  label=${entry%%:*} in=${entry#*:}
  for f in "" "-c" "-m"; do
    add "$label.meta.$(fname "$f")" "$in" meta $f
  done
  for mode in single single-iso; do
    for f in "" "-c" "-m" "-n" "-g 4"; do
      add "$label.$mode.$(fname "$f")" "$in" $mode $f
    done
  done
  add "$label.model" "$in" model
  add "$label.model-iso" "$in" model-iso
done
run_jobs "$O/jobs" "$SELF" --case

ncase=0 nfail=0 nknown=0
while read -r name; do
  ncase=$((ncase + 1)) show=
  r=$(cat "$O/cases/$name/RESULT" 2>/dev/null) || r="FAIL $name: did not run"
  if matches "$name" ${expect_diff[@]+"${expect_diff[@]}"}; then
    case $r in
      DIFF*) nknown=$((nknown + 1)) r="$r (expected)" ;;
      OK*) nfail=$((nfail + 1)) r="$r; now identical: remove it from expect_diff" ;;
      *) nfail=$((nfail + 1)) ;;
    esac
  else
    case $r in
      OK*) ;;
      DIFF*) nfail=$((nfail + 1)) show=1 ;;
      *) nfail=$((nfail + 1)) ;;
    esac
  fi
  echo "$r"
  [ -z "$show" ] || head -4 "$O/cases/$name/diff" | cut -c1-240 | sed 's/^/     /'
done <"$O/cases.list"
echo "$((ncase - nfail - nknown))/$ncase cases identical, $nknown expected to differ"
[ $nfail = 0 ]
