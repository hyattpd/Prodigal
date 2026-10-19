#!/usr/bin/env bash
# Library API ≡ CLI on the same build.
#
# For each case, runs the CLI (GFF and -a output) and tests/api_dump (the
# public API, with proteins) with the same options and compares their output
# in the form of tests/api_gff.h: every gene's fields and its protein. Single
# mode is checked both training on the input, as the CLI does (the model the
# API trains must also equal the CLI's -t file), and (the "model" cases)
# predicting from a CLI-written training file.
#
# Usage: tests/api_compare.sh
# Environment:
#   PRODIGAL     CLI under test (default ./prodigal in the repo root)
#   ECOLI        optional FASTA of E. coli K-12 MG1655; adds the genome and
#                its 2 kbp fragments (tests/make_corpus.sh)
#   COMPARE_OUT  scratch directory (default tests/out/compare)
#   JOBS         parallel jobs (default: number of CPUs)
# Exit status: 0 all identical, 1 a case differs, 2 setup error.

set -u
. "$(dirname "$0")/lib.sh"
D=$ROOT/tests/api_dump

# Worker mode: one case. Writes RESULT (one line) into the case directory.
# Arguments: CLI, case directory, input, mode (single, meta, or model: single
# mode from a training file the CLI writes first), then options for both.
if [ "${1:-}" = --case ]; then
  P=$2 dir=$3 in=$4 mode=$5; shift 5
  mkdir -p "$dir" && cd "$dir" || exit 1
  name=${dir##*/} api_opts=
  if [ "$mode" != meta ]; then
    "$P" -i "$in" "$@" -t train.bin -q 2>cli.err ||
      { echo "FAIL $name: CLI training failed: $(head -1 cli.err)" >RESULT; exit 0; }
    if [ "$mode" = model ]; then
      mode=single
      set -- "$@" -t "$dir/train.bin"
    else
      api_opts="--write-model api.bin"
    fi
  fi
  if ! "$P" -i "$in" -p "$mode" "$@" -f gff -q -o cli.gff -a cli.faa 2>cli.err; then
    echo "FAIL $name: CLI failed: $(head -1 cli.err)" >RESULT
  elif ! "$D" -p "$mode" "$@" $api_opts --proteins "$in" >api.txt 2>api.err; then
    echo "FAIL $name: API failed: $(head -1 api.err)" >RESULT
  elif ! "$D" --from-cli -p "$mode" --faa cli.faa cli.gff >cli.txt 2>>api.err; then
    echo "FAIL $name: cannot translate the CLI's GFF: $(tail -1 api.err)" >RESULT
  elif ! grep -q '^#' cli.txt; then
    echo "FAIL $name: no contigs in the CLI's GFF" >RESULT
  elif ! diff cli.txt api.txt >diff; then
    echo "DIFF $name: $(grep -c '^<' diff) of $(wc -l <cli.txt | tr -d ' ') lines differ" >RESULT
  elif [ -n "$api_opts" ] && ! cmp train.bin api.bin >diff; then
    echo "DIFF $name: the API's model is not the CLI's -t file" >RESULT
  else
    echo "OK   $name ($(grep -vc '^#' cli.txt) genes)" >RESULT
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

add() {  # add <case> <input> <mode> <options...>
  local name=$1; shift
  echo "$name" >>"$O/cases.list"
  echo "$P $O/cases/$name $*" >>"$O/jobs"
}
for entry in "${corpus[@]}" "${fork_corpus[@]}"; do
  label=${entry%%:*} in=${entry#*:}
  for f in "" "-c" "-m"; do
    add "$label.meta.$(fname "$f")" "$in" meta $f
  done
  for f in "" "-c" "-m" "-n" "-g 4" "-g 32"; do
    add "$label.single.$(fname "$f")" "$in" single $f
  done
  add "$label.model" "$in" model
done
run_jobs "$O/jobs" "$SELF" --case

ncase=0 nfail=0
while read -r name; do
  ncase=$((ncase + 1))
  r=$(cat "$O/cases/$name/RESULT" 2>/dev/null) || r="FAIL $name: did not run"
  echo "$r"
  case $r in
    OK*) ;;
    DIFF*) nfail=$((nfail + 1)); head -4 "$O/cases/$name/diff" | cut -c1-240 | sed 's/^/     /' ;;
    *) nfail=$((nfail + 1)) ;;
  esac
done <"$O/cases.list"
echo "$((ncase - nfail))/$ncase cases identical"
[ $nfail = 0 ]
