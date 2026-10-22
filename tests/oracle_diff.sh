#!/usr/bin/env bash
# Compare this fork's CLI with upstream Prodigal c1e2d36 (the oracle).
#
# Every case runs both binaries with the same arguments and compares the exit
# status, stderr, and every output file. Text is compared after normalising
# each side's version string (tests/lib.sh); training files byte for byte.
# Both binaries are built on this machine, with the same compiler and flags,
# so the comparison holds on any platform.
#
# Usage: tests/oracle_diff.sh [--stamp-key]
#   --stamp-key  print the oracle build's cache key and exit (used by CI)
# Environment:
#   PRODIGAL      binary under test (default ./prodigal in the repo root)
#   ORACLE        oracle binary; when unset, c1e2d36 + ORACLE_PATCH is built
#                 into ORACLE_DIR (default tests/out/oracle-c1e2d36) with $CC
#                 (default gcc), -ffp-contract=off and $CFLAGS, and rebuilt
#                 when its stamp (key below) changes. Needs full history (not
#                 a shallow clone).
#   ORACLE_PATCH  patch applied to the oracle (default
#                 tests/oracle-c1e2d36.patch); empty builds upstream unpatched,
#                 whose meta mode crashes on macOS
#   ECOLI         optional FASTA of E. coli K-12 MG1655 (NC_000913.3); adds the
#                 genome and its 2 kbp fragments to the corpus
#   ORACLE_OUT    scratch directory (default tests/out/oracle)
#   JOBS          parallel jobs (default: number of CPUs)
# Exit status: 0 all cases identical, 1 a case differs, 2 setup error.

set -u

# Worker mode: run one binary on one case. Arguments are the binary, the
# case directory, the input ("stdin:PATH" redirects PATH to stdin, "pipe:PATH"
# pipes it) and the CLI arguments. It runs inside the case directory, so
# relative output names land there (as does the CLI's temporary copy of piped
# input); inputs must be absolute.
if [ "${1:-}" = --job ]; then
  bin=$2 dir=$3 input=$4; shift 4
  mkdir -p "$dir" && cd "$dir" || exit 1
  case $input in
    stdin:*) "$bin" "$@" <"${input#stdin:}" ;;
    pipe:*)  cat "${input#pipe:}" | "$bin" "$@" ;;
    *)       "$bin" -i "$input" "$@" ;;
  esac >stdout 2>stderr
  echo $? >status
  exit 0
fi

. "$(dirname "$0")/lib.sh"

# The oracle build's stamp: what it was built from and how. The directory is
# ours only if it holds one. "recipe" changes whenever the build steps do.
OCOMMIT=c1e2d361479cc1b18175ea79ebd8ff10411c46cb
patch=${ORACLE_PATCH-$ROOT/tests/oracle-c1e2d36.patch}
patch=${patch:+$(abspath "$patch")}
cc=${CC:-gcc}
# The fork's Makefile builds with -ffp-contract=off; without it the oracle
# would differ on ARM64 by contraction alone.
ocflags="-ffp-contract=off${CFLAGS:+ $CFLAGS}"
key="recipe=3 $OCOMMIT patch=$([ -n "$patch" ] && cksum <"$patch" || echo none)"
key="$key cc=$cc ($("$cc" --version 2>/dev/null | head -1)) cflags=$ocflags"
[ "${1:-}" = --stamp-key ] && { echo "$key"; exit 0; }

SELF=$(abspath "$0")
P=$(prodigal_bin) || exit 2
O=$(abspath "${ORACLE_OUT:-$ROOT/tests/out/oracle}")
OD=$(abspath "${ORACLE_DIR:-$ROOT/tests/out/oracle-c1e2d36}")
ORACLE=${ORACLE:+$(abspath "$ORACLE")}
ECOLI=${ECOLI:+$(abspath "$ECOLI")}
cd "$ROOT" || exit 2
xargs_safe "$ROOT" "$P" "$O" "$OD" "$ORACLE" "$ECOLI"

if [ -z "$ORACLE" ]; then
  ORACLE=$OD/prodigal
  if [ "$(cat "$OD/.oracle-stamp" 2>/dev/null)" != "$key" ] || [ ! -x "$ORACLE" ]; then
    if [ -e "$OD" ]; then
      [ -f "$OD/.oracle-stamp" ] ||
        { echo "error: $OD exists but is not an oracle build; remove it" >&2; exit 2; }
      rm -rf "$OD"
    fi
    git cat-file -e "$OCOMMIT^{commit}" 2>/dev/null ||
      { echo "error: commit $OCOMMIT not available (shallow clone?)" >&2; exit 2; }
    echo "building oracle c1e2d36 in $OD (cc=$cc)"
    # training.c is 25 MB of constant data: -O0 changes no output and cuts
    # the build from about a minute to seconds.
    mkdir -p "$OD" && echo building >"$OD/.oracle-stamp" &&
      git archive "$OCOMMIT" | tar -x -f - -C "$OD" &&
      { [ -z "$patch" ] || patch -s -p1 -d "$OD" <"$patch"; } &&
      ( unset MAKEFLAGS MFLAGS   # "make oracle CFLAGS=..." must not reach it
        make -C "$OD" CC="$cc" CFLAGS=-O0 training.o &&
        CFLAGS="$ocflags" make -C "$OD" CC="$cc" ) >"$OD/build.log" 2>&1 ||
      { echo "error: oracle build failed; see $OD/build.log" >&2; exit 2; }
    echo "$key" >"$OD/.oracle-stamp"
  fi
fi
[ -x "$ORACLE" ] || { echo "error: oracle $ORACLE not found" >&2; exit 2; }

mkdir -p "$O" || exit 2
rm -rf "$O/cases" "$O"/jobs.* "$O/cases.list"
corpus_setup "$O/corpus" || exit 2

# One job line per (case, binary), with $quiet appended (never a trailing
# blank, see run_jobs). The case list, not the directories found afterwards,
# drives the comparison, so a case that never ran fails. A case must succeed
# on the oracle unless it matches a pattern in expect_error.
add() {  # add <case> <input> <args...>
  local name=$1; shift
  echo "$name" >>"$O/cases.list"
  echo "$P $O/cases/$name/fork $*${quiet:+ $quiet}" >>"$O/jobs.$phase"
  echo "$ORACLE $O/cases/$name/oracle $*${quiet:+ $quiet}" >>"$O/jobs.$phase"
}
expect_error=("trained.*.n")   # -n with a training file is a usage error
outs="-o out -a faa -d fna -s starts"
flagsets=("" "-c" "-m" "-n" "-g 4" "-g 25")
# No -g in meta mode: the fork refuses all but 11 there, which changes
# nothing (a deviation cli_regress checks)
metasets=("" "-c" "-m" "-n")
trainsets=("" "-c" "-m" "-n" "-g 4")

phase=1 quiet=-q
for entry in "${corpus[@]}"; do
  label=${entry%%:*} in=${entry#*:}
  for mode in single meta; do
    if [ $mode = meta ]; then sets=("${metasets[@]}"); else sets=("${flagsets[@]}"); fi
    for f in "${sets[@]}"; do
      add "$label.$mode.$(fname "$f")" "$in" -p $mode $f -f gff $outs
    done
    add "$label.$mode.gbk" "$in" -p $mode -f gbk -o out
    add "$label.$mode.sco" "$in" -p $mode -f sco -o out
  done
  for f in "${trainsets[@]}"; do
    add "train.$label.$(fname "$f")" "$in" -t train.bin $f
  done
done
for entry in "gbk:$C/formats.gbk" "embl:$C/formats.embl" "gz:$C/anthus.fa.gz"; do
  label=${entry%%:*} in=${entry#*:}
  for mode in single meta; do
    for f in "" "-m"; do
      add "$label.$mode.$(fname "$f")" "$in" -p $mode $f -f gff $outs
    done
  done
done
for mode in single meta; do
  add "stdin.$mode" "stdin:$ROOT/anthus_aco.fas" -p $mode -f gff $outs
  add "pipe.$mode" "pipe:$ROOT/anthus_aco.fas" -p $mode -f gff $outs
done
quiet=   # progress messages on stderr
add loud.single "$ROOT/anthus_aco.fas" -p single -f gff -o out
add loud.meta "$ROOT/anthus_aco.fas" -p meta -f gff -o out
add loud.train "$ROOT/anthus_aco.fas" -t train.bin
run_jobs "$O/jobs.1" "$SELF" --job

# Phase 2 predicts from the training files the oracle wrote in phase 1. If
# the oracle wrote none, its training case has already failed; the
# prediction is skipped, since -t naming a missing file would make both
# binaries write it.
phase=2 quiet=-q nskip=0
for entry in "${corpus[@]}"; do
  label=${entry%%:*} in=${entry#*:}
  for f in "${trainsets[@]}"; do
    n=$(fname "$f") t=$O/cases/train.$label.$n/oracle/train.bin
    if [ -f "$t" ]; then add "trained.$label.$n" "$in" -t "$t" $f -f gff $outs
    else echo "SKIP trained.$label.$n: the oracle wrote no training file"; nskip=$((nskip + 1)); fi
  done
done
quiet=
t=$O/cases/loud.train/oracle/train.bin
[ -f "$t" ] && add loud.trained "$ROOT/anthus_aco.fas" -t "$t" -f gff -o out
run_jobs "$O/jobs.2" "$SELF" --job

ncase=0 nfail=0
while read -r name; do
  dir=$O/cases/$name bad=""
  ncase=$((ncase + 1))
  # Agreeing on a failure is not agreement on output; a crash fails either
  # way (a fork crash also shows up as a status difference).
  st=$(cat "$dir/oracle/status" 2>/dev/null || echo 999)
  if matches "$name" "${expect_error[@]}"; then
    [ "$st" -gt 0 ] && [ "$st" -lt 128 ] || bad="$bad oracle-status($st)"
  else
    [ "$st" = 0 ] || bad="$bad oracle-status($st)"
  fi
  for k in "$dir"/fork/*; do
    [ -e "$dir/oracle/${k##*/}" ] || bad="$bad ${k##*/}(extra)"
  done
  for o in "$dir"/oracle/*; do
    f=${o##*/} k=$dir/fork/$f
    if [ ! -f "$k" ]; then bad="$bad $f(missing)"; continue; fi
    case $f in
      *.bin) cmp -s "$o" "$k" || bad="$bad $f" ;;
      *) diff_norm "$o" "$k" "$dir/$f.diff" || bad="$bad $f" ;;
    esac
  done
  if [ -n "$bad" ]; then
    nfail=$((nfail + 1))
    echo "DIFF $name:$bad"
    for f in $bad; do
      [ -s "$dir/$f.diff" ] && head -6 "$dir/$f.diff" | sed 's/^/     /'
    done
  fi
done <"$O/cases.list"
echo "$((ncase - nfail))/$ncase cases identical, $nskip skipped (oracle: $ORACLE)"
[ $nfail = 0 ] && [ $nskip = 0 ]
