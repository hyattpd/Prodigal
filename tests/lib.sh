# Shared helpers for tests/*.sh, sourced as
#   . "$(dirname "$0")/lib.sh"
# Callers resolve user-supplied paths with abspath, then cd "$ROOT".

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd) || exit 2

abspath() { case $1 in /*) echo "$1" ;; *) echo "$PWD/$1" ;; esac; }

# prodigal_bin: absolute path of $PRODIGAL (default: this tree's build).
# Use as P=$(prodigal_bin) || exit 2.
prodigal_bin() {
  local p
  p=$(abspath "${PRODIGAL:-$ROOT/prodigal}")
  [ -x "$p" ] || { echo "error: $p not found; run make first" >&2; exit 2; }
  echo "$p"
}

# Parallel job runners (oracle_diff.sh, api_compare.sh): a job file holds one
# job per line, run as "SCRIPT FLAG <line words>" by xargs -L 1, which splits
# lines on blanks and interprets quotes and backslashes, and joins a line
# that ends in a blank to the next.
JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}
xargs_safe() {  # xargs_safe PATH...: exit 2 unless every path is safe in a job line
  local p
  for p in "$@"; do
    case $p in *[[:space:]\'\"\\]*)
      echo "error: path contains whitespace, a quote or a backslash: $p" >&2; exit 2 ;;
    esac
  done
}
run_jobs() {  # run_jobs JOBFILE SCRIPT FLAG
  [ -s "$1" ] || return 0   # GNU xargs runs once on empty input
  xargs -P "$JOBS" -L 1 "$2" "$3" <"$1" ||
    { echo "error: job runner failed" >&2; exit 2; }
}

# corpus_setup DIR: generate the test corpus (tests/make_corpus.sh) into DIR,
# set C to DIR and corpus to its "label:path" entries. An absolute $ECOLI
# adds the genome and its 2 kbp fragments. fork_corpus holds the inputs only
# the fork runs (upstream overflows or exits on them).
corpus_setup() {
  "$ROOT/tests/make_corpus.sh" "$1" ${ECOLI:+"$ECOLI"} || return 1
  C=$1
  corpus=("anthus:$ROOT/anthus_aco.fas" "synthetic:$C/synthetic.fa" "big:$C/big.fa")
  fork_corpus=("lowcomplex:$C/lowcomplex.fa" "manymasks:$C/manymasks.fa")
  [ -z "${ECOLI:-}" ] || corpus+=("ecoli:$ECOLI" "ecoli2k:$C/ecoli_2k.fa")
}

# fname FLAGS: the case-name word for a flag set: "" -> default, "-g 4" -> g4.
fname() { local f=${1:-default}; f=${f// /}; echo "${f#-}"; }

# matches NAME PATTERN...: does NAME match one of the glob patterns?
matches() {
  local n=$1 p; shift
  for p in "$@"; do case $n in $p) return 0 ;; esac; done
  return 1
}

# The ground truth and the oracle print upstream's version; the fork prints
# PRODIGAL_VERSION_STRING. Each side is normalised for its own version only,
# so pointing a script at the wrong binary shows up as a difference.
UPSTREAM_VERSION=2.6.3
FORK_VERSION=$(sed -n 's/^#define PRODIGAL_VERSION_STRING "\(.*\)"$/\1/p' "$ROOT/prodigal.h")
[ -n "$FORK_VERSION" ] || { echo "error: PRODIGAL_VERSION_STRING not found" >&2; exit 2; }

# On Windows, text-mode fopen writes CRLF. Strip only a line-final CR: the
# ground truth also carries mid-line CRs (upstream keeps the CR of a CRLF
# FASTA header in seqhdr="..."), and those must still compare.
case $(uname -s) in MINGW* | MSYS* | CYGWIN*) EOL_CR=$'\r' ;; *) EOL_CR= ;; esac

# norm FILE VERSION: print FILE with this version's strings replaced:
# "Prodigal.vVERSION" -> "Prodigal.vX", and the stderr banner
# "PRODIGAL vVERSION [date]" -> "PRODIGAL vX [date]". The dots in VERSION
# stay regex wildcards; the digits still have to match.
norm() {
  sed -e "s/$EOL_CR\$//" -e "s/\(Prodigal[._]v\)$2/\1X/g" \
      -e "s/^PRODIGAL v$2 \[[^]]*\]/PRODIGAL vX [date]/" "$1"
}

# diff_norm UPSTREAM_FILE FORK_FILE DIFF: diff the two, each normalised for
# its own version, into DIFF.
diff_norm() {
  diff <(norm "$1" "$UPSTREAM_VERSION") <(norm "$2" "$FORK_VERSION") >"$3"
}
