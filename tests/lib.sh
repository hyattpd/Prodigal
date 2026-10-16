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
