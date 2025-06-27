#!/bin/bash
# compare.sh: compare two versions of prodigal
# ./compare.sh <old_binary> <input_file>

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <old_binary> <input_file> [additional opts]"
    echo "Example: $0 ./prodigal.windows.exe anthus_aco.fas -t test/anthus_aco.train"
    exit 1
fi

OLD_BINARY="$1"
INPUT_FILE="$2"
OUT_PREFIX="test/$(basename "$INPUT_FILE")"
shift 2

if ! command -v "$OLD_BINARY" >/dev/null; then
    echo "Error: Old binary '$OLD_BINARY' is not executable. Did you forget a './' prefix?"
    exit 1
fi
if [[ ! -f $INPUT_FILE ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

mkdir -p "test"
rm -f "$OUT_PREFIX".old.* "$OUT_PREFIX".new.* "$OUT_PREFIX".*.diff

time "$OLD_BINARY" -i "$INPUT_FILE" -f gbk -o "$OUT_PREFIX".old.gbk -d "$OUT_PREFIX".old.ffn -a "$OUT_PREFIX".old.faa "$@"
ret=$?
if ((ret)); then
    echo "Error: Old binary failed to run with $ret."
fi
# rerun if empty (first training creation)
if [[ ! -s "$OUT_PREFIX".old.gbk || ! -s "$OUT_PREFIX".old.ffn || ! -s "$OUT_PREFIX".old.faa ]]; then
    echo "Warning: Old binary produced empty output files, retrying..."
    time "$OLD_BINARY" -i "$INPUT_FILE" -f gbk -o "$OUT_PREFIX".old.gbk -d "$OUT_PREFIX".old.ffn -a "$OUT_PREFIX".old.faa "$@"
    ret=$?
    if ((ret)); then
        echo "Error: Old binary failed to run with $ret on retry."
    fi
fi

time ./prodigal -i "$INPUT_FILE" -f gbk -o "$OUT_PREFIX".new.gbk -d "$OUT_PREFIX".new.ffn -a "$OUT_PREFIX".new.faa "$@"
ret=$?
if ((ret)); then
    echo "Error: New binary failed to run with $ret."
fi

# normalize line endings
sed -i -e 's/\r//g' "$OUT_PREFIX".{old,new}.*

if ! diff -u "$OUT_PREFIX".old.gbk "$OUT_PREFIX".new.gbk >"$OUT_PREFIX".gbk.diff; then
    echo "GBK files differ."
else
    echo "GBK files are identical."
fi

if ! diff -u "$OUT_PREFIX".old.ffn "$OUT_PREFIX".new.ffn >"$OUT_PREFIX".ffn.diff; then
    echo "FFN files differ."
else
    echo "FFN files are identical."
fi

if ! diff -u "$OUT_PREFIX".old.faa "$OUT_PREFIX".new.faa >"$OUT_PREFIX".faa.diff; then
    echo "FAA files differ."
else
    echo "FAA files are identical."
fi
