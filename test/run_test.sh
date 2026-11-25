#!/usr/bin/env bash
# Runs dms_data_lint on the bundled test data and compares exit status and message with the
# expected ones (those of fq lint 0.12.0 on the same files; read errors by kind).
# Usage (from the repository root, after building): bash test/run_test.sh [path/to/dms_data_lint]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="${1:-$HERE/../dms_data_lint}"
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
cd "$HERE"
"$BIN" --version
fail=0

# 5,000 read pairs of the nf-core/deepmutscan test data: passes, highest base quality Q40 in both files
mkdir -p out
"$BIN" --max-base-qual out/max_base_qual.tsv GID1A_test_1.fastq.gz GID1A_test_2.fastq.gz
if [ "$(cat out/max_base_qual.tsv)" != "$(printf '40\t40')" ]; then echo "FAIL max-base-qual: $(cat out/max_base_qual.tsv)"; fail=1; fi

# the small cases, one pair each
cd cases
n=0
while IFS=$'\t' read -r name want_rc want_msg; do
    [ "$name" = case ] && continue
    stem=${name%%.*}; ext=.${name#*.}
    rc=0
    msg=$("$BIN" "${stem}_1${ext}" "${stem}_2${ext}" 2>&1 >/dev/null | grep -v '^dms_data_lint: ' | head -n 1) || true
    "$BIN" "${stem}_1${ext}" "${stem}_2${ext}" >/dev/null 2>&1 || rc=$?
    case "$msg" in "Error: I/O error"*) msg="Error: I/O error";; esac
    if [ "$rc" != "$want_rc" ] || [ "$msg" != "$want_msg" ]; then
        echo "FAIL $name: exit $rc, '$msg' (expected exit $want_rc, '$want_msg')"; fail=1
    fi
    n=$((n + 1))
done < expected.tsv
cd ..
if [ "$fail" -ne 0 ]; then exit 1; fi
echo "OK: the test pair and all $n cases as expected"
