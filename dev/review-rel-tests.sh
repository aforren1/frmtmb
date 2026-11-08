#!/usr/bin/env bash
# Release review 0.80.0: run test files one per R process, in order, from
# a job list of lines "<package> <test file> [<mutant lib>]", gated tier.
#   bash dev/review-rel-tests.sh <jobs file> <out dir>
set -u
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$2"
while read -r p f lib; do
  [ -z "$p" ] && continue
  o="$2/$(basename "$f" .R)${lib:+-$(basename "$lib")}.log"
  NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_RELLIB="${lib:-}" \
    "$R" dev/release/run-tests.R "$p" "$f" > "$o" 2>&1
  echo "$(grep -h '^RESULT' "$o" | tail -1) | $(grep -h '^lib:' "$o" | head -1) [$o]"
done < "$1"
echo "TESTS END"
