#!/bin/sh
# Run test files one per R process, N at a time, against the lane build.
# Usage: sh dev/famuni-sample-suite.sh <tag> <N>
# The job list: every frmtmb.sample test file, gated ones gated.
# Output: dev/famuni-suite-<tag>/<pkg>--<file>.txt, and a summary of the
# RESULT lines in dev/famuni-suite-<tag>/SUMMARY.txt.
tag="$1"; N="${2:-8}"; gated=yes
WT=/c/Users/adf44/source/r/frmtmb-wt-famuni
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$WT/dev/famuni-suite-$tag"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/stan-cache
export NOT_CRAN=true
[ -n "$gated" ] && export FRMTMB_BRMS_FIT_TESTS=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"
jobs_file="$OUT/JOBS.txt"
: > "$jobs_file"
pat='.'
for f in; do
  echo "frmtmb $f" >> "$jobs_file"
done
for e in sample; do
  for f in $(grep -l "$pat" "$WT"/extensions/frmtmb.$e/tests/testthat/test-*.R); do
    echo "frmtmb.$e $f" >> "$jobs_file"
  done
done
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$(dirname "$f")" && "$RS" "$WT/dev/famuni-test1.R" "$p" "$f" > "$o" 2>&1 ) &
done < "$jobs_file"
wait
want=$(wc -l < "$jobs_file")
got=$(grep -l "^RESULT" "$OUT"/*--*.txt | wc -l)
{
  echo "jobs: $want, files with a RESULT line: $got"
  for o in "$OUT"/*--*.txt; do
    r=$(grep "^RESULT" "$o")
    l=$(grep "^lib:" "$o" | cut -c1-70)
    if [ -z "$r" ]; then echo "NORESULT $(basename "$o")"; else echo "$r | $l"; fi
  done
} > "$OUT/SUMMARY.txt"
