#!/usr/bin/env bash
# Punch round 1: the gated tier (FRMTMB_BRMS_FIT_TESTS=true) of every
# frmtmb.sample and frmtmb.spline test file and of core's brms-suite
# files that reach predict(), one R process per file, N at a time.
#   bash dev/gpby-p1-gated.sh [N]
set -u
N="${1:-12}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-gpby
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/gpby-p1-gated"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
rm -rf "$OUT"; mkdir -p "$OUT"
{
  for e in sample spline; do
    for f in "$ROOT"/extensions/frmtmb.$e/tests/testthat/test-*.R; do
      echo "frmtmb.$e $f"
    done
  done
  for f in "$ROOT"/tests/testthat/test-brms-suite-*.R; do
    echo "frmtmb $f"
  done
} > "$OUT/jobs.txt"
cd "$ROOT"
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" dev/gpby-runtest.R "$p" "$f" > "$o" 2>&1 ) &
done < "$OUT/jobs.txt"
wait
for p in frmtmb frmtmb.sample frmtmb.spline; do
  pa=0; fa=0; er=0; sk=0; wa=0; n=0
  for o in "$OUT"/$p--*.txt; do
    r=$(grep -a "^RESULT" "$o" | tail -1)
    n=$((n + 1))
    if [ -z "$r" ]; then echo "  NO RESULT $(basename "$o")"; continue; fi
    pa=$((pa + $(echo "$r" | sed 's/.*pass=\([0-9]*\).*/\1/')))
    fa=$((fa + $(echo "$r" | sed 's/.*fail=\([0-9]*\).*/\1/')))
    er=$((er + $(echo "$r" | sed 's/.*error=\([0-9]*\).*/\1/')))
    sk=$((sk + $(echo "$r" | sed 's/.*skip=\([0-9]*\).*/\1/')))
    wa=$((wa + $(echo "$r" | sed 's/.*warn=\([0-9]*\).*/\1/')))
    case "$r" in *fail=0\ error=0*) ;; *) echo "  BAD $r" ;; esac
  done
  echo "$p files=$n pass=$pa fail=$fa error=$er skip=$sk warn=$wa"
done | tee "$OUT/SUMMARY.txt"
