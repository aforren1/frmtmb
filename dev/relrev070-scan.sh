#!/usr/bin/env bash
# Reviewer of the 0.70.0 consolidation: run a list of test files, one R
# process per file, through dev/ciharden-run1.R, writing only under
# dev/relrev070-log/ (the review may not write elsewhere in the tree).
#
#   bash dev/relrev070-scan.sh <run> <ref|ob26> <base lib> <list> [N]
#
# <list>: "<package> <test file>" lines, paths relative to TREE (default
# the release tree).
set -u
run="$1"; cfg="$2"; base="$3"; list="$4"; N="${5:-12}"
ROOT=C:/Users/adf44/source/r/frmtmb-wt-release
TREE="${TREE:-$ROOT}"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/relrev070-out/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export ROUND_BASE="$base"
unset FRMTMB_SAMPLER_GATES
: "${TMP:?TMP is unset}"
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
       unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS ;;
  ob26) RS="$ROOT/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 ;;
  *) echo "unknown config $cfg"; exit 1 ;;
esac
OUT="$ROOT/dev/relrev070-log/$run"
rm -rf "$OUT"; mkdir -p "$OUT"
cd "$TREE"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$ROOT/dev/ciharden-run1.R" base "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
S="$OUT.sum"; : > "$S"; want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  lib=$(grep -a '^lib:' "$o" | head -1)
  hit=$(grep -a '^RESULT' "$o" | tail -1)
  if [ -n "$hit" ]; then ran=$((ran + 1)); echo "$p $hit | $lib" >> "$S"
  else echo "$p NO RESULT LINE $f" >> "$S"; fi
done < "$list"
echo "ran $ran of $want" >> "$S"
