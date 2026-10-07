#!/usr/bin/env bash
# Fragility scan: run every test file of all eight suites, one R process
# per file, every gate on, under one BLAS and thread configuration.
# Run it under several configurations and compare them with
# dev/ciharden-scan-diff.R; a test whose result moves with the
# configuration depends on platform rounding.
#
#   bash dev/ciharden-scan.sh <run> <config> <lib dir|base> [N] [list]
#
# <config> is "ref" (this R, reference BLAS and LAPACK, as on the
# Windows runner) or "ob<version>-<blas|lapack>-t<threads>", for
# example ob0.3.32-lapack-t4 (ubuntu-26.04's OpenBLAS for BLAS and
# LAPACK, 4 threads as on the 4-core runner). Build the R copy first:
#   bash dev/ciharden-openblas.sh 0.3.32 lapack
# <lib dir|base>: a library put ahead of rellib-r6, or "base" for
# rellib-r6 alone. Point ROUND_BASE at another base library to rerun
# the scan on a merged tree (see dev/ciharden-run1.R). Point TREE at
# another checkout to run ITS test files, for example an export of the
# base commit (git archive HEAD | tar -x -C <dir>, plus dev/brms-suite)
# for a "before" arm whose tests are the base's own.
# [N] processes at once (default 16). [list]: "<package> <test file>"
# lines; default every test-*.R of core and the seven extensions.
#
# Output: dev/ciharden-log/<run>-suite/<pkg>--<file>.txt per file and
# dev/ciharden-log/<run>.sum, one RESULT line per file or NO RESULT
# LINE, then "ran X of Y". Never tail a per-file log while it runs.
set -u
run="$1"; cfg="$2"; lib="$3"; N="${4:-16}"; list="${5:-}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TREE="${TREE:-$ROOT}"
cd "$TREE"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
export FRMTMB_SCALE_TESTS=true
unset FRMTMB_SAMPLER_GATES
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
case "$cfg" in
  ref)
    RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
    unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS ;;
  ob*-t*)
    d="${cfg%-t*}"; d="Rob${d#ob}"; t="${cfg##*-t}"
    RS="$ROOT/dev/ciharden-out/$d/bin/x64/Rscript.exe"
    export OPENBLAS_NUM_THREADS="$t" OMP_NUM_THREADS="$t" ;;
  *) echo "unknown config $cfg"; exit 1 ;;
esac
[ -x "$RS" ] || { echo "no Rscript at $RS; build it first"; exit 1; }
LOGD="$ROOT/dev/ciharden-rev-log"
mkdir -p "$LOGD"
if [ -z "$list" ]; then
  list="$LOGD/$run.jobs"
  : > "$list"
  for f in tests/testthat/test-*.R; do echo "frmtmb $f" >> "$list"; done
  for e in coupling eam latent learn ode sample spline; do
    for f in extensions/frmtmb.$e/tests/testthat/test-*.R; do
      echo "frmtmb.$e $f" >> "$list"
    done
  done
fi
OUT="$LOGD/$run-suite"
SUM="$LOGD/$run.sum"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$SUM"
# longest files first when an earlier run left timings, so the slow
# fits do not start last
order="$list"
if [ -f "$LOGD/timings.tsv" ]; then
  order="$LOGD/$run.order"
  awk 'NR==FNR { t[$1 " " $2] = $3; next }
       { k = $1 " " $2; print ((k in t) ? t[k] : 999999) "\t" $0 }' \
    "$LOGD/timings.tsv" "$list" | sort -t "$(printf '\t')" -k1,1nr |
    cut -f2- > "$order"
fi
echo "config $cfg Rscript $RS threads ${OPENBLAS_NUM_THREADS:-default} tree $TREE lib $lib" \
  > "$OUT/CONFIG"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( s=$(date +%s); "$RS" "$ROOT/dev/ciharden-run1.R" "$lib" "$p" "$f" > "$o" 2>&1
    echo "ELAPSED $(( $(date +%s) - s ))" >> "$o" ) &
done < "$order"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -n 1)
  if [ -n "$hit" ] && ! echo "$hit" | grep -q LOADERROR; then
    ran=$((ran + 1)); echo "$p $hit" >> "$SUM"
  else
    echo "$p NO RESULT LINE for $f" >> "$SUM"; tail -n 5 "$o" >> "$SUM"
  fi
done < "$list"
echo "ran $ran of $want" >> "$SUM"
tail -n 1 "$SUM"
