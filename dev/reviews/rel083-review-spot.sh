#!/usr/bin/env bash
# Release review of 0.83.0: spot reruns of test files on rellib-r22, every
# gate on, one R process per file, N at once. Output in
# dev/reviews/rel083-review-out/<run>/ (gitignored, *-out/).
#   bash dev/reviews/rel083-review-spot.sh <run> <ref|ob0.3.32-lapack-t4> <N> <list>
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel083
cd "$ROOT" || exit 1
run="$1"; cfg="$2"; N="$3"; list="$4"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r22
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
unset FRMTMB_SAMPLER_GATES
case "$cfg" in
  ref) RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"; unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS ;;
  ob*-t*) d="${cfg%-t*}"; d="Rob${d#ob}"; t="${cfg##*-t}"
    RS="$ROOT/dev/ciharden-out/$d/bin/x64/Rscript.exe"; export OPENBLAS_NUM_THREADS="$t" OMP_NUM_THREADS="$t" ;;
esac
OUT=$ROOT/dev/reviews/rel083-review-out/$run
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$ROOT/dev/ciharden-run1.R" base "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
grep -h '^RESULT\|^lib:\|^blas:' "$OUT"/*.txt > "$OUT/SUM.txt"
echo "SPOT END $run"
