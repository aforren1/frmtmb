#!/usr/bin/env bash
# Reviewer batch 1: the two pinning test files on both builds, the
# lane's own comparison scripts on the lane build, and brms's update().
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drop1raw
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/drop1raw-rev-out"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
cd "$ROOT"
for arm in lane base; do
  "$RS" dev/drop1raw-runtest.R $arm frmtmb tests/testthat/test-refit-data.R \
    > "$OUT/t-refit-$arm.txt" 2>&1 &
  "$RS" dev/drop1raw-runtest.R $arm frmtmb.latent \
    extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R \
    > "$OUT/t-hmm-$arm.txt" 2>&1 &
done
"$RS" dev/drop1raw-compare.R lane > "$OUT/compare-lane.txt" 2>&1 &
"$RS" dev/drop1raw-influence.R lane > "$OUT/influence-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-brms2.R > "$OUT/brms2.txt" 2>&1 &
wait
echo done
