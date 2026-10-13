#!/usr/bin/env bash
# Reviewer, final check: run every test file of the named packages, one R
# process per file, N at a time, through dev/release/run-tests.R with
# the lane's library first; logs under dev/moparity-rev3-log/<tier>/.
#   bash dev/moparity-rev3-runtests.sh <tier> <N> <pkg> ...
set -u
tier="$1"; N="$2"; shift 2
ROOT=C:/Users/adf44/source/r/frmtmb-wt-moparity
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/moparity-rev3-log/$tier"
mkdir -p "$OUT"
export TMP="C:/Users/adf44/AppData/Local/Temp/1" TEMP="C:/Users/adf44/AppData/Local/Temp/1"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/moparity-rev3-stan-cache"
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-moparity-lib
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
for p in "$@"; do
  for f in "$ROOT"/extensions/$p/tests/testthat/test-*.R; do
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    ( "$RS" "$ROOT/dev/release/run-tests.R" "$p" "$f" > "$OUT/$p--$(basename "$f" .R).txt" 2>&1 ) &
  done
done
wait
echo DONE
