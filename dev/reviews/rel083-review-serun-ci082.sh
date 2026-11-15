#!/usr/bin/env bash
# (copy of rel083-review-serun.sh for lane ci082, which 8159568e holds)
# Release review of 0.83.0: the 23 core test files whose SE records moved
# (dev/rel083-log/se-cmp.txt), traced on each lane library, N at once.
#   bash dev/reviews/rel083-review-serun.sh [N]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel083
cd "$ROOT" || exit 1
N="${1:-3}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win NOT_CRAN=true
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
unset FRMTMB_BRMS_FIT_TESTS FRMTMB_DRMTMB_FIT_TESTS FRMTMB_FUZZ FRMTMB_SCALE_TESTS
FILES=$(grep 'keyed records differ' dev/rel083-log/se-cmp.txt | sed 's/^frmtmb--//; s/\.trace:.*//' | sort -u)
port=25100
for l in ci082; do
  OUT=$ROOT/dev/reviews/rel083-review-out/se-$l
  rm -rf "$OUT"; mkdir -p "$OUT"
  for f in $FILES; do
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    b="$OUT/frmtmb--$f"; port=$((port + 1))
    ( REL083_SE_LIB=C:/Users/adf44/source/r/wt-$l-lib SELIMIT2_TRACE="$b.trace" R_PARALLEL_PORT=$port \
        "$RS" dev/reviews/rel083-review-setraced.R frmtmb "tests/testthat/$f.R" > "$b.txt" 2>&1 ) &
  done
done
wait
echo "SERUN CI082 DONE"
