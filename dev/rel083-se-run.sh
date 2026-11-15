#!/usr/bin/env bash
# Lane selimit2's standard-error comparison over every core test file of
# the merged tree, on one build: one traced file per R process, N at once
# (dev/rel083-se-traced.R). NOT_CRAN only, as the lane ran it.
#
#   bash dev/rel083-se-run.sh <r22|r21> [N]
#
# Output: dev/rel083-log/se-<arm>/<pkg>--<file>.trace and .txt.
set -u
arm="$1"; N="${2:-3}"
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel083
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win NOT_CRAN=true
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
unset FRMTMB_BRMS_FIT_TESTS FRMTMB_DRMTMB_FIT_TESTS FRMTMB_FUZZ FRMTMB_SCALE_TESTS
case "$arm" in
  r22) export REL083_SE_LIB=C:/Users/adf44/source/r/rellib-r22;;
  r21) export REL083_SE_LIB=C:/Users/adf44/source/r/rellib-r21;;
  *) echo "arm r22 or r21"; exit 1;;
esac
OUT=$ROOT/dev/rel083-log/se-$arm
rm -rf "$OUT"; mkdir -p "$OUT"
port=24100
for f in tests/testthat/test-*.R; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$OUT/frmtmb--$(basename "$f" .R)"
  port=$((port + 1))
  ( SELIMIT2_TRACE="$b.trace" R_PARALLEL_PORT=$port \
      "$RS" dev/rel083-se-traced.R frmtmb "$f" > "$b.txt" 2>&1 ) &
done
wait
echo "SE RUN DONE $arm $(ls $OUT/*.trace | wc -l) traces"
