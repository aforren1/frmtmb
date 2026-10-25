#!/usr/bin/env bash
# Reviewer's runner: one test file per R process, N at a time, against
# the lane library (FRMTMB_RELLIB) or the base library (arm "base").
#   bash dev/review-brmsnames-run.sh <arm lane|base> <tag> <N> <pkg file>...
set -u
arm="$1"; tag="$2"; N="$3"; shift 3
ROOT=/c/Users/adf44/source/r/frmtmb-wt-brmsnames
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/review-brmsnames-log/$tag"
mkdir -p "$OUT"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS="${FRMTMB_BRMS_FIT_TESTS:-true}"
if [ "$arm" = lane ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-brmsnames-lib
else
  export FRMTMB_RELLIB=
fi
cd "$ROOT"
while [ $# -gt 1 ]; do
  p="$1"; f="$2"; shift 2
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  ( "$RS" dev/release/run-tests.R "$p" "$f" > "$OUT/$p--$(basename "$f" .R).txt" 2>&1 ) &
done
wait
for o in "$OUT"/*.txt; do
  echo "$(basename "$o"): $(grep -a '^lib:' "$o" | head -1)"
  grep -a "^RESULT" "$o" | tail -1
done
