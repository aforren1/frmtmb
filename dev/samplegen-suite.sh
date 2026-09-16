#!/usr/bin/env bash
# One test file per process, WORKERS at a time (default 2).
#   WORKERS=3 bash dev/samplegen-suite.sh <runner.R> <outdir> <file>...
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd /c/Users/adf44/source/r/frmtmb-wt-samplegen || exit 1
RUNNER="$1"; OUT="$2"; shift 2
W="${WORKERS:-2}"
mkdir -p "$OUT"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/stan-cache"
for f in "$@"; do
  while [ "$(jobs -rp | wc -l)" -ge "$W" ]; do sleep 1; done
  b=$(basename "$f" .R)
  "$RS" "$RUNNER" "$f" > "$OUT/$b.txt" 2>&1 &
done
wait
echo SUITEDONE
