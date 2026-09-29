#!/usr/bin/env bash
# Every test file the change can reach, ONE per R process, up to $JOBS
# at a time. Results land in dev/thresrefit-testlog/.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-thresrefit
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/thresrefit-testlog"
JOBS=${JOBS:-7}
mkdir -p "$OUT"

FILES="
test-thres-refit.R
test-thres.R
test-sratio-thresholds.R
test-ordinal.R
test-ordinal-fitted.R
test-boot.R
test-influence-plot.R
test-multiple-pooling.R
test-simulate-ergonomics.R
test-simulate-density.R
test-simulate-newdata.R
test-autoscale.R
test-frame.R
test-unpinned-seams.R
test-data2.R
test-importance.R
test-methods-audit.R
test-numerical-robustness.R
test-review-fixes.R
test-v14.R
test-v15.R
test-verbose.R
test-setprior.R
test-brms-shapes.R
test-brms-shapes-punch.R
test-ce-bands.R
test-compat.R
test-bracket-access.R
test-structure.R
test-predfix.R
test-api-spellings.R
test-arg-refusal.R
test-diagnostics-ux.R
test-input-validation.R
test-spectral.R
test-brms-families.R
test-review-v29.R
test-osa-inference.R
test-case-studies.R
test-edgecases.R
test-v07.R
test-tabular-inputs.R
test-mv-gaps.R
test-brms-agreement.R
"

n=0
for f in $FILES; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  n=$((n + 1))
  ( cd "$ROOT" && "$RS" dev/thresrefit-runtest.R "$f" \
      > "$OUT/$f.log" 2>&1 ) &
done
wait
echo "launched $n files"
grep -h "^RESULT" "$OUT"/*.log | sort
echo "--- files with no RESULT line (aborted) ---"
for f in $FILES; do
  grep -q "^RESULT" "$OUT/$f.log" 2>/dev/null || echo "$f"
done
