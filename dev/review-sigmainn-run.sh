#!/bin/bash
# Reviewer's driver: one test file per R process, N at a time, on the
# lane library (FRMTMB_RELLIB) or on rellib-r14 (arm "before").
# Usage: bash dev/review-sigmainn-run.sh <arm> <jobs> <tag> [N] [gated]
arm=$1; jobs=$2; tag=$3; N=${4:-4}; gated=${5:-}
WT=C:/Users/adf44/source/r/frmtmb-wt-sigmainn
O=$WT/dev/review-sigmainn-out/$tag
mkdir -p "$O"
export NOT_CRAN=true R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
[ -n "$gated" ] && export FRMTMB_BRMS_FIT_TESTS=true
if [ "$arm" = after ]; then export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-sigmainn-lib; else unset FRMTMB_RELLIB; fi
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  ( cd $WT && "$RS" $WT/dev/release/run-tests.R "$p" "$f" > "$O/$p--$(basename "$f" .R).txt" 2>&1 ) &
done < "$jobs"
wait
grep -H "^lib:\|^RESULT" "$O"/*.txt
