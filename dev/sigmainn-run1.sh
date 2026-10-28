#!/bin/bash
# One test file in one R process, on the lane lib ("after") or on
# rellib-r14 ("before"). Extra environment (FRMTMB_BRMS_FIT_TESTS) is
# passed through. Usage: bash dev/sigmainn-run1.sh <arm> <pkg> <file> <out>
arm=$1; p=$2; f=$3; o=$4
WT=C:/Users/adf44/source/r/frmtmb-wt-sigmainn
export NOT_CRAN=true R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
if [ "$arm" = after ]; then export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-sigmainn-lib; else unset FRMTMB_RELLIB; fi
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" $WT/dev/release/run-tests.R "$p" "$f" > "$o" 2>&1
