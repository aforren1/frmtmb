#!/bin/sh
# Reviewer round 2: reruns on the trial merge of gpby and fixes, gated,
# one R process per file, library gpby-rev-lib (merged core, spline,
# sample; other extensions from rellib-r5).
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
T="/c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/gpby-rev2-merge/tree"
cd /c/Users/adf44/source/r/frmtmb-wt-gpby
mkdir -p dev/gpby-rev2-mtests
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
run() { "$R" dev/gpby-rev2-mrun.R $1 "$T/$2/$3" lane > dev/gpby-rev2-mtests/$1-$(basename $3 .R).txt 2>&1; }
C=tests/testthat; P=extensions/frmtmb.spline/tests/testthat; M=extensions/frmtmb.sample/tests/testthat
for f in test-fd-chain.R test-predict-re-uncertainty.R test-emmeans.R test-brms-suite-emmeans.R test-brms-suite-methods.R test-prior-compat.R test-brms-priors.R test-get-prior-route.R test-gp-by.R test-gp-multidim.R test-ce-aterms.R test-ordinal-thres-names.R test-smooths.R test-nl.R test-compat.R test-brms-likelihood.R test-brms-agreement.R test-conditional-smooths.R; do run frmtmb $C $f & done
wait
for f in test-curve.R test-deriv.R test-difference.R test-new-levels.R; do run frmtmb.spline $P $f & done
for f in test-gp-by-draws.R test-incl-thres-draws.R test-draws-spellings.R; do run frmtmb.sample $M $f & done
run frmtmb.eam extensions/frmtmb.eam/tests/testthat test-family.R &
wait
echo ALLDONE
