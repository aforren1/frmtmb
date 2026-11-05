#!/usr/bin/env bash
# The final evidence: every suite, gated, one file per process, on this tree's
# core with rellib-r16 behind it; then the failing form on rellib-r16.
set -u

export FRMTMB_BRMS_FIT_TESTS=true
rm -rf dev/ncgp-out/gated
bash dev/ncgp-par.sh dev/ncgp-out/gated frmtmb frmtmb.sample frmtmb.spatial frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.survival
grep -h '^RESULT' dev/ncgp-out/gated/*.txt > dev/ncgp-out/gated.txt
echo GATED-DONE
# seen failing on the base build
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
NCGP_TLIB="C:/Users/adf44/source/r/rellib-r16" "$RS" dev/ncgp-runtest.R frmtmb tests/testthat/test-gp-noncentered.R > dev/ncgp-out/seen-fail.txt 2>&1
echo SEEN-DONE
