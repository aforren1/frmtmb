#!/usr/bin/env bash
# Reviewer of lane seflags (check of punch 2b): the merge build op4's
# own test files (and frmtmb.eam's test-defects.R) under both BLAS.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
S=/c/Users/adf44/AppData/Local/Temp/3/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
T=$S/seflags-rev4-tree-op/tests/testthat
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=2
L=dev/seflags-rev4-log/merge; mkdir -p $L
for c in ref ob; do
  R="$RS"; [ $c = ob ] && R="$OB"
  for f in test-se-flags.R test-se-check.R test-se-tier.R test-optrap.R test-optima.R test-diagnostics-ux.R test-autoscale.R test-se-lost-re-predict.R test-brms-postfit.R; do
    "$R" dev/seflags-rev-testfile.R op4 $T/$f < /dev/null > $L/$c-$f.txt 2>&1 &
  done
  "$R" dev/seflags-rev-testfile.R op4 extensions/frmtmb.eam/tests/testthat/test-defects.R frmtmb.eam < /dev/null > $L/$c-eam-test-defects.R.txt 2>&1 &
  wait
done
echo all done
