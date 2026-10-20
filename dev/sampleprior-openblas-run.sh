#!/usr/bin/env bash
# Lane sampleprior: the prior-density test files under OpenBLAS 0.3.26
# with LAPACK (dev/ciharden-openblas.sh), one file per process.
#   bash dev/sampleprior-openblas-run.sh <out dir name>
WT=/c/Users/adf44/source/r/frmtmb-wt-sampleprior
cd "$WT"
export OPENBLAS_NUM_THREADS=4 FRMTMB_BRMS_FIT_TESTS=true
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}" TEMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}"
RS=dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe
out="dev/sampleprior-suite/$1"
mkdir -p "$out"
for job in "frmtmb tests/testthat/test-prior-support.R" \
           "frmtmb tests/testthat/test-prior-compat.R" \
           "frmtmb tests/testthat/test-brms-priors.R" \
           "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-default-priors-dpar.R" \
           "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-default-priors-bounds.R" \
           "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-default-priors-brms.R" \
           "frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-draws-sd-names.R"; do
  set -- $job
  "$RS" dev/sampleprior-runtest.R lane "$1" "$2" > "$out/$1--$(basename "$2" .R).txt" 2>&1
  grep -h "^RESULT" "$out/$1--$(basename "$2" .R).txt"
done
