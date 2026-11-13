#!/usr/bin/env bash
# Lane boundfam, consolidation item 3 (no warning on a NaN prediction
# under OpenBLAS): test-se-limit.R and prediction files under OpenBLAS
# 0.3.26 and 0.3.32 (LAPACK), the core suite on the reference BLAS, and
# R CMD check --as-cran of core.
set -u
cd "$(dirname "$0")/.."
ob() {
  bash dev/boundfam-ob.sh "$1" 3 frmtmb:tests/testthat/test-se-limit.R \
    frmtmb:tests/testthat/test-se-tier.R \
    frmtmb:tests/testthat/test-predict-newdata.R \
    frmtmb:tests/testthat/test-predict-re-uncertainty.R \
    frmtmb:tests/testthat/test-linpred-layers.R \
    frmtmb:tests/testthat/test-predict-lp-basis.R
}
rm -rf dev/boundfam-out/ob-0.3.26 dev/boundfam-out/ob-0.3.32
ob 0.3.26
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/boundfam-out/openblas-build-sep32.txt 2>&1
ob 0.3.32
rm -rf dev/ciharden-out
export NOT_CRAN=true
bash dev/boundfam-par2.sh lane dev/boundfam-jobs-core318.txt 3
rm -rf dev/boundfam-suite6-lane; mv dev/boundfam-suite2-lane dev/boundfam-suite6-lane
unset NOT_CRAN
bash dev/boundfam-check.sh . > dev/boundfam-out/check-core-sep.txt 2>&1
echo "fixrun5 done"
