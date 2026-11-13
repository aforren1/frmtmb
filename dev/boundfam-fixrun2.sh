#!/usr/bin/env bash
# Lane boundfam fix round: OpenBLAS 0.3.26 (LAPACK) on the changed
# numeric test files, then R CMD check --as-cran of core.
set -u
cd "$(dirname "$0")/.."
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/boundfam-out/openblas-build-fix.txt 2>&1
rm -rf dev/boundfam-out/ob-0.3.26
bash dev/boundfam-ob.sh 0.3.26 3 frmtmb:tests/testthat/test-se-limit.R \
  frmtmb:tests/testthat/test-se-tier.R frmtmb:tests/testthat/test-latres.R \
  frmtmb:tests/testthat/test-cens-trunc.R \
  frmtmb:tests/testthat/test-diagnostics-ux.R \
  frmtmb:tests/testthat/test-gp-start.R \
  frmtmb:tests/testthat/test-review-v29.R \
  frmtmb:tests/testthat/test-predfix.R \
  frmtmb.eam:extensions/frmtmb.eam/tests/testthat/test-defects.R \
  frmtmb.learn:extensions/frmtmb.learn/tests/testthat/test-engine.R
rm -rf dev/ciharden-out
bash dev/boundfam-check.sh . > dev/boundfam-out/check-core-fix.txt 2>&1
echo "fixrun2 done"
