#!/usr/bin/env bash
# Reviewer of lane optrap: the lane's changed test files on its build
# under OpenBLAS 0.3.26 as BLAS and LAPACK, 4 threads.
cd /c/Users/adf44/source/r/frmtmb-wt-optrap
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-optrap/dev/stan-cache
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
OB=dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe
mkdir -p dev/optrap-rev-log/tf-ob26
for f in optrap grad-verdict se-tier se-check bcm-gaussian glmm-gaussian autocor id-kron brms-suite-methods fuzz; do
  $OB dev/ciharden-run1.R C:/Users/adf44/source/r/wt-optrap-lib frmtmb tests/testthat/test-$f.R > dev/optrap-rev-log/tf-ob26/lane-$f.txt 2>&1 &
done
wait
grep -ah "^RESULT" dev/optrap-rev-log/tf-ob26/*.txt
grep -ah "^blas" dev/optrap-rev-log/tf-ob26/lane-optrap.txt
grep -ah "^DETAIL" dev/optrap-rev-log/tf-ob26/*.txt | grep -v SKIP | cut -c1-200
