#!/usr/bin/env bash
# Lane illse2, fix round 1 resumed: the 44 files that read the
# covariance, the SE check, the optimizer's stop, REML, nonlinear
# bodies or bounds, under the reference BLAS and OpenBLAS 0.3.26 and
# 0.3.32 LAPACK, one after another, on the final install.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-illse2
LIB=C:/Users/adf44/source/r/wt-illse2-lib
J=dev/illse2-log/ob-subset.jobs
bash dev/illse2-scan.sh f4-ref ref $LIB 6 $J
bash dev/illse2-scan.sh f4-ob26 ob0.3.26-lapack-t4 $LIB 6 $J
bash dev/illse2-scan.sh f4-ob32 ob0.3.32-lapack-t4 $LIB 6 $J
