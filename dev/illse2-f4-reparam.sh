#!/usr/bin/env bash
# Lane illse2, resumed: test-se-reparam.R (its refusal test now on the
# A = 0.05 body) under the three BLAS, and with the guard off.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-illse2
LIB=C:/Users/adf44/source/r/wt-illse2-lib
J=dev/illse2-log/f4-reparam.jobs
bash dev/illse2-scan.sh f4b-ref ref $LIB 1 $J &
bash dev/illse2-scan.sh f4b-ob26 ob0.3.26-lapack-t4 $LIB 1 $J &
bash dev/illse2-scan.sh f4b-ob32 ob0.3.32-lapack-t4 $LIB 1 $J &
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/illse2-fx-mock.R noguard > dev/illse2-log/f4/mock-noguard-b.txt 2>&1 &
wait
