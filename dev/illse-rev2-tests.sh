#!/usr/bin/env bash
# Reviewer of lane illse, round 2: the SE and optimizer test files on
# install 13, reference BLAS and OpenBLAS 0.3.32, one file per process,
# at most 4 at once.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-illse
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe"
export NOT_CRAN=true R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
L=dev/illse-rev2-log/tests
mkdir -p $L
N=4
jw() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
for f in test-se-flags.R test-se-check.R test-se-tier.R test-optrap.R test-grad-verdict.R test-se-reparam.R test-se-conv.R test-optima.R; do
  for c in ref ob; do
    jw
    if [ $c = ob ]; then OPENBLAS_NUM_THREADS=2 "$OB" dev/illse-rev-tf.R tests/testthat/$f > $L/${f%.R}-$c.txt 2>&1 &
    else "$RS" dev/illse-rev-tf.R tests/testthat/$f > $L/${f%.R}-$c.txt 2>&1 & fi
  done
done
wait
grep -h RESULT $L/*.txt > $L/summary.txt
