#!/usr/bin/env bash
# Lane seflags, punch round 2b: the merge build's test files (the merge
# tree's own copies, optrap's changed tests included), both BLAS builds,
# one file per process. Waits for every job.
#   bash dev/seflags-p2b-merge-tests.sh
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
T=dev/seflags-merge-out/src/tests/testthat
L=dev/seflags-log/p2b/merge
mkdir -p $L
for f in se-flags se-check se-tier optrap optima diagnostics-ux autoscale \
         se-lost-re-predict brms-postfit; do
  "$RS" dev/seflags-testfile.R merge $T/test-$f.R > $L/$f-ref.txt 2>&1 &
  OPENBLAS_NUM_THREADS=2 "$OB" dev/seflags-testfile.R merge $T/test-$f.R \
    > $L/$f-ob.txt 2>&1 &
done
E=extensions/frmtmb.eam/tests/testthat/test-defects.R
"$RS" dev/seflags-testfile.R merge $E frmtmb.eam > $L/eam-defects-ref.txt 2>&1 &
OPENBLAS_NUM_THREADS=2 "$OB" dev/seflags-testfile.R merge $E frmtmb.eam \
  > $L/eam-defects-ob.txt 2>&1 &
wait
grep -H RESULT $L/*.txt
