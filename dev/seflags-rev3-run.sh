#!/usr/bin/env bash
# Reviewer of lane seflags (final check): B4 attacks, illpred, the
# degree-7 prediction, the bitwise check, on lane, base, the lane's
# merge build (m9) and the merge with optrap's tree of 23:49 (op3),
# reference BLAS and OpenBLAS 0.3.26 LAPACK. Waits for every job.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
L=dev/seflags-rev3-log
for a in lane op3 base; do
  "$RS" dev/seflags-rev3-illid.R $a $L/illid-ref-$a.tsv > $L/illid-ref-$a.txt 2>&1 &
done
for a in lane op3; do
  "$OB" dev/seflags-rev3-illid.R $a $L/illid-ob-$a.tsv > $L/illid-ob-$a.txt 2>&1 &
done
wait
for a in lane op3 m9 base; do
  "$RS" dev/seflags-rev3-illpred.R $a > $L/illpred-ref-$a.txt 2>&1 &
  "$OB" dev/seflags-rev3-illpred.R $a > $L/illpred-ob-$a.txt 2>&1 &
done
wait
for a in lane op3 m9 base; do
  "$RS" dev/seflags-rev3-deg7pred.R $a > $L/deg7pred-ref-$a.txt 2>&1 &
  "$OB" dev/seflags-rev3-deg7pred.R $a > $L/deg7pred-ob-$a.txt 2>&1 &
done
"$RS" dev/seflags-rev3-deg7exact.R lane > $L/deg7exact-ref-lane.txt 2>&1 &
"$RS" dev/seflags-rev3-deg7exact.R op3 > $L/deg7exact-ref-op3.txt 2>&1 &
wait
for a in lane base op3; do
  "$RS" dev/seflags-rev3-bitwise.R $a $L/bitwise-ref-$a.rds > $L/bitwise-ref-$a.txt 2>&1 &
done
"$OB" dev/seflags-rev3-bitwise.R lane $L/bitwise-ob-lane.rds > $L/bitwise-ob-lane.txt 2>&1 &
"$OB" dev/seflags-rev3-bitwise.R base $L/bitwise-ob-base.rds > $L/bitwise-ob-base.txt 2>&1 &
wait
echo all done
