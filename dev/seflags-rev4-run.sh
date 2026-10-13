#!/usr/bin/env bash
# Reviewer of lane seflags (check of punch 2b): B4/B5 scripts, the
# routes and whole-covariance scripts, the bitwise set, on lane, base and
# op4, reference BLAS and OpenBLAS 0.3.26 LAPACK. Waits for every job.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
L=dev/seflags-rev4-log
for a in lane op4; do
  for c in ref ob; do
    R="$RS"; [ $c = ob ] && R="$OB"
    "$R" dev/seflags-rev4-illpred.R $a > $L/illpred-$c-$a.txt 2>&1 &
    "$R" dev/seflags-rev4-illid.R $a $L/illid-$c-$a.tsv > $L/illid-$c-$a.txt 2>&1 &
    "$R" dev/seflags-rev4-illfull.R $a > $L/illfull-$c-$a.txt 2>&1 &
  done
  wait
done
"$RS" dev/seflags-rev4-routes.R op4 > $L/routes-ref-op4.txt 2>&1 &
"$OB" dev/seflags-rev4-routes.R lane > $L/routes-ob-lane.txt 2>&1 &
"$OB" dev/seflags-rev4-whole.R lane > $L/whole-ob-lane.txt 2>&1 &
"$RS" dev/seflags-rev4-whole.R op4 > $L/whole-ref-op4.txt 2>&1 &
"$RS" dev/seflags-rev4-deg7pred.R lane > $L/deg7pred-ref-lane.txt 2>&1 &
"$RS" dev/seflags-rev4-deg7pred.R op4 > $L/deg7pred-ref-op4.txt 2>&1 &
wait
for a in lane base; do
  "$RS" dev/seflags-rev4-bitwise.R $a $L/bitwise-ref-$a.rds > $L/bitwise-ref-$a.txt 2>&1 &
  "$OB" dev/seflags-rev4-bitwise.R $a $L/bitwise-ob-$a.rds > $L/bitwise-ob-$a.txt 2>&1 &
done
wait
"$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-ref-lane.rds $L/bitwise-ref-base.rds > $L/bitcmp-ref.txt 2>&1
"$RS" dev/seflags-rev4-bitcmp.R $L/bitwise-ob-lane.rds $L/bitwise-ob-base.rds > $L/bitcmp-ob.txt 2>&1
echo all done
