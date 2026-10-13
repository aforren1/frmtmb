#!/usr/bin/env bash
# Reviewer of lane seflags (re-check): the new attacks (alias variants
# under OpenBLAS, the span threshold scan, illcond, QR cost) and the
# curved-ridge reruns. Waits for every job before it exits.
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
L=dev/seflags-rev2-log/B
mkdir -p $L
for a in lane base op2; do
  "$OB" dev/seflags-rev2-alias.R $a > $L/alias-ob-$a.txt 2>&1 &
done
"$RS" dev/seflags-rev2-alias.R op2 > $L/alias-ref-op2.txt 2>&1 &
for a in lane op2; do
  "$RS" dev/seflags-rev2-thresh.R $a $L/thresh-ref-$a.tsv > $L/thresh-ref-$a.txt 2>&1 &
  "$RS" dev/seflags-rev2-ill.R $a > $L/ill-ref-$a.txt 2>&1 &
  "$RS" dev/seflags-rev2-aexpb.R $a > $L/aexpb-ref-$a.txt 2>&1 &
  "$RS" dev/seflags-rev2-c0k-sweep.R $a > $L/c0k-ref-$a.txt 2>&1 &
done
"$OB" dev/seflags-rev2-thresh.R lane $L/thresh-ob-lane.tsv > $L/thresh-ob-lane.txt 2>&1 &
"$OB" dev/seflags-rev2-ill.R lane > $L/ill-ob-lane.txt 2>&1 &
"$OB" dev/seflags-rev2-aexpb.R lane > $L/aexpb-ob-lane.txt 2>&1 &
"$OB" dev/seflags-rev2-c0k-sweep.R lane > $L/c0k-ob-lane.txt 2>&1 &
wait
for a in lane base; do
  "$RS" dev/seflags-rev2-qrcost.R $a > $L/qrcost-ref-$a.txt 2>&1 &
done
wait
echo all done
