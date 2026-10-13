#!/usr/bin/env bash
# Lane seflags, punch round 2b: the reviewer's illpred, illid (p2b copy),
# illfull (p2b copy), deg7exact and deg7pred scripts and the p2b scan
# on one arm, reference BLAS and OpenBLAS 0.3.26 LAPACK. Waits for every
# job it starts.
#   bash dev/seflags-p2b-run.sh <lane|m9|base> <log dir>
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
A=$1
L=$2
mkdir -p "$L"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
for b in ref ob; do
  if [ $b = ref ]; then X="$RS"; else X="$OB"; fi
  "$X" dev/seflags-rev3-illpred.R $A > $L/illpred-$b-$A.txt 2>&1 &
  "$X" dev/seflags-p2b-illid.R $A $L/illid-$b-$A.tsv > $L/illid-$b-$A.txt 2>&1 &
  "$X" dev/seflags-p2b-illfull.R $A > $L/illfull-$b-$A.txt 2>&1 &
  "$X" dev/seflags-rev3-deg7pred.R $A > $L/deg7pred-$b-$A.txt 2>&1 &
done
wait
"$RS" dev/seflags-rev3-deg7exact.R $A > $L/deg7exact-ref-$A.txt 2>&1 &
"$RS" dev/seflags-p2b-scan.R $A > $L/scan-ref-$A.txt 2>&1 &
"$OB" dev/seflags-p2b-scan.R $A > $L/scan-ob-$A.txt 2>&1 &
wait
