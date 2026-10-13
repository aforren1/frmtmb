#!/usr/bin/env bash
# Reviewer of lane optrap, item 2: the lane's 300-seed fuzz sweeps
# (dev/optrap-fuzzsweep.R) rerun on the lane's final build, reference
# BLAS and OpenBLAS 0.3.26 as BLAS and LAPACK (4 threads).
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-optrap
cd "$WT"
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
O=dev/optrap-rev-log/sweep; mkdir -p "$O"
declare -A RS
RS[ref]="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RS[ob26]="$WT/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
for v in spec mo; do
  for b in ref ob26; do
    if [ "$b" = ob26 ]; then
      export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
    else
      unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS
    fi
    "${RS[$b]}" dev/optrap-fuzzsweep.R C:/Users/adf44/source/r/wt-optrap-lib \
      C:/Users/adf44/source/r/rellib-r7 "$v" 1 300 \
      "$O/$v-$b-lane.tsv" > "$O/$v-$b-lane.log" 2>&1 &
  done
done
wait
echo SWEEP DONE
