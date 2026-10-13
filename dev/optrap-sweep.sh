#!/usr/bin/env bash
# Lane optrap: dev/optrap-fuzzsweep.R under the reference BLAS and
# OpenBLAS 0.3.26 as BLAS and LAPACK (4 threads), rellib-r7 against the
# lane's library, for each variant given.
#   bash dev/optrap-sweep.sh <k from> <k to> <variant> [variant ...]
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-optrap
cd "$WT"
k1="$1"; k2="$2"; shift 2
O=dev/optrap-log/sweep; mkdir -p "$O"
declare -A RS
RS[ref]="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RS[ob26]="$WT/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
for v in "$@"; do
  for b in ref ob26; do
    for arm in r7 lane; do
      case $arm in
        r7) lib=base;;
        lane) lib=C:/Users/adf44/source/r/wt-optrap-lib;;
      esac
      if [ "$b" = ob26 ]; then
        export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
      else
        unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS
      fi
      "${RS[$b]}" dev/optrap-fuzzsweep.R "$lib" \
        C:/Users/adf44/source/r/rellib-r7 "$v" "$k1" "$k2" \
        "$O/$v-$b-$arm.tsv" > "$O/$v-$b-$arm.log" 2>&1 &
    done
  done
done
wait
echo SWEEP DONE
