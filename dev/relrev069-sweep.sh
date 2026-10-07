#!/usr/bin/env bash
# Reviewer of the 0.69.0 consolidation: dev/relrev069-fuzzsweep.R over
# three BLAS builds (reference; OpenBLAS 0.3.26 and 0.3.32 as BLAS and
# LAPACK, 4 threads), rellib-r7 (0.69.0) and rellib-r6 (0.68.1), and
# rellib-r7 with the candidate fix, for each variant given.
#   bash dev/relrev069-sweep.sh <k from> <k to> <variant> [variant ...]
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
cd "$REL"
k1="$1"; k2="$2"; shift 2
O=dev/relrev069-log/sweep; mkdir -p "$O"
declare -A RS
RS[ref]="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RS[ob26]="$REL/dev/relrev069-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
RS[ob32]="$REL/dev/ciharden-out/Rob0.3.32-lapack/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
for v in "$@"; do
  for b in ref ob26 ob32; do
    for arm in r7 r6 r7fix; do
      case $arm in
        r7) lib=C:/Users/adf44/source/r/rellib-r7; fx="";;
        r6) lib=C:/Users/adf44/source/r/rellib-r6; fx="";;
        r7fix) lib=C:/Users/adf44/source/r/rellib-r7; fx=fix;;
      esac
      "${RS[$b]}" dev/relrev069-fuzzsweep.R base "$lib" "$v" "$k1" "$k2" \
        "$O/$v-$b-$arm.tsv" $fx > "$O/$v-$b-$arm.log" 2>&1 &
    done
  done
done
wait
echo SWEEP DONE
