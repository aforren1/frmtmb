#!/usr/bin/env bash
# Lane nbrobust: the repro grid of dev/nbrobust-repro.R, two R processes
# at a time (one chain per case).
#   bash dev/nbrobust-repro-all.sh <arm> <blas: ref|ob326>
cd "$(dirname "$0")/.."
ARM="$1"; BL="$2"
if [ "$BL" = ref ]; then RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"; CT200=""; else
  RS="dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"; CT200=Sandybridge; fi
export OPENBLAS_NUM_THREADS=1
chain() {
  case="$1"; ct="$2"
  for cse in warning ignore; do
    OPENBLAS_CORETYPE="$ct" "$RS" dev/nbrobust-repro.R "$ARM" "$case" 0 20 "$cse" \
      > "dev/nbrobust-out/$ARM-$BL-$case-$cse.txt" 2>&1
  done
}
chain nb200 "$CT200" &
chain nb300 "" &
wait
