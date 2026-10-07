#!/usr/bin/env bash
# Reviewer of lane setier: rerun (re-check) lane nanse review attack set on the
# setier build (ref and OpenBLAS) and on rellib-r6 (ref).
#   bash dev/setier-rev-attacks.sh
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-setier
cd "$W"
L=dev/setier-rev2-log/attacks
mkdir -p "$L"
REF="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$W/dev/setier-out/Rob/bin/x64/Rscript.exe"
export OPENBLAS_NUM_THREADS=4
run() { # rscript tag script args...
  local rs="$1" tag="$2" s="$3"; shift 3
  "$rs" "$s" "$@" > "$L/$(basename "$s" .R)-$tag.txt" 2>&1
}
for s in dev/setier-rev-a-rev-cases.R dev/setier-rev-a-rev-spread2.R \
         dev/setier-rev-a-rev-grby.R dev/setier-rev-a-rev-pred.R; do
  run "$REF" lane-ref "$s" lane &
  run "$OB" lane-ob "$s" lane &
  run "$REF" base-ref "$s" base &
done
for s in dev/setier-rev-a-rev2-b1.R dev/setier-rev-a-rev2-b1-rho.R \
         dev/setier-rev-a-rev3-falseloss-sweep.R; do
  run "$REF" lane-ref "$s" lane &
  run "$OB" lane-ob "$s" lane &
  run "$REF" base-ref "$s" merge &
done
run "$REF" lane-ref dev/setier-rev-a-rev2-ridge-re.R lane &
run "$OB" lane-ob dev/setier-rev-a-rev2-ridge-re.R lane &
run "$REF" base-ref dev/setier-rev-a-rev2-ridge-re.R base &
wait
echo done
