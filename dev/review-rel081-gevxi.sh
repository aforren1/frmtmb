#!/usr/bin/env bash
# Release review of 0.81.0: dev/review-rel081-gevxi.R over a grid, 3 R
# processes at a time. Logs dev/review-rel081-gevxi-<case>_<xi0>_<n>_<seed>.txt
#   bash dev/review-rel081-gevxi.sh <lib>
cd /c/Users/adf44/source/r/frmtmb-wt-rel081
. dev/rel081-env.sh
export NOT_CRAN=true
LIB=${1:-C:/Users/adf44/source/r/rellib-r20}
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
run() {
  set -- $1
  "$RS" dev/review-rel081-gevxi.R "$LIB" "$1" "$4" "$2" "$3" \
    > "dev/review-rel081-gevxi-$1_$2_$3_$4.txt" 2>&1
}
export -f run
export LIB RS
for c in xix xi1 sc; do
  for xi in -0.3 0.2; do
    for s in 1 2 3; do echo "$c $xi 100 $s"; done
  done
done | xargs -P 3 -I{} bash -c 'run "{}"'
echo GRID DONE
