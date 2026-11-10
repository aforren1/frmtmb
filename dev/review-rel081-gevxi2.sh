#!/usr/bin/env bash
# Release review of 0.81.0: dev/review-rel081-gevxi.R on small n and on
# a censored scalar xi (blocked from the tmp_xi chart), 3 at a time.
#   bash dev/review-rel081-gevxi2.sh <lib>
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
for s in 1 2 3; do
  echo "cens -0.3 100 $s"; echo "cens -0.6 30 $s"; echo "xi1 -0.6 30 $s"
done | xargs -P 3 -I{} bash -c 'run "{}"'
echo GRID2 DONE
