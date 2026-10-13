#!/usr/bin/env bash
# Reviewer of lane optrap, overlap with lane seflags: seflags' SE attack
# set (dev/seflags-attacks.sh, its "lane" arm only) on the trial merge of
# optrap and seflags (library optrap-rev-lib), reference BLAS. The
# scripts are seflags' retargeted copies with the lane library replaced
# by optrap-rev-lib (dev/optrap-rev-att/scripts/). Compare with seflags'
# own lane logs in frmtmb-wt-seflags/dev/seflags-log/att-ref3/.
#   bash dev/optrap-rev-att.sh [N]
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-optrap
cd "$W"
N="${1:-10}"
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
export NOT_CRAN=true
unset OPENBLAS_NUM_THREADS OMP_NUM_THREADS
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LL=C:/Users/adf44/source/r/optrap-rev-lib
L="$W/dev/optrap-rev-log/att"
S="$W/dev/optrap-rev-att/scripts"
mkdir -p "$L"
export SEFLAGS_OUT="$L"
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
"$RS" "$S/setier-singular.R" "$LL" "$L/singular-lane.tsv" \
  > "$L/singular-lane.txt" 2>&1 &
"$RS" "$S/setier-rev-smallsd.R" "$LL" "$L/smallsd-lane.tsv" 40 \
  > "$L/smallsd-lane.txt" 2>&1 &
"$RS" "$S/setier-sep.R" "$LL" "$L/sep-lane.tsv" > "$L/sep-lane.txt" 2>&1 &
for s in setier-rev-a-rev-cases setier-rev-a-rev-spread2 \
         setier-rev-a-rev-grby setier-rev-a-rev-pred setier-rev-a-rev2-b1 \
         setier-rev-a-rev2-b1-rho setier-rev-a-rev3-falseloss-sweep \
         setier-rev-a-rev2-ridge-re setier-rev2-window setier-rev-probe \
         setier-rev-c0k setier-rev-c0k-sweep setier-rev-aexpb \
         setier-rev-deg7 setier-rev2-bnd setier-rev-smooth2 \
         setier-rev-fixes-table; do
  extra=""
  [ "$s" = setier-rev-deg7 ] && extra=7
  arm=lane
  [ "$s" = setier-rev-fixes-table ] && arm=merge
  jobs_wait
  "$RS" "$S/$s.R" "$arm" $extra > "$L/$s-lane.txt" 2>&1 &
done
for s in setier-rev2-trap setier-rev2-short setier-rev2-scale2 \
         setier-rev2-item8 setier-rev-rs20 setier-rev3-trapq setier-rev3-quad; do
  jobs_wait
  "$RS" "$S/$s.R" "$LL" > "$L/$s-lane.txt" 2>&1 &
done
jobs_wait
"$RS" "$S/setier-rev2-sepmiss.R" "$LL" budget few > "$L/sepmiss-lane.txt" 2>&1 &
jobs_wait
"$RS" "$S/setier-rev2-sepmiss.R" "$LL" bigp > "$L/sepmiss-bigp-lane.txt" 2>&1 &
wait
echo ATTACKS DONE
