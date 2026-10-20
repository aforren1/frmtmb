#!/usr/bin/env bash
# Review of lane illse2, item 6: the false-alarm battery on fresh seeds
# (101:115, outside the lane's 1:100), reference BLAS, lane against
# rellib-r9, at most N R processes at once.
#   bash dev/illse2-rev-sf.sh
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-illse2
cd "$W"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
L=dev/illse2-rev-log/sf
mkdir -p $L
N=${N:-4}
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
for g in nl lin poly mixed bounded; do
  for a in lane base; do
    jobs_wait
    "$RS" dev/illse2-sf-rev3-falarm.R $a $g $L/$g-ref-$a.tsv "101:115" \
      < /dev/null > $L/$g-ref-$a.txt 2>&1 &
  done
done
wait
ILLSE_CFGS="ref" "$RS" dev/illse2-sf-falarm-cmp.R $L lane $L base > $L/cmp.txt 2>&1
echo done > $L/done.txt
