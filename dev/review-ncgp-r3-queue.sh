#!/usr/bin/env bash
# Short check: the runs after the fresh tables, at most 5 at once, cost last and alone.
cd /c/Users/adf44/source/r/frmtmb-wt-ncgp
O=dev/review-ncgp-out
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export TMP="C:/Users/adf44/AppData/Local/Temp" TEMP="C:/Users/adf44/AppData/Local/Temp" R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
until grep -q DONE $O/r3-gauss-lane.txt && grep -q DONE $O/r3-r2fresh-lane.txt && grep -q DONE $O/r3-gauss-base.txt && grep -q DONE $O/r3-r2fresh-base.txt; do sleep 20; done
"$R" dev/ncgp-r2arm.R > $O/r3-arm.txt 2>&1 &
"$R" dev/ncgp-polish.R > $O/r3-polish.txt 2>&1 &
"$R" dev/review-ncgp-r2-claim.R > $O/r3-claim.txt 2>&1 &
"$R" dev/review-ncgp-r2-start.R 7702 > $O/r3-start.txt 2>&1 &
wait
"$R" dev/review-ncgp-nfpred.R > $O/r3-nfpred.txt 2>&1 &
"$R" dev/review-ncgp-lanerows.R > $O/r3-lanerows.txt 2>&1 &
for m in control nug11 nug13 sd1e4 sd1e2 rho1e4 flip zdens; do "$R" dev/review-ncgp-r3-row10.R $m > $O/r3-row10-$m.txt 2>&1 & while [ "$(jobs -r | wc -l)" -ge 5 ]; do sleep 5; done; done
wait
for m in control nopolish quiet_headroom nowarm noorder norowcheck; do "$R" dev/ncgp-mut.R $m 2>&1 | grep "^NCGPM" | cut -c1-170 >> $O/r3-mut.txt & while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 5; done; done
wait
"$R" dev/ncgp-readers.R > $O/r3-readers-lane.txt 2>&1 &
powershell.exe -NoProfile -File dev/review-ncgp-regress.ps1 > $O/regress-driver.txt 2>&1 &
wait
bash dev/review-ncgp-check.sh . frmtmb > $O/r3-check-core.txt 2>&1 &
bash dev/review-ncgp-check.sh extensions/frmtmb.sample frmtmb.sample > $O/r3-check-sample.txt 2>&1 &
wait
"$R" dev/review-ncgp-cost.R base 500 3 2>&1 | grep "^RCO" > $O/r3-cost.txt
"$R" dev/review-ncgp-cost.R lane 500 3 2>&1 | grep "^RCO" >> $O/r3-cost.txt
"$R" dev/review-ncgp-mem.R lane 300 2 > $O/r3-mem-lane.txt 2>&1
echo QUEUEDONE > $O/r3-queue.done
