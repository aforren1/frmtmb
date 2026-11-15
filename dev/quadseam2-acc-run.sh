#!/usr/bin/env bash
# The 120-fit accuracy table and the two tight sweeps of 0.82.0 on the
# lane build (base logs: frmtmb-wt-rel082/dev/rel082-log/kt/quadseam/).
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
LIB=C:/Users/adf44/source/r/wt-quadseam2-lib
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$WT"
O=dev/quadseam2-log
"$RS" dev/quadseam2-acc3.R lane bern500,pois200 6:35 > $O/acc3-lane-a.txt 2>&1 &
"$RS" dev/quadseam2-acc3.R lane bern2sd3,poissd2 6:35 > $O/acc3-lane-b.txt 2>&1 &
"$RS" dev/quadseam2-sweep.R $LIB bern5,bern500,bern20big,bern2sd3,bin500,bin50,binz,pois8,pois200,poissd2 1:5 1e-4,1e-5 > $O/sweep-lib.txt 2>&1 &
wait
echo done
