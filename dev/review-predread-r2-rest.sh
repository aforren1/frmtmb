#!/usr/bin/env bash
# Short check of fix round 2: the lane's claims read again on its build, one process at a
# time, each log under r2/. mefam writes a fixed name, so it is moved in by a copy.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
O=dev/review-predread-out/r2
export REV_LIB="C:/Users/adf44/source/r/wt-predread-lib"
"$RS" dev/predread-r1-mutant.R C:/Users/adf44/source/r/wt-predread-lib > "$O/r1-mutant.log" 2>&1
"$RS" dev/review-predread-r2-meslope.R > "$O/meslope.log" 2>&1
"$RS" dev/review-predread-r1-check1.R > "$O/check1.log" 2>&1
"$RS" dev/review-predread-r1-check3.R > "$O/check3.log" 2>&1
"$RS" dev/review-predread-r1-check5.R > "$O/check5.log" 2>&1
"$RS" dev/review-predread-mefam.R > "$O/mefam.log" 2>&1
cp dev/review-predread-out/mefam.tsv "$O/mefam.tsv"
echo finished > "$O/rest.done"
