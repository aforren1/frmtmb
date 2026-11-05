#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
"$R" dev/review-latres-m17b.R > dev/review-latres-out/m17b.log 2>&1 &
"$R" dev/review-latres-m6b.R > dev/review-latres-out/m6b.log 2>&1 &
"$R" dev/review-latres-m2b.R > dev/review-latres-out/m2b.log 2>&1 &
wait
