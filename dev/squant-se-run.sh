#!/bin/sh
# Lane squant: run dev/squant-se.R over its designs, two R processes at
# a time. Run: sh dev/squant-se-run.sh
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd C:/Users/adf44/source/r/frmtmb-wt-squant
"$R" dev/squant-se.R homo 300 400 > dev/squant-log/se-homo-400.log 2>&1 &
"$R" dev/squant-se.R hetero 300 400 > dev/squant-log/se-hetero-400.log 2>&1
wait
"$R" dev/squant-se.R skew 300 400 > dev/squant-log/se-skew-400.log 2>&1 &
"$R" dev/squant-se.R homo 300 100 > dev/squant-log/se-homo-100.log 2>&1
wait
"$R" dev/squant-se.R homo 200 2000 > dev/squant-log/se-homo-2000.log 2>&1
