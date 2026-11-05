#!/usr/bin/env bash
# Every case of dev/nlfix-saddle.R, one R process each, 6 at once. A process
# that dies leaves its log and no row.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
mkdir -p dev/nlfix-saddle-logs
"$R" dev/nlfix-saddle.R list 2>/dev/null | tr -d '\r' > dev/nlfix-saddle-ids.txt
xargs -P 6 -I{} bash -c "timeout 1500 \"$R\" dev/nlfix-saddle.R {} > dev/nlfix-saddle-logs/{}.log 2>&1 || echo FAILED {} exit \$? >> dev/nlfix-saddle-logs/{}.log" < dev/nlfix-saddle-ids.txt
echo done > dev/nlfix-saddle.done
