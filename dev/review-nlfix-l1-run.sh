#!/usr/bin/env bash
# bash dev/review-nlfix-l1-run.sh <spec> ...: each spec, one R process, 5 at once.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
out=dev/review-nlfix-l1-out
mkdir -p "$out/logs"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
printf '%s\n' "$@" | xargs -P 5 -I{} bash -c \
  "timeout 1200 \"$R\" dev/review-nlfix-l1.R {} > \"$out/logs/{}.log\" 2>&1 || echo FAILED {} >> \"$out/logs/{}.log\""
for s in "$@"; do
  if [ -f "$out/$s.txt" ]; then cat "$out/$s.txt"; echo; else echo "$s: NO RESULT"; tail -3 "$out/logs/$s.log" | cut -c1-240; fi
done
