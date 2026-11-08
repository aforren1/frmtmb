#!/usr/bin/env bash
# Lane osa: run a job list ("<package> <test file>" per line) through
# dev/osa-test1.R, one R process per file, at most three at once.
#   OSA_ARM=before|after bash dev/osa-par.sh <list> <outdir>
# Writes <outdir>/<package>--<file>.txt per file and <outdir>/summary.log
# with each file's RESULT line, or NO RESULT LINE and its last lines.
set -u
list="$1"; out="$2"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-osa/dev/stan-cache"
export OSA_ARM="${OSA_ARM:-after}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd /c/Users/adf44/source/r/frmtmb-wt-osa || exit 1
mkdir -p "$out"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "${OSA_N:-3}" ]; do sleep 1; done
  o="$out/$p--$(basename "$f" .R).txt"
  ( "$RS" dev/osa-test1.R "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
: > "$out/summary.log"
while read -r p f; do
  [ -n "$p" ] || continue
  o="$out/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  lib=$(grep -a "^lib:" "$o" | head -1)
  if [ -n "$hit" ]; then
    echo "$hit | $lib" >> "$out/summary.log"
  else
    echo "NO RESULT LINE for $f" >> "$out/summary.log"
    tail -5 "$o" >> "$out/summary.log"
  fi
done < "$list"
cat "$out/summary.log"
