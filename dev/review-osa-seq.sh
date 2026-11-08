#!/usr/bin/env bash
# Review lane osa: run test files one after another, one R process each.
#   REV_OSA_LIB=<lib> bash dev/review-osa-seq.sh <outdir> <pkg> <file>...
set -u
out="$1"; pkg="$2"; shift 2
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-osa/dev/stan-cache"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd /c/Users/adf44/source/r/frmtmb-wt-osa || exit 1
mkdir -p "$out"
for f in "$@"; do
  o="$out/$pkg--$(basename "$f" .R).txt"
  "$RS" dev/review-osa-test1.R "$pkg" "$f" > "$o" 2>&1
  echo "exit $?" >> "$o"
  grep -a "^lib:\|^RESULT" "$o"
done
