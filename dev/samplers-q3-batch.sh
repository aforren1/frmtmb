#!/bin/sh
# Question 3 batch: one design per R process, one at a time.
# An entry "id:--laplace" runs that design with laplace = TRUE.
cd /c/Users/adf44/source/r/frmtmb-wt-samplers
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
for e in "$@"; do
  id="${e%%:*}"; fl=""; tag="$id"
  case "$e" in *:--laplace) fl="--laplace"; tag="$id-laplace";; esac
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/samplers-q3.R "$id" dev/samplers-q3-out $fl > "dev/samplers-q3-$tag.log" 2>&1
  echo "$tag rc=$?"
done
