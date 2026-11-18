#!/bin/sh
# Question 1 large designs without the mock backend, one at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-samplers
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
for id in "$@"; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/samplers-q1raw.R "$id" > "dev/samplers-q1-$id.log" 2>&1
  echo "$id rc=$?"
done
