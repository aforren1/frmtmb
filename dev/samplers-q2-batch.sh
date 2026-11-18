#!/bin/sh
# Question 2 batch, then the asym_laplace design of question 3 through
# the formula route (dev/samplers-q3f.R). One R process at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-samplers
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
for id in p075 p124 p027 p113; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/samplers-q2.R "$id" > "dev/samplers-q2-$id.log" 2>&1
  echo "q2 $id rc=$?"
done
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/samplers-q3f.R ald dev/samplers-q3-out > dev/samplers-q3-ald.log 2>&1
echo "q3 ald rc=$?"
