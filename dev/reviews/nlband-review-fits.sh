#!/usr/bin/env bash
# Review of lane nlband: dev/reviews/nlband-review-fits.R on both arms, one
# process per arm, run side by side (2 R processes).
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=dev/reviews/nlband-review-out
SECS="${1:-seeds,edge,inside,dense,two,visits,jitter,student,poisson,bstep,bband}"
TAG="${2:-fits}"
NLBAND_ARM=before "$RS" dev/reviews/nlband-review-fits.R "$SECS" \
  > $OUT/$TAG-before.txt 2>&1 &
NLBAND_ARM=after "$RS" dev/reviews/nlband-review-fits.R "$SECS" \
  > $OUT/$TAG-after.txt 2>&1 &
wait
echo "FITS DONE"
