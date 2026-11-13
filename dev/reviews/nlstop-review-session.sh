#!/bin/sh
# One core session over the review's subset (dev/reviews/nlstop-review-session):
# nlstop-review-session.sh <order>
cd C:/Users/adf44/source/r/frmtmb-wt-nlstop
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_PARALLEL_PORT=$((12000 + $$ % 1000))
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop-lib
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/release/run-session.R frmtmb $1 dev/reviews/nlstop-review-session > dev/reviews/nlstop-review-out/session-$1.txt 2>&1
echo "EXIT $?" >> dev/reviews/nlstop-review-out/session-$1.txt
