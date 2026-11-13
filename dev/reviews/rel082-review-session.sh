#!/usr/bin/env bash
# Release review of 0.82.0: frmtmb.sample's suite in one session, reverse.
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel082
export PATH="/c/Program Files/R/R-4.6.1/bin:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_LIBS_USER= FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
Rscript.exe --vanilla dev/release/run-session.R frmtmb.sample reverse > dev/reviews/rel082-review-out/session-sample-reverse.log 2>&1
echo "exit $?" >> dev/reviews/rel082-review-out/session-sample-reverse.log
