#!/usr/bin/env bash
# Lane nlband: the lane's test file and the nonlinear files, with the
# helpers and setup.R, in one session each way (dev/release/run-session.R),
# one session at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
for o in forward reverse preload; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/release/run-session.R \
    frmtmb "$o" dev/nlband-session-tests > "dev/nlband-out/session-$o.txt" 2>&1
done
echo "SESSIONS DONE"
