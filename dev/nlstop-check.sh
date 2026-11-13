#!/usr/bin/env bash
# Lane nlstop: R CMD check --as-cran of core, built and checked in
# dev/nlstop-check/.
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop || exit 1
rm -rf dev/nlstop-check && mkdir -p dev/nlstop-check && cd dev/nlstop-check
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-nlstop-lib;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlstop/dev/stan-cache
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build ../.. > build.log 2>&1 || { tail -20 build.log; exit 1; }
tb=$(ls frmtmb_*.tar.gz)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|ERROR|WARNING|NOTE" check.log
