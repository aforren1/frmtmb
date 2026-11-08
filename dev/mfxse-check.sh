#!/usr/bin/env bash
# R CMD build and R CMD check --as-cran of core from the lane's worktree,
# inside dev/mfxse-check/, with pandoc and TinyTeX on PATH and R_LIBS the
# lane's library, then rellib-r18, then the user library.
#   bash dev/mfxse-check.sh
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-mfxse
R="/c/Program Files/R/R-4.6.1/bin"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="C:/Users/adf44/source/r/wt-mfxse-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
mkdir -p "$WT/dev/mfxse-check"
cd "$WT/dev/mfxse-check"
rm -rf frmtmb_*.tar.gz frmtmb.Rcheck
"$R/R.exe" CMD build "$WT"
"$R/R.exe" CMD check --as-cran frmtmb_*.tar.gz
