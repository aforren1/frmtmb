#!/usr/bin/env bash
# Lane slopetrig: R CMD build and R CMD check --as-cran of core, inside
# dev/slopetrig-check/, against the lane library, then rellib-r22, then
# the user library.
# Usage: bash dev/slopetrig-check.sh
WT=/c/Users/adf44/source/r/frmtmb-wt-slopetrig
OUT="$WT/dev/slopetrig-check"
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-slopetrig-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT" || exit 1
"$R/R.exe" CMD build "$WT" > build.log 2>&1 || { echo "build failed"; exit 1; }
tb=$(ls frmtmb_*.tar.gz)
"$R/R.exe" CMD check --as-cran "$tb" > check.log 2>&1
grep -a "^Status" check.log
