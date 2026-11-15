#!/usr/bin/env bash
# Lane slopescale: R CMD build and R CMD check --as-cran of core, inside
# dev/slopescale-check/, against the lane library, then rellib-r21, then
# the user library.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-slopescale
D="$WT/dev/slopescale-check"
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-slopescale-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export R_PARALLEL_PORT=11875
export NOT_CRAN=true
rm -rf "$D"; mkdir -p "$D"; cd "$D"
"$R/R.exe" CMD build "$WT" > build.log 2>&1 || {
  echo "BUILD FAILED"; tail -20 build.log; exit 1; }
tb=$(ls frmtmb_*.tar.gz)
"$R/R.exe" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|WARNING|ERROR|NOTE" check.log
