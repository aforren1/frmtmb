#!/usr/bin/env bash
# Lane drawdir: R CMD build and R CMD check --as-cran of core, inside
# dev/drawdir-check/, with the lane library, rellib-r22 and the user
# library, pandoc and TinyTeX on PATH.
set -u
unset NOT_CRAN FRMTMB_BRMS_FIT_TESTS FRMTMB_DRMTMB_FIT_TESTS
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drawdir
C="$ROOT/dev/drawdir-check"
rm -rf "$C"; mkdir -p "$C"; cd "$C"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-drawdir-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE R_PARALLEL_PORT=11849
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$ROOT" > build.txt 2>&1 || { echo build failed; tail build.txt; exit 1; }
tb=$(ls frmtmb_*.tar.gz)
"$R" CMD check --as-cran "$tb" > check.txt 2>&1
grep -E "^Status|WARNING|NOTE|ERROR" check.txt
echo check done
