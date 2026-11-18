#!/usr/bin/env bash
# Lane sparsese: R CMD build and R CMD check --as-cran of core, inside
# dev/sparsese-check/, with the lane library first.
set -u
ROOT="C:/Users/adf44/source/r/frmtmb-wt-sparsese"
RB="/c/Program Files/R/R-4.6.1/bin"
RT=/c/rtools45
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-sparsese-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS_USER="C:/Users/adf44/source/r/wt-sparsese-lib"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
CK="$ROOT/dev/sparsese-check"
rm -rf "$CK"; mkdir -p "$CK"; cd "$CK"
"$RB/R.exe" CMD build "$ROOT" > build.txt 2>&1 || { tail -20 build.txt; exit 1; }
tb=$(ls frmtmb_*.tar.gz)
"$RB/R.exe" CMD check --as-cran "$tb" > check.txt 2>&1
grep -a "^Status" check.txt
