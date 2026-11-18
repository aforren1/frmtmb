#!/usr/bin/env bash
# Lane quadesc: R CMD build and check --as-cran of core, inside
# dev/quadesc-check/, with the lane library first.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-quadesc
cd "$WT"
mkdir -p dev/quadesc-check
cd dev/quadesc-check
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-quadesc-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$WT" > build.txt 2>&1
TB=$(ls -t frmtmb_*.tar.gz | head -1)
"$R" CMD check --as-cran "$TB" > check.txt 2>&1
grep "^Status" check.txt
