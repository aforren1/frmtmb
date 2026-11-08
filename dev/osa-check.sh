#!/usr/bin/env bash
# Lane osa: R CMD build and R CMD check --as-cran of core, inside
# dev/osa-check/, against wt-osa-lib, then rellib-r18, then the user
# library.
#   bash dev/osa-check.sh
set -u
WT="C:/Users/adf44/source/r/frmtmb-wt-osa"
OUT="$WT/dev/osa-check"
R="/c/Program Files/R/R-4.6.1/bin"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="C:/Users/adf44/source/r/wt-osa-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT" || exit 1
"$R/R.exe" CMD build "$WT" > build.log 2>&1 || { tail -20 build.log; exit 1; }
tb=$(ls frmtmb_*.tar.gz | head -1)
"$R/R.exe" CMD check --as-cran "$tb" > check.log 2>&1
grep -n "^Status" check.log
