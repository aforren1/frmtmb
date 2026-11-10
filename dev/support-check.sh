#!/usr/bin/env bash
# Lane support: R CMD build and R CMD check --as-cran of core, inside
# dev/support-check/, against the lane's library, then rellib-r19, then
# the user library (lane-env-r081.md, "Deliverables" item 5).
#
#   bash dev/support-check.sh
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-support
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-support-lib;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
rm -rf "$WT/dev/support-check"
mkdir -p "$WT/dev/support-check"
cd "$WT/dev/support-check"
"$R" CMD build "$WT"
tb=$(ls frmtmb_*.tar.gz)
"$R" CMD check --as-cran "$tb"
grep "^Status" frmtmb.Rcheck/00check.log
