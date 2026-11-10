#!/usr/bin/env bash
# Lane support: R CMD build and R CMD check --as-cran of one extension,
# inside dev/support-<package>-check/, against the lane's library, then
# rellib-r19, then the user library (lane-env-r081.md, "Deliverables"
# item 5). The check log is copied to dev/support-out/00check-<package>.log.
# Core's check is dev/support-check.sh.
#
#   bash dev/support-check-pkg.sh <package>
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-support
PKG="$1"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-support-lib;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
D="$WT/dev/support-$PKG-check"
rm -rf "$D"
mkdir -p "$D"
cd "$D"
"$R" CMD build "$WT/extensions/$PKG"
tb=$(ls "$PKG"_*.tar.gz)
"$R" CMD check --as-cran "$tb" || true
cp "$PKG.Rcheck/00check.log" "$WT/dev/support-out/00check-$PKG.log"
grep "^Status" "$PKG.Rcheck/00check.log"
