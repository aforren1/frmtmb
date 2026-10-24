#!/usr/bin/env bash
# R CMD build and R CMD check --as-cran for one package of the lane,
# inside dev/splitsurv-check/<pkg>/ (lane-env.md, Deliverables 5).
#   bash dev/splitsurv-check.sh <frmtmb | extension name>
set -u
p="$1"
ROOT=C:/Users/adf44/source/r/frmtmb-wt-splitsurv
LIB=C:/Users/adf44/source/r/wt-splitsurv-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r11;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
if [ "$p" = frmtmb ]; then src="$ROOT"; else src="$ROOT/extensions/$p"; fi
out="$ROOT/dev/splitsurv-check/$p"
rm -rf "$out"; mkdir -p "$out"; cd "$out"
"$R" CMD build "$src" > build.log 2>&1 || {
  echo "BUILD FAILED"; tail -20 build.log; exit 1; }
tb=$(ls ${p}_*.tar.gz)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|WARNING|ERROR|NOTE" check.log | head -40
