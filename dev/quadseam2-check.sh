#!/usr/bin/env bash
# R CMD check --as-cran of core, or of an extension (argument: its
# directory under extensions/), built and checked inside
# dev/quadseam2-check/ (lane-env r083, Deliverables item 5). Log
# dev/quadseam2-check/check-<pkg>.txt.
#   bash dev/quadseam2-check.sh [frmtmb.latent]
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
OUT=$WT/dev/quadseam2-check
PKG=${1:-}
SRC=$WT
[ -n "$PKG" ] && SRC=$WT/extensions/$PKG
NAME=${PKG:-frmtmb}
R="/c/Program Files/R/R-4.6.1/bin"
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-quadseam2-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=false
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
mkdir -p "$OUT/$NAME"; cd "$OUT/$NAME"
rm -rf ./*
{
  "$R/R.exe" CMD build "$SRC" || exit 1
  tb=$(ls ${NAME}_*.tar.gz)
  "$R/R.exe" CMD check --as-cran "$tb"
} > "$OUT/check-$NAME.txt" 2>&1
echo "exit $?"
grep -E "^Status|ERROR|WARNING|NOTE" "$OUT/check-$NAME.txt" | tail -12
