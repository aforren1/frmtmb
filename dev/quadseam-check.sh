#!/usr/bin/env bash
# R CMD check --as-cran of core, built and checked inside
# dev/quadseam-check/ (lane-env r082, Deliverables item 5). Log
# dev/quadseam-check/check.txt.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam
OUT=$WT/dev/quadseam-check
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-quadseam-lib;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=false
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
{
  "$R/R.exe" CMD build "$WT" || exit 1
  tb=$(ls frmtmb_*.tar.gz)
  "$R/R.exe" CMD check --as-cran "$tb"
} > "$OUT/check.txt" 2>&1
echo "exit $?"
grep -E "^Status|ERROR|WARNING|NOTE" "$OUT/check.txt" | tail -12
