#!/usr/bin/env bash
# R CMD build + check --as-cran of one package, inside dev/gpby-check/.
#   bash dev/gpby-check.sh core|frmtmb.sample|frmtmb.spline
set -u
pkg="$1"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-gpby
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
src="$ROOT"; [ "$pkg" = core ] || src="$ROOT/extensions/$pkg"
OUT="$ROOT/dev/gpby-check/$pkg"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-gpby-lib;C:/Users/adf44/source/r/rellib-r5;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=false
"$R" CMD build "$src" > build.log 2>&1
tgz=$(ls *.tar.gz | head -1)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1
grep -E "^Status|ERROR|WARNING|NOTE" check.log | tail -20
