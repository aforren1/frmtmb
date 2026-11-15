#!/usr/bin/env bash
# Lane silent2: R CMD build and R CMD check --as-cran of one package
# inside dev/silent2-check/<package>/, against the lane library, then
# rellib-r21, then the user library.
#   bash dev/silent2-check.sh frmtmb|frmtmb.sample
set -u
: "${TMP:?TMP is unset}"
pkg="$1"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT"
[ "$pkg" = frmtmb ] || SRC="$ROOT/extensions/$pkg"
OUT="$ROOT/dev/silent2-check/$pkg"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
export PATH="$PANDOC:$TEX:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-silent2-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
"$R" CMD build "$SRC" > build.log 2>&1 || {
  tail -20 build.log; exit 1; }
tb=$(ls "${pkg}"_*.tar.gz)
"$R" CMD check --as-cran "$tb" > check.log 2>&1 || true
grep -E "^Status|ERROR|WARNING|NOTE|checking tests" check.log | head -40
