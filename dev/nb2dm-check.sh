#!/bin/sh
# R CMD build and check --as-cran of core and frmtmb.sample, inside
# dev/nb2dm-check/, against the lane library, then rellib-r10, then the
# user library.   sh dev/nb2dm-check.sh [frmtmb|frmtmb.sample]
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
OUT="$ROOT/dev/nb2dm-check"
mkdir -p "$OUT"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-nb2dm-lib;C:/Users/adf44/source/r/rellib-r10;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
for p in "$@"; do
  if [ "$p" = frmtmb ]; then src="$ROOT"; else src="$ROOT/extensions/$p"; fi
  d="$OUT/$p"
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  "$RC" CMD build "$src" > build.txt 2>&1
  tgz=$(ls "$p"_*.tar.gz)
  "$RC" CMD check --as-cran "$tgz" > check.txt 2>&1 || true
  grep -h "^Status" check.txt "$p.Rcheck/00check.log" 2>/dev/null | tail -1
  cd "$ROOT"
done
