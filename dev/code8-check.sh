#!/usr/bin/env bash
# Lane code8: R CMD check --as-cran of core, built and checked inside
# dev/code8-check/.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/dev/code8-check"
mkdir -p "$OUT"; cd "$OUT"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-code8-lib;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
rm -f frmtmb_*.tar.gz
"$R" CMD build "$ROOT" > build.log 2>&1
tb=$(ls frmtmb_*.tar.gz)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
grep -a "Status:" check.log
