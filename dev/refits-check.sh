#!/usr/bin/env bash
# R CMD build and check --as-cran of one package of the worktree, inside
# dev/refits-check/. Usage: bash dev/refits-check.sh <pkg dir rel. to WT>
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-refits
PKG="$1"
OUT="$WT/dev/refits-check/$(basename "$PKG")"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="C:/Users/adf44/source/r/wt-refits-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=false
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
"$R" CMD build "$WT/$PKG" > build.log 2>&1
tb=$(ls *.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|WARNING|ERROR|NOTE" check.log
