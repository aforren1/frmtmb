#!/bin/bash
# R CMD check --as-cran of one package of the lane worktree, built and
# checked inside dev/priortab-check/.
#   bash dev/priortab-check.sh .                     # core
#   bash dev/priortab-check.sh extensions/frmtmb.sample
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-priortab
PKG="$1"
R="/c/Program Files/R/R-4.6.1/bin"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
export PATH="$PANDOC:$TEX:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export R_LIBS="C:/Users/adf44/source/r/wt-priortab-lib;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
name=$(basename "$(cd "$WT/$PKG" && pwd)")
[ "$name" = frmtmb-wt-priortab ] && name=frmtmb
D="$WT/dev/priortab-check/$name"
rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1
"$R/R.exe" CMD build "$WT/$PKG" > build.log 2>&1 || {
  echo "build failed"; tail -5 build.log; exit 1; }
tb=$(ls ${name}_*.tar.gz | head -1)
"$R/R.exe" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|WARNING|ERROR|NOTE" check.log | head -20
