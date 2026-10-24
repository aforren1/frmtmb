#!/usr/bin/env bash
# R CMD check --as-cran of one package of this worktree, built and
# checked inside dev/spatial-check/<pkg>/, against the lane library.
#   bash dev/spatial-check.sh <core|frmtmb.x>
set -u
p="$1"
ROOT="$(cd "$(dirname "$0")/.." && pwd -W)"
LIB=C:/Users/adf44/source/r/wt-spatial-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r11;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
if [ "$p" = core ]; then src="$ROOT"; name=frmtmb; else src="$ROOT/extensions/$p"; name=$p; fi
W="$ROOT/dev/spatial-check/$name"
rm -rf "$W"; mkdir -p "$W"; cd "$W"
"$R" CMD build "$src" > build.log 2>&1 || { echo "BUILD FAILED $name"; tail -20 build.log; exit 1; }
tb=$(ls "$name"_*.tar.gz)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|ERROR|WARNING|NOTE" check.log | head -20
