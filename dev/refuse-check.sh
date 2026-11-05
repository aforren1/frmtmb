#!/usr/bin/env bash
# R CMD build and check --as-cran, inside dev/refuse-check/.
#   bash dev/refuse-check.sh <subdir of the worktree, "." for core>
set -u
ROOT="$(pwd)"
SRC="$1"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export R_LIBS="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
mkdir -p dev/refuse-check
cd dev/refuse-check
"$R" CMD build --no-manual "$ROOT/$SRC" > "build-$(basename "$(cd "$ROOT/$SRC" && pwd)").log" 2>&1
pkg=$(basename "$(cd "$ROOT/$SRC" && pwd)")
[ "$SRC" = "." ] && pkg=frmtmb
tb=$(ls -t ${pkg}_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > "check-$pkg.log" 2>&1
grep -a '^Status' "check-$pkg.log" || tail -3 "check-$pkg.log"
