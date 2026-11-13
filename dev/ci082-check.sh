#!/usr/bin/env bash
# R CMD check --as-cran of core from this worktree, built and checked in
# dev/ci082-check/. Rtools goes at the END of PATH (lane-env-r083.md).
#   bash dev/ci082-check.sh
set -eu
ROOT=/c/Users/adf44/source/r/frmtmb-wt-ci082
RT=/c/rtools45
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PANDOC:$TEX:$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:${PANDOC#/c}"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
LIBS="C:/Users/adf44/source/r/wt-ci082-lib"
LIBS="$LIBS;C:/Users/adf44/source/r/rellib-r21"
export R_LIBS="$LIBS;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
: "${TMP:?TMP is unset}"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
out="$ROOT/dev/ci082-check"
rm -rf "$out"; mkdir -p "$out"; cd "$out"
"$R" CMD build "$ROOT" > build.log 2>&1
tgz=$(ls frmtmb_*.tar.gz)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1 || true
grep -a "^Status" check.log || tail -5 check.log
echo "CHECK DONE"
