#!/usr/bin/env bash
# Lane fixes: R CMD build and check --as-cran of one package, inside
# dev/fixes-check/.   bash dev/fixes-check.sh <package dir> <name>
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PKG="$1"; NAME="$2"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
PD="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TT="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
RT=/c/rtools45
export PATH="$PD:$TT:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="$PD"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-fixes-lib;C:/Users/adf44/source/r/rellib-r5;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
: "${TMP:?TMP is unset}"
W="$ROOT/dev/fixes-check/$NAME"
rm -rf "$W"; mkdir -p "$W"; cd "$W"
"$R" CMD build "$ROOT/$PKG" > build.log 2>&1
tb=$(ls *.tar.gz | head -n 1)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
grep -a "^Status" check.log
