#!/usr/bin/env bash
# R CMD check --as-cran of one package from this worktree, built and
# checked inside dev/cifix-check/<name>/.
#   bash dev/cifix-check.sh <name> <package directory>
set -eu
name="$1"; src="$2"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-cifix
RT=/c/rtools45
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
MW="$RT/x86_64-w64-mingw32.static.posix/bin"
export PATH="$PANDOC:$TEX:$RT/usr/bin:$MW:$PATH"
export RSTUDIO_PANDOC="C:${PANDOC#/c}"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
LIBS="C:/Users/adf44/source/r/wt-cifix-lib"
LIBS="$LIBS;C:/Users/adf44/source/r/rellib-r6"
export R_LIBS="$LIBS;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
: "${TMP:?TMP is unset}"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
out="$ROOT/dev/cifix-check/$name"
rm -rf "$out"; mkdir -p "$out"; cd "$out"
"$R" CMD build "$ROOT/$src" > build.log 2>&1
tgz=$(ls *.tar.gz)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1 || true
grep -a "^Status" check.log || tail -5 check.log
