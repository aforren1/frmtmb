#!/usr/bin/env bash
# Reviewer: R CMD build and R CMD check --as-cran of frmtmb and
# frmtmb.spline from a copy of the worktree, with the OpenBLAS 0.3.26 R
# the lane built (dev/cifix-openblas.sh). Everything in
# dev/cifixrev-check/.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-cifix
CK="$ROOT/dev/cifixrev-check"
RT=/c/rtools45
PAN="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PAN:$TEX:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/cifixrev-lib;C:/Users/adf44/source/r/rellib-r6;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
export OPENBLAS_NUM_THREADS=4
: "${TMP:?TMP is unset}"
R="$ROOT/dev/cifix-out/Rob/bin/x64/R.exe"
pkg="$1"   # core or spline
if [ "$pkg" = core ]; then
  d="$CK/core"; srcd="$CK/src"
else
  d="$CK/spline"; srcd="$CK/src/extensions/frmtmb.spline"
fi
mkdir -p "$d"; cd "$d"
"$R" CMD build "$srcd" > build.log 2>&1
tb=$(ls -t *.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.log 2>&1
echo "DONE $pkg"
