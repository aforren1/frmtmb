#!/bin/sh
# R CMD check --as-cran, once per changed package, built and checked in
# dev/optima-rev-check/<package>/ (gitignored). Reviewer copy of
# dev/nanse-check.sh.
#
#   sh dev/optima-rev-check.sh frmtmb
#   sh dev/optima-rev-check.sh frmtmb.eam
#
# pandoc and TinyTeX on PATH (or a bogus pdflatex ERROR), no
# --no-manual, and the remote CRAN-incoming check off (frmtmb is not on
# CRAN, a pre-existing WARNING).
set -e
W=/c/Users/adf44/source/r/frmtmb-wt-optima
cd "$W"
PKG="$1"
case "$PKG" in
  frmtmb) SRC="$W" ;;
  *) SRC="$W/extensions/$PKG" ;;
esac
OUT="$W/dev/optima-rev-check/$PKG"
rm -rf "$OUT"
mkdir -p "$OUT"
export R_LIBS="C:/Users/adf44/source/r/wt-optima-lib;C:/Users/adf44/source/r/rellib-r6;C:/Users/adf44/AppData/Local/R/win-library/4.6"
pd="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
tex="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$pd:$tex:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/optima-rev-out/stan-cache"
: "${TMP:?TMP is unset}"
cd "$OUT"
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD build "$SRC" > build.log 2>&1
TAR=$(sed -n "s/^.*building '\(.*\)'.*$/\1/p" build.log)
if [ -z "$TAR" ]; then echo "BUILD FAILED"; cat build.log; exit 1; fi
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD check --as-cran \
  --output="$OUT" "$TAR" > check.log 2>&1 || true
grep -h "^Status" "$OUT"/*.Rcheck/00check.log
