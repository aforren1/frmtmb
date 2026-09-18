#!/bin/sh
# R CMD check --as-cran, ONCE per package at the end of the lane.
#
#   sh dev/adefects-check.sh frmtmb
#   sh dev/adefects-check.sh frmtmb.sample
#
# pandoc and TinyTeX must be on PATH or the check reports a bogus
# pdflatex ERROR, and --no-manual is NOT passed: it skips the HTML and
# PDF manual sections, which is where an unescaped `%` in Rd surfaces.
# _R_CHECK_CRAN_INCOMING_REMOTE_ is off because the CRAN-incoming
# WARNING about frmtmb not being on CRAN is pre-existing and remote.
set -e
cd /c/Users/adf44/source/r/frmtmb-wt-adefects
PKG="$1"
if [ -z "$PKG" ]; then echo "usage: adefects-check.sh <package>"; exit 1; fi
case "$PKG" in
  frmtmb) SRC=. ;;
  *) SRC="extensions/$PKG" ;;
esac
OUT="/c/Users/adf44/source/r/adefects-check"
mkdir -p "$OUT"
lane="C:/Users/adf44/source/r/adefects-lib"
ref="C:/Users/adf44/source/r/rellib-r3"
usr="C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS="$lane;$ref;$usr"
pd="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
tex="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
rsp="C:/Program Files/RStudio/resources/app/bin/quarto/bin"
export RSTUDIO_PANDOC="$rsp/tools"
export PATH="$pd:$tex:/c/rtools45/usr/bin:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
# build and check from OUT, so the tarball does not land in the tree
SRC=$(cd "$SRC" && pwd)
cd "$OUT"
TAR=$("/c/Program Files/R/R-4.6.1/bin/R.exe" CMD build "$SRC" 2>&1 |
      tee /dev/stderr | sed -n "s/^.*building '\\(.*\\)'.*$/\\1/p")
if [ -z "$TAR" ]; then echo "BUILD FAILED"; exit 1; fi
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD check --as-cran \
  --output="$OUT" "$TAR"
echo "=== 00check.log ==="
cat "$OUT/$PKG.Rcheck/00check.log"
