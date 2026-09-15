#!/usr/bin/env bash
# R CMD check --as-cran, ONCE, on the final pass.  Mirrors
# dev/release/run-check.ps1 so the result is comparable to the release
# harness's own.
#
#   bash dev/generics-cran.sh <package dir> <outdir>
#
# Built WITH vignettes: a tarball with vignette sources and no inst/doc
# manufactures two WARNINGs and a NOTE that are not findings.  No
# --no-manual either: it skips the HTML-manual section, which is where
# an Rd defect surfaces and where this machine's V8 NOTE comes from.
set -u
PKG="$1"
OUT="$2"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export R_LIBS="C:/Users/adf44/source/r/generics-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-generics/dev/stan-cache"
rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT" || exit 1
"$R" CMD build "$PKG" || exit 1
TAR=$(ls -t *.tar.gz | head -1)
echo "BUILT $TAR"
"$R" CMD check --as-cran "$TAR"
echo "CHECK EXIT $?"
