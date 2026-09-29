#!/usr/bin/env bash
# Lane wt-arcovsample: R CMD check --as-cran, ONCE per changed package.
# Both packages changed, so both run. Not --no-manual: the manual
# sections are where an Rd defect surfaces (lane-rules.md).
#
#   bash dev/arcovsample-check.sh <core|sample>
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-arcovsample
LIB="C:/Users/adf44/source/r/wt-arcovsample-lib"
USER="C:/Users/adf44/AppData/Local/R/win-library/4.6"
RBIN="/c/Program Files/R/R-4.6.1/bin"
# The lane library holds frmtmb and frmtmb.sample alone
# (dependencies = FALSE), and frmtmb.sample SUGGESTS three siblings, so
# --as-cran stops at "Packages suggested but not available". The
# reference build supplies them, READ-ONLY and AFTER the lane library,
# so the frmtmb and frmtmb.sample under test are this lane's.
# drmTMB, which core suggests, is in neither: see the findings.
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r3;$USER"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"

case "${1:-}" in
  core)   PKG="$W" ; OUT="$W/dev/arcovsample-log/check-core" ;;
  sample) PKG="$W/extensions/frmtmb.sample"
          OUT="$W/dev/arcovsample-log/check-sample" ;;
  *) echo "usage: arcovsample-check.sh <core|sample>"; exit 2 ;;
esac
mkdir -p "$OUT"
cd "$OUT"
# built WITH vignettes: a tarball with vignette sources and no inst/doc
# manufactures two WARNINGs and a NOTE (dev/release/run-check.ps1)
"$RBIN/R.exe" CMD build "$PKG"
TAR=$(ls -t ./*.tar.gz | head -1)
echo "built $TAR"
"$RBIN/R.exe" CMD check --as-cran "$TAR"
