#!/bin/bash
# R CMD check --as-cran on frmtmb.ode, once, on the final pass.
#
# pandoc and TinyTeX go on PATH or the check reports a bogus pdflatex
# ERROR. The CRAN-incoming remote check is off because "not on CRAN" is
# pre-existing here and is not this lane's to answer.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-check.sh
set -u
SRC=C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode
OUT=C:/Users/adf44/source/r/nss-check
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="C:/Users/adf44/source/r/nsslib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT" || exit 1
"$R/R.exe" CMD build "$SRC" > "$OUT/build.log" 2>&1
TAR=$(ls -1 frmtmb.ode_*.tar.gz | head -1)
echo "built $TAR"
"$R/R.exe" CMD check --as-cran "$TAR" > "$OUT/check.log" 2>&1
tail -40 "$OUT/check.log"
