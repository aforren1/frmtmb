#!/bin/sh
# R CMD build then R CMD check --as-cran, ONCE, on the final pass.
set -e
cd /c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad
rm -rf frailty-check
mkdir -p frailty-check
cd frailty-check

PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TINYTEX="/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="$PANDOC:$TINYTEX:$PATH"
export R_LIBS="C:/Users/adf44/source/r/frailtylib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE

R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build /c/Users/adf44/source/r/frmtmb-wt-frailty/extensions/frmtmb.spline
TARBALL=$(ls frmtmb.spline_*.tar.gz)
echo "built $TARBALL"
"$R" CMD check --as-cran --no-multiarch "$TARBALL"
echo "CHECK DONE"
