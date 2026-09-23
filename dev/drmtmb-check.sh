#!/bin/sh
# R CMD build (with vignettes) and check --as-cran of the worktree's
# frmtmb, in the session scratchpad; the logs are copied to
# dev/drmtmb-log. drmtmb-lib comes first so the new Suggests entry is
# satisfied, pinlib before the user library for StanHeaders.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-drmtmb
S="$1"
export R_LIBS="C:/Users/adf44/source/r/drmtmb-lib;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS_USER="$R_LIBS"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
unset NOT_CRAN
R="C:/Program Files/R/R-4.6.1/bin/R.exe"
cd "$S" || exit 1
rm -rf frmtmb.Rcheck frmtmb_*.tar.gz
"$R" CMD build "$WT" > build.log 2>&1
echo "build exit $?"
TB=$(ls frmtmb_*.tar.gz)
"$R" CMD check --as-cran "$TB" > check.log 2>&1
echo "check exit $?"
cp build.log "$WT/dev/drmtmb-log/check-build.log"
cp check.log "$WT/dev/drmtmb-log/check-console.log"
cp frmtmb.Rcheck/00check.log "$WT/dev/drmtmb-log/check-00check.log"
