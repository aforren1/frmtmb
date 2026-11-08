#!/usr/bin/env bash
# Release review 0.80.0: R CMD build and check --as-cran of core from the
# release tree, outside it (rev-rel-src/check), on rellib-r19.
#   bash dev/review-rel-check.sh
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel080
D=/c/Users/adf44/source/r/rev-rel-src/check
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
cd "$D" || exit 1
echo "build $(date +%T)"
"$R" CMD build "$ROOT" > build.log 2>&1
tgz=$(ls -t frmtmb_*.tar.gz | head -1)
echo "check $tgz $(date +%T)"
"$R" CMD check --as-cran "$D/$tgz" > check.log 2>&1
grep -E "^Status:" check.log
echo "CHECK END $(date +%T)"
