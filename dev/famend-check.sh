#!/usr/bin/env bash
# Lane famend: R CMD build and R CMD check --as-cran of core, inside
# dev/famend-check/, with pandoc and TinyTeX on PATH and the lane's
# library first.
#   bash dev/famend-check.sh
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-famend
cd "$WT/dev/famend-check"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-famend-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=false
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$WT" > build.txt 2>&1
TB=$(ls -t frmtmb_*.tar.gz | head -1)
"$R" CMD check --as-cran "$TB" > check.txt 2>&1 || true
grep -a "^Status" check.txt || tail -20 check.txt
