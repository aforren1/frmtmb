#!/usr/bin/env bash
# R CMD build and R CMD check --as-cran of core, inside dev/fitedge-check/,
# with pandoc and TinyTeX on PATH and R_LIBS = the lane's library, then
# rellib-r18, then the user library.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-fitedge
OUT=$WT/dev/fitedge-check
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PANDOC:$TEX:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-fitedge-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/stan-cache
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
stray=$(cd "$WT" && find . -path ./dev -prune -o \( -name "testthat-problems.rds" \
  -o -name "_problems" -o -name "Rplots.pdf" -o -name "*.Rcheck" \) -print)
if [ -n "$stray" ]; then echo "STRAY FILES, not building:"; echo "$stray"; exit 1; fi
"$R" CMD build "$WT" > build.log 2>&1
tb=$(ls frmtmb_*.tar.gz | head -1)
"$R" CMD check --as-cran "$OUT/$tb" > check.log 2>&1
grep -E "^Status|ERROR|WARNING|NOTE" check.log | head -30
