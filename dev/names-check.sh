#!/usr/bin/env bash
# Lane names: R CMD build and R CMD check --as-cran of one package, inside
# dev/names-check/. Usage: bash dev/names-check.sh <core|frmtmb.sample>
cd /c/Users/adf44/source/r/frmtmb-wt-names || exit 1
p="$1"
if [ "$p" = core ]; then src=/c/Users/adf44/source/r/frmtmb-wt-names; else
  src=/c/Users/adf44/source/r/frmtmb-wt-names/extensions/$p; fi
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-names-lib;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-names/dev/stan-cache"
out=/c/Users/adf44/source/r/frmtmb-wt-names/dev/names-check/$p
rm -rf "$out"; mkdir -p "$out"; cd "$out" || exit 1
"$R/R.exe" CMD build "$src" > build.log 2>&1 || {
  echo "BUILD FAIL"; tail -20 build.log; exit 1; }
tb=$(ls *.tar.gz | head -1)
"$R/R.exe" CMD check --as-cran "$tb" > check.log 2>&1
grep -a "^Status" check.log
grep -a -B1 -A6 "WARNING\|ERROR\|NOTE" check.log | head -80
