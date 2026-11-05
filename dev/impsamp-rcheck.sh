#!/usr/bin/env bash
# R CMD check --as-cran of one package, built and checked inside
# dev/impsamp-check/:  bash dev/impsamp-rcheck.sh <package dir> <label>
set -u
src="$1"
label="$2"
ROOT="$(pwd)"
RB="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export R_LIBS="C:/Users/adf44/source/r/wt-impsamp-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
out="$ROOT/dev/impsamp-check/$label"
mkdir -p "$out"
cd "$out"
"$RB/R.exe" CMD build "$ROOT/$src" > build.txt 2>&1
tb=$(ls -t *.tar.gz | head -1)
"$RB/R.exe" CMD check --as-cran "$tb" > check.txt 2>&1
tail -3 check.txt
echo "CHECK-DONE $label"
