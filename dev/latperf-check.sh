#!/usr/bin/env bash
# Lane latperf: R CMD build and R CMD check --as-cran of core, inside
# dev/latperf-check/, against the lane library, rellib-r18 and the user
# library. Usage: bash dev/latperf-check.sh
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RB="/c/Program Files/R/R-4.6.1/bin"
CK="$ROOT/dev/latperf-check"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PANDOC:$TEX:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export R_LIBS="C:/Users/adf44/source/r/wt-latperf-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
rm -rf "$CK"
mkdir -p "$CK"
cd "$CK"
"$RB/R.exe" CMD build "$ROOT"
TB=$(ls frmtmb_*.tar.gz)
"$RB/R.exe" CMD check --as-cran "$TB"
