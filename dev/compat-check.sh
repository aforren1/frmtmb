#!/usr/bin/env bash
# R CMD check --as-cran once, on core, inside dev/compat-check/.
#   bash dev/compat-check.sh
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RCMD="/c/Program Files/R/R-4.6.1/bin/R.exe"
out="$ROOT/dev/compat-check"
mkdir -p "$out"
export PANDOC_DIR="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PANDOC_DIR:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-compat-lib;C:/Users/adf44/source/r/rellib-r19;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
cd "$out" || exit 1
"$RCMD" CMD build "$ROOT" > build.log 2>&1 || { tail -20 build.log; exit 1; }
tb=$(ls -t frmtmb_*.tar.gz | head -1)
"$RCMD" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|^\* DONE|ERROR|WARNING|NOTE" check.log | head -20
