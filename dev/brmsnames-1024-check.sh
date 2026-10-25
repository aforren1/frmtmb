#!/usr/bin/env bash
# The lane's one R CMD check --as-cran per changed package, without
# --no-manual, built and checked inside dev/brmsnames-1024-check/<pkg>.
# pandoc and TinyTeX go on PATH or the check reports a bogus pdflatex
# ERROR (dev/lane-rules.md).
#
#   bash dev/brmsnames-1024-check.sh frmtmb|frmtmb.sample
set -u
pkg="$1"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [ "$pkg" = frmtmb ]; then SRC="$ROOT"; else SRC="$ROOT/extensions/$pkg"; fi
OUT="$ROOT/dev/brmsnames-1024-check/$pkg"
mkdir -p "$OUT"; cd "$OUT" || exit 1
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="C:/Users/adf44/source/r/wt-brmsnames-lib;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$SRC" > build.log 2>&1
tar=$(ls -t "$pkg"_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tar" > check.log 2>&1
grep -h "^Status" check.log
