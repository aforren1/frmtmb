#!/usr/bin/env bash
# R CMD build and R CMD check --as-cran of one package of this worktree,
# inside dev/readers-check/<pkg>, against the lane library, then rellib-r18,
# then the user library.
#   bash dev/readers-check.sh <pkg dir relative to the worktree> <name>
set -eu
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/$1"; NAME="$2"
OUT="$ROOT/dev/readers-check/$NAME"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="C:/Users/adf44/source/r/wt-readers-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$SRC" > build.txt 2>&1
tgz=$(ls "$NAME"_*.tar.gz)
"$R" CMD check --as-cran "$tgz" > check.txt 2>&1 || true
grep -a "^Status" check.txt || tail -20 check.txt
