#!/usr/bin/env bash
# Lane boundfam: R CMD build and check --as-cran of one package, inside
# dev/boundfam-check/<pkg>/.
#   bash dev/boundfam-check.sh <package dir relative to the worktree>
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT/$1"
P=$(basename "$SRC"); [ "$1" = "." ] && P=frmtmb
OUT="$ROOT/dev/boundfam-check/$P"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
RT=/c/rtools45
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="$PANDOC"
export R_LIBS="C:/Users/adf44/source/r/wt-boundfam-lib;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$SRC" > build.log 2>&1 || { echo "build failed"; tail build.log; exit 1; }
TB=$(ls *.tar.gz | head -1)
"$R" CMD check --as-cran "$TB" > check.log 2>&1
grep "^Status" check.log || tail -5 check.log
