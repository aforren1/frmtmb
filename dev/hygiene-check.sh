#!/usr/bin/env bash
# R CMD check --as-cran of one package of this worktree, built and
# checked inside dev/hygiene-check/<pkg>/ (lane-env-r083.md, item 5).
#   bash dev/hygiene-check.sh core|frmtmb.ode|frmtmb.sample
set -u
p="$1"
ROOT="$(cd "$(dirname "$0")/.." && pwd -W)"
if [ "$p" = core ]; then src="$ROOT"; name=frmtmb; else
  src="$ROOT/extensions/$p"; name="$p"; fi
out="$ROOT/dev/hygiene-check/$name"
rm -rf "$out"; mkdir -p "$out"; cd "$out" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PANDOC:$TEX:$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-hygiene-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
: "${TMP:?TMP is unset}"
"$R" CMD build "$src" > build.log 2>&1 || {
  echo "BUILD FAILED"; tail -20 build.log; exit 1; }
tgz=$(ls "$name"_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1
grep -E "^Status|ERROR|WARNING|NOTE" check.log | tail -20
