#!/usr/bin/env bash
# Lane postfix: R CMD build and R CMD check --as-cran of one package of
# this worktree, inside dev/postfix-check/<pkg>/, against the lane
# library, then rellib-r8, then the user library.
#   bash dev/postfix-check.sh <package dir relative to the worktree>
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
src="$ROOT/$1"
name="$(basename "$src")"
[ "$1" = "." ] && name=frmtmb
out="$ROOT/dev/postfix-check/$name"
rm -rf "$out"; mkdir -p "$out"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
L1=C:/Users/adf44/source/r/wt-postfix-lib
L2=C:/Users/adf44/source/r/rellib-r8
L3=C:/Users/adf44/AppData/Local/R/win-library/4.6
export R_LIBS="$L1;$L2;$L3"
export NOT_CRAN=true _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
PAN="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="$PAN"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
RT=/c/rtools45
MG=$RT/x86_64-w64-mingw32.static.posix/bin
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$TEX:$RT/usr/bin:$MG:$PATH"
cd "$out"
"$R" CMD build "$src" > build.log 2>&1
tgz="$(ls -t *.tar.gz 2>/dev/null | head -1)"
[ -n "$tgz" ] || { echo "NO TARBALL"; tail -20 build.log; exit 1; }
"$R" CMD check --as-cran "$tgz" > check.log 2>&1
grep -a "^Status:" check.log || tail -20 check.log
