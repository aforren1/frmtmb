#!/usr/bin/env bash
# R CMD build (with vignettes) and R CMD check --as-cran of one package
# of the worktree, inside dev/mvnlatent-check/<name>/, against the
# lane's library, then rellib-r11, then the user library.
#
#   bash dev/mvnlatent-check.sh <package dir relative to the worktree>
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
src="$ROOT/$1"
name=$(basename "$src"); [ "$1" = "." ] && name=frmtmb
OUT="$ROOT/dev/mvnlatent-check/$name"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
PD="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TT=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PD:$TT:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-mvnlatent-lib;C:/Users/adf44/source/r/rellib-r11;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-C:/Users/adf44/source/r/frmtmb-wt-mvnlatent/dev/stan-cache}"
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
"$R" CMD build "$src" > build.txt 2>&1
tgz=$(ls -t *.tar.gz 2>/dev/null | head -1)
if [ -z "$tgz" ]; then echo "NO TARBALL for $name"; tail -20 build.txt; exit 1; fi
"$R" CMD check --as-cran "$tgz" > check.txt 2>&1
grep -a "^Status" check.txt | sed "s/^/$name: /"
