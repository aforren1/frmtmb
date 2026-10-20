#!/usr/bin/env bash
# Lane prilabel: R CMD build and check --as-cran of one package, inside
# dev/prilabel-check/<pkg>/. Arg 1: core or an extension name.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-prilabel
p="${1:-core}"
if [ "$p" = core ]; then d="$WT"; n=frmtmb; else d="$WT/extensions/$p"; n="$p"; fi
OUT="$WT/dev/prilabel-check/$n"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
PD="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TT=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PD:$TT:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-prilabel-lib;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=false
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
unset FRMTMB_BRMS_FIT_TESTS FRMTMB_RELLIB
"$R" CMD build "$d" > build.txt 2>&1
tb=$(ls -t "$n"_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.txt 2>&1
grep -a "^Status" check.txt
