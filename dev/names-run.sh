#!/usr/bin/env bash
# Lane names: run one R script from the worktree root with the round's
# toolchain (Rtools at the END of PATH, pandoc set). Usage:
#   bash dev/names-run.sh <script.R> <log> [args...]
cd /c/Users/adf44/source/r/frmtmb-wt-names || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-names/dev/stan-cache"
export NOT_CRAN=true
s="$1"; log="$2"; shift 2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$s" "$@" > "$log" 2>&1
echo "EXIT $?" >> "$log"
