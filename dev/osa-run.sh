#!/usr/bin/env bash
# Lane osa: run one R script with the round's toolchain on PATH.
#   bash dev/osa-run.sh <script.R> <log> [args...]
script="$1"; log="$2"; shift 2
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-osa/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
cd /c/Users/adf44/source/r/frmtmb-wt-osa || exit 1
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$script" "$@" > "$log" 2>&1
echo "exit $?" >> "$log"
