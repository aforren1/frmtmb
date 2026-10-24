#!/bin/bash
# Usage: review-spatial-par.sh <jobs file> <libs> <nproc> <logdir>
# Each job line: <pkg> <test file> <log name>
jobs=$1; libs=$2; np=$3; out=$4
mkdir -p "$out"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-spatial/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
cat "$jobs" | xargs -P "$np" -L 1 bash -c '"$0" C:/Users/adf44/source/r/frmtmb-wt-spatial/dev/review-spatial-run.R "'"$libs"'" "$1" "$2" > "'"$out"'/$3.log" 2>&1' "$R"
grep -h "^RESULT\|^LOADERROR" "$out"/*.log
