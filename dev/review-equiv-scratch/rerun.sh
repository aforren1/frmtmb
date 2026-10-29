#!/usr/bin/env bash
# Reviewer: warm rerun of the lane's per-pair measurement on a subset,
# into a scratch output and a copied compile cache, so nothing of the
# lane's record is touched. equiv-run.sh is NOT called: its final step
# (equiv-summary.R) rewrites dev/equiv-data/results.tsv and the findings
# whatever outdir it is given.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-equiv
jobs=${1:-6}
out=dev/review-equiv-out2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$(cd dev/review-equiv-cache && pwd -W)"
export NOT_CRAN=true
mkdir -p "$out/logs"
ids="$(cat dev/review-equiv-scratch/rerun-ids.txt)"
start=$(date +%s)
printf '%s\n' $ids | xargs -P "$jobs" -I{} bash -c \
  "timeout 2400 \"$R\" dev/equiv-measure.R {} \"$out\" \
     > \"$out/logs/{}.log\" 2>&1 ||
   echo \"FAILED {} exit \$?\" >> \"$out/logs/{}.log\""
end=$(date +%s)
echo "measured in $((end - start)) s"
