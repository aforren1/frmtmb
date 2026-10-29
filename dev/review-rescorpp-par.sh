#!/usr/bin/env bash
# Review rescorpp: one test file per R process, N at a time.
#   bash dev/review-rescorpp-par.sh <lane|base> <jobs> <outdir> [N]
set -u
arm="$1"; jobs="$2"; out="$3"; N="${4:-5}"
WT=/c/Users/adf44/source/r/frmtmb-wt-rescorpp
mkdir -p "$out"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
while read -r pkg f; do
  [ -z "$pkg" ] && continue
  echo "$pkg $f"
done < "$jobs" | xargs -P "$N" -L 1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" '"$WT"'/dev/rescorpp-runtest.R '"$arm"' "$0" "$1" > '"$out"'/"$0--$1.txt" 2>&1'
grep -h "^RESULT\|^lib:" "$out"/*.txt
