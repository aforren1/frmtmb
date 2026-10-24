#!/usr/bin/env bash
# Reviewer's driver: bash dev/review-mvnlatent-par.sh <jobs> <pkg:file>...
# Each file in its own Rscript, at most <jobs> at once; logs go to
# dev/review-mvnlatent-out/<pkg>__<file>.txt
WT=/c/Users/adf44/source/r/frmtmb-wt-mvnlatent
OUT=$WT/dev/review-mvnlatent-out
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$OUT"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-mvnlatent/dev/stan-cache
J=$1; shift
for spec in "$@"; do
  pkg=${spec%%:*}; f=${spec#*:}
  if [ "$pkg" = frmtmb ]; then d=$WT/tests/testthat; else d=$WT/extensions/$pkg/tests/testthat; fi
  while [ "$(jobs -rp | wc -l)" -ge "$J" ]; do wait -n; done
  ( cd "$d" && "$RS" "$WT/dev/review-mvnlatent-run1.R" "$pkg" "$f" \
      > "$OUT/${pkg}__${f%.R}.txt" 2>&1 ) &
done
wait
grep -h "^RESULT" "$OUT"/*.txt
