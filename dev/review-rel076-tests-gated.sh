#!/usr/bin/env bash
# Review: the same files with every gate on (no session), 4 at a time.
# Review: the brief's test files on rellib-r14, one file per process,
# 6 at a time, then one core session in reverse order. Rtools is put
# LAST on PATH (the brief: Rtools first gave the wrong library before).
#   bash dev/review-rel076-tests.sh
REL=/c/Users/adf44/source/r/frmtmb-wt-rel076
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$REL/dev/review-rel076-out/tests-gated
export FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$REL/dev/review-rel076-out/stan-cache
export ROUND_BASE=C:/Users/adf44/source/r/rellib-r14
T=tests/testthat
S=extensions/frmtmb.sample/tests/testthat
P=extensions/frmtmb.spatial/tests/testthat
jobs_list=(
"frmtmb $T/test-mo-coding.R" "frmtmb $T/test-autocor-loglik.R"
"frmtmb $T/test-autocor-seam.R" "frmtmb $T/test-priors-autocor-classes.R"
"frmtmb $T/test-emmeans.R" "frmtmb $T/test-se-check.R"
"frmtmb $T/test-compat.R" "frmtmb $T/test-message-uniqueness.R"
"frmtmb $T/test-guardrail1.R" "frmtmb $T/test-extension-home-mixed.R"
"frmtmb.sample $S/test-mo-coding.R"
"frmtmb.sample $S/test-default-priors-autocor.R"
"frmtmb.sample $S/test-default-priors-bounds.R"
"frmtmb.sample $S/test-default-priors-fit.R"
"frmtmb.sample $S/test-loglik-autocor.R" "frmtmb.sample $S/test-loo.R"
"frmtmb.spatial $P/test-sar.R" "frmtmb.spatial $P/test-guardrail.R"
"BIN2 dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R"
"BIN2 dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-posterior_predict.R"
)
cd "$REL" || exit 1
for j in "${jobs_list[@]}"; do
  set -- $j
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do sleep 2; done
  o="$OUT/$1--$(basename "$2" .R).log"
  if [ "$1" = BIN2 ]; then
    ( "$RS" "$REL/dev/review-rel076-bin2run.R" frmtmb "$REL/$2" > "$o" 2>&1 ) &
  else
    ( "$RS" "$REL/dev/ciharden-run1.R" base "$1" "$REL/$2" > "$o" 2>&1 ) &
  fi
done
wait
for f in "$OUT"/*--*.log; do
  echo "$(basename "$f") | $(grep -a '^RESULT' "$f" | tail -1) | $(grep -a '^lib:' "$f" | head -1)"
done > "$OUT/summary.txt"
echo ALLDONE >> "$OUT/summary.txt"
