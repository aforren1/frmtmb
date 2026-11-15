#!/usr/bin/env bash
# Review of lane hygiene, item 3631: every core test file that reads the
# joint-density translator of helper-brms.R, gated, on the final tree.
#   bash dev/reviews/hygiene-review-gated.sh > dev/reviews/hygiene-review-gated.log
cd /c/Users/adf44/source/r/frmtmb-wt-hygiene || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=/c/Users/adf44/source/r/frmtmb-wt-hygiene/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
grep -l -E "brms_lp_check|stan_pars_from_fit|brms_group_pars" \
  tests/testthat/test-*.R | awk '{print NR, $0}' | xargs -P 3 -L 1 bash -c '
  port=$((12210 + $0))
  t0=$(date +%s)
  out=$(HYG_ARM=mine R_PARALLEL_PORT=$port "'"$RS"'" \
        dev/hygiene-run-tests.R frmtmb $1 2>&1)
  echo "$out" | grep -E "^(RESULT|DETAIL|lib:)"
  echo "TIME $(basename $1) $(( $(date +%s) - t0 )) s"
'
