#!/bin/sh
# Reviewer of lane loglik: the two bin-2 files on the lane build and on
# rellib-r13, one process each, sequentially (one staging directory).
cd /c/Users/adf44/source/r/frmtmb-wt-loglik || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export FRMTMB_STAN_CACHE=/c/Users/adf44/source/r/frmtmb-wt-loglik/dev/stan-cache
for arm in lane base; do
  if [ $arm = lane ]; then
    export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/wt-loglik-lib
    export R_LIBS_SITE=C:/Users/adf44/source/r/rellib-r13
  else
    export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r13
    unset R_LIBS_SITE
  fi
  for topic in log_lik posterior_predict; do
    f=dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-$topic.R
    "$R" dev/review-loglik-bin2-run.R frmtmb "$f" \
      "dev/review-loglik-log/rec-$arm-$topic.tsv" \
      > "dev/review-loglik-log/run-$arm-$topic.txt" 2>&1
    grep -h "^lib:\|^RESULT" "dev/review-loglik-log/run-$arm-$topic.txt"
  done
done
rm -rf dev/review-loglik-stage
