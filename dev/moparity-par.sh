#!/usr/bin/env bash
# Run test files one per process through dev/release/run-par.sh.
# Usage: bash dev/moparity-par.sh <tier> <lane|base> <pkg>:<file> ...
# lane puts wt-moparity-lib ahead of rellib-r7 (FRMTMB_RELLIB); base
# runs rellib-r7 alone. Output: dev/release/moparity-<tier>-files/.
tier=$1; arm=$2; shift 2
ROOT=/c/Users/adf44/source/r/frmtmb-wt-moparity
list=$ROOT/dev/release/moparity-$tier.jobs
: > "$list"
for j in "$@"; do
  p=${j%%:*}; f=${j#*:}
  if [ "$p" = "frmtmb" ]; then d=$ROOT/tests/testthat; else d=$ROOT/extensions/$p/tests/testthat; fi
  echo "$p $d/$f" >> "$list"
done
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
if [ "$arm" = "lane" ]; then export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-moparity-lib; else unset FRMTMB_RELLIB; fi
bash "$ROOT/dev/moparity-run-par.sh" "moparity-$tier" "$list" "${PAR:-16}"
grep -h "RESULT\|NO RESULT" "$ROOT/dev/release/moparity-$tier.log"
