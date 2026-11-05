#!/usr/bin/env bash
# The lane's final runs, in order. Logs under dev/nlfix-*.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
bash dev/nlfix-install.sh > dev/nlfix-install-final.log 2>&1
rm -f dev/nlfix-suite-ext/frmtmb.sample--test-nl-body-draws.txt
ONLY='test-nl-body-draws' bash dev/nlfix-suite.sh ext 0 1 frmtmb.sample > dev/nlfix-suite-ext.log 2>&1
rm -f dev/nlfix-suite-core/frmtmb--test-compat*.txt dev/nlfix-suite-core/frmtmb--test-nl*.txt \
  dev/nlfix-suite-core/frmtmb--test-special-names.txt
ONLY='test-compat|test-nl|test-special-names' bash dev/nlfix-suite.sh core 0 8 frmtmb > dev/nlfix-suite-core.log 2>&1
# before-arm evidence on rellib-r16, the new files only
mkdir -p dev/nlfix-suite-r16
( export NLFIX_LIB=C:/Users/adf44/source/r/rellib-r16
  "$R" dev/nlfix-testrun.R frmtmb tests/testthat/test-nl-family.R dev/nlfix-suite-r16/family.txt > /dev/null 2>&1 &
  "$R" dev/nlfix-testrun.R frmtmb tests/testthat/test-nl-stan-funs.R dev/nlfix-suite-r16/stanfuns.txt > /dev/null 2>&1 &
  FRMTMB_BRMS_FIT_TESTS=true "$R" dev/nlfix-testrun.R frmtmb tests/testthat/test-nl-brms.R dev/nlfix-suite-r16/nlbrms.txt > /dev/null 2>&1 &
  "$R" dev/nlfix-testrun.R frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-nl-body-draws.R dev/nlfix-suite-r16/sample.txt > /dev/null 2>&1 &
  wait )
# the logs the report reads
export NLFIX_LIB="C:/Users/adf44/source/r/wt-nlfix-lib;C:/Users/adf44/source/r/rellib-r16"
"$R" dev/nlfix-d1check.R > dev/nlfix-d1check.log 2>&1
"$R" dev/nlfix-startcheck.R > dev/nlfix-startcheck.log 2>&1
"$R" dev/nlfix-extfam.R > dev/nlfix-extfam.log 2>&1
"$R" dev/nlfix-asym.R > dev/nlfix-asym.log 2>&1
"$R" dev/nlfix-sar.R > dev/nlfix-sar.log 2>&1
"$R" dev/nlfix-arnl2.R > dev/nlfix-arnl2.log 2>&1
"$R" dev/nlfix-cost.R > dev/nlfix-cost.log 2>&1
"$R" dev/nlfix-u955.R > dev/nlfix-u955.log 2>&1
unset NLFIX_LIB
# gated tier: every file that skipped ungated
gated=$(grep -l "SKIP=[1-9]" dev/nlfix-suite-core/*.txt | sed -E 's/.*frmtmb--(.*)[.]txt/\1/' | paste -sd'|')
mkdir -p dev/nlfix-suite-gated
ONLY="($gated)[.]R" bash dev/nlfix-suite.sh gated 1 8 frmtmb > dev/nlfix-suite-gated.log 2>&1
gateds=$(grep -l "SKIP=[1-9]" dev/nlfix-suite-ext/frmtmb.sample--*.txt | sed -E 's/.*frmtmb.sample--(.*)[.]txt/\1/' | paste -sd'|')
if [ -n "$gateds" ]; then
  ONLY="($gateds)[.]R" bash dev/nlfix-suite.sh gatedsample 1 8 frmtmb.sample > dev/nlfix-suite-gatedsample.log 2>&1
fi
# the sweep and the equiv pairs
mkdir -p dev/nlfix-equiv-out
FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-nlfix-lib;C:/Users/adf44/source/r/rellib-r16" \
  "$R" dev/equiv-nl-sweep.R dev/nlfix-equiv-out > dev/nlfix-sweep.log 2>&1
base=/c/Users/adf44/source/r/frmtmb-wt-r16base/dev/r16base-equiv-out
cp "$base"/p*.tsv dev/nlfix-equiv-out/
cp "$base"/p*.rds dev/nlfix-equiv-out/ 2>/dev/null
for id in $(cat dev/nlfix-equiv-ids.txt); do rm -f dev/nlfix-equiv-out/$id.tsv dev/nlfix-equiv-out/$id.rds; done
FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-nlfix-lib;C:/Users/adf44/source/r/rellib-r16" \
  bash dev/equiv-run.sh 8 dev/nlfix-equiv-out > dev/nlfix-equiv-run.log 2>&1
"$R" dev/nlfix-equiv-compare.R > dev/nlfix-equiv-compare.log 2>&1
echo FINAL1 DONE > dev/nlfix-final1.done
