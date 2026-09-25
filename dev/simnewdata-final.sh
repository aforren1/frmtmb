#!/bin/sh
# Re-run every lane-arm measurement against the FINAL lane build, so
# each log postdates the last source edit. The base arms are not re-run:
# rellib-r3 does not change.
#   sh dev/simnewdata-final.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-simnewdata
cd "$ROOT" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/simnewdata-log
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$R" dev/simnewdata-probe2.R > $L/probe2-lane.txt 2>&1
"$R" dev/simnewdata-invariant.R > $L/invariant.txt 2>&1
"$R" dev/simnewdata-bitwise.R lane > $L/bitwise-lane.txt 2>&1
"$R" dev/simnewdata-bitwise.R compare > $L/bitwise-compare.txt 2>&1
"$R" dev/simnewdata-callers.R lane > $L/callers-lane.txt 2>&1
"$R" dev/simnewdata-callers.R compare > $L/callers-compare.txt 2>&1
sh dev/simnewdata-rplots.sh lane test-draws-methods.R test-brms-suite-methods.R \
  > /dev/null 2>&1
FRMTMB_PORT_ROOT=$ROOT FRMTMB_PORT_LIB=C:/Users/adf44/source/r/simnewdata-lib \
  FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache sh dev/brmsport-record.sh \
  > /dev/null 2>&1
"$R" dev/brmsport-ledger.R > $L/ledger-final.txt 2>&1
{
  echo "# RESULT lines of the lane record, from dev/brmsport-log/run-*.txt"
  for f in dev/brmsport-log/run-*.txt; do
    printf "%s  " "$(date -r "$f" '+%Y-%m-%d %H:%M:%S')"
    grep -h "^RESULT" "$f" || echo "NO RESULT $f"
  done
} > $L/record.txt
echo FINAL DONE
