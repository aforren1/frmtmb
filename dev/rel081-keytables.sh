#!/usr/bin/env bash
# The lanes' key tables on the merged build (rellib-r20), each by the
# lane's or reviewer's own script, ONE R process at a time, so this runs
# beside the survey (2 processes) under the round's limit of three.
# Lane libraries in a script read rellib-r20 through dev/rel081-rerun.R;
# a "before" arm keeps rellib-r19. Logs: dev/rel081-log/kt/.
#
#   bash dev/rel081-keytables.sh [code8|support|hollow|all]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
O=dev/rel081-log/kt
mkdir -p $O dev/code8-out dev/review-code8-out dev/support-out
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/stan-cache
export NOT_CRAN=true
run() { local n="$1"; shift; echo "$(date +%T) $n" >> $O/progress.txt
        "$RS" dev/rel081-rerun.R "$@" > "$O/$n.txt" 2>&1; }
note() { echo "$(date '+%F %T') KT $*" >> dev/rel081-log/chain.txt; }
what="${1:-all}"
if [ "$what" = code8 ] || [ "$what" = all ]; then
  note "code8 start"
  CODE8_ARM=after run c8-battery-base dev/code8-battery.R after-base
  CODE8_ARM=after run c8-battery-ywx dev/code8-battery.R after-ywx ywx
  CODE8_ARM=after run c8-short dev/code8-short.R after
  CODE8_ARM=after run c8-short-x dev/code8-short.R after-x 1,2,3,4,5,6,7,8 1e-3,3e-3,0.03
  run c8-summarize dev/code8-summarize.R
  run c8-short-bound dev/code8-short-bound.R after after-x
  CODE8_ARM=after run c8-ridge dev/code8-ridge-exp.R none rel
  CODE8_ARM=after run c8-fix1 dev/code8-fix1-draws.R rel
  RV_ARM=after run c8-ridge-seeds dev/review-code8-ridge-seeds.R after
  RV_ARM=after run c8-dnull3-after dev/review-code8-fix-dnull3.R after
  RV_ARM=before run c8-dnull3-before dev/review-code8-fix-dnull3.R before
  note "code8 end"
fi
if [ "$what" = support ] || [ "$what" = all ]; then
  note "support start"
  export FRMTMB_PORT_LIB=$L/rellib-r20
  for s in is is2 nq; do run sp-rev-$s dev/review-support-$s.R; done
  # the bitwise fits, pinned to core 0 as the lane pinned both its arms
  # (dev/rtmb-pitfalls.md, item 21)
  cmd //c "dev\\support-pin.bat dev\\support-fitsame-re.R $(cygpath -w $O)\\sp-fitsame-re.txt"
  cmd //c "dev\\support-pin.bat dev\\support-fitsame.R $(cygpath -w $O)\\sp-fitsame.txt"
  run sp-ndt dev/support-ndt-probe.R
  run sp-gev349 dev/support-gev349-probe.R
  unset FRMTMB_PORT_LIB
  note "support end"
fi
if [ "$what" = hollow ] || [ "$what" = all ]; then
  note "hollow start"
  run ho-probe6 dev/review-hollow-probe6.R
  note "hollow end"
fi
echo "KEYTABLES END" >> $O/progress.txt
