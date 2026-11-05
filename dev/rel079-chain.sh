#!/usr/bin/env bash
# The verification of the final build, one step after another, about 8 processes
# at a time. Each step's log is named below, and dev/rel079-log/chain.txt has a line
# per step. Do not stop this shell: the surviving R session then cannot start a
# child R, with status 322 (dev/reviews/2026-11-02-nlfix.md, item 10).
#
#   bash dev/rel079-chain.sh [step ...]    # default all
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel079
cd "$ROOT" || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
L=dev/rel079-log
C=$L/chain.txt
note() { echo "$(date '+%F %T') $*" >> "$C"; }
steps="${*:-plain ref sessions focus-base obbuild ob check floor}"
for s in $steps; do
  note "$s start"
  case "$s" in
    plain|ref) bash dev/rel079-verify.sh "$s" 8;;
    sessions)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.79.0 4 > $L/sessions-plain.txt 2>&1
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.79.0g 6 > $L/sessions-gated.txt 2>&1;;
    focus-base) bash dev/rel079-verify.sh base 8;;
    obbuild)
      bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/openblas-build.txt 2>&1
      bash dev/ciharden-openblas.sh 0.3.32 lapack >> $L/openblas-build.txt 2>&1;;
    ob) cp $L/focus.jobs $L/ob.jobs; bash dev/rel079-verify.sh ob 2;;
    check) bash dev/rel079-check.sh frmtmb frmtmb.sample frmtmb.learn frmtmb.eam frmtmb.survival frmtmb.spatial > $L/check-1.txt 2>&1;;
    floor) bash dev/rel079-floor.sh 8 > $L/floor-1.txt 2>&1;;
    *) note "unknown step $s";;
  esac
  note "$s end, exit $?"
done
note "CHAIN END"
