#!/usr/bin/env bash
# The second verification, after the release review's fixes. Each step writes a line to
# dev/rel079-log/chain2.txt; the wait is on the line "CHAIN2 END <steps>", printed at the
# end of this script. Do not stop this shell.
#
#   bash dev/rel079-chain2.sh <step ...>
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel079
cd "$ROOT" || exit 1
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=dev/rel079-log
C=$L/chain2.txt
LIBS=C:/Users/adf44/source/r
note() { echo "$(date '+%F %T') $*" >> "$C"; }
for s in "$@"; do
  note "$s start"
  case "$s" in
    install) bash dev/rel079-install.sh $LIBS/rellib-r17 > $L/install7.txt 2>&1;;
    survey)
      rm -rf dev/rel079-equiv-rerun-out
      FRMTMB_EQUIV_LIB="$LIBS/rellib-r17;$LIBS/rellib-r16" bash dev/equiv-run.sh 8 dev/rel079-equiv-rerun-out > $L/survey2.txt 2>&1
      "$RS" dev/rel079-equiv-diff.R dev/rel079-equiv-out dev/rel079-equiv-rerun-out > $L/survey2-diff.txt 2>&1;;
    keytables)
      xargs -P 8 -L 1 bash dev/rel079-rerun.sh < dev/rel079-rerun-final.jobs > $L/rerun-final-driver.txt 2>&1;;
    ledger) bash dev/rel079-ledger.sh > $L/ledger-run-final.txt 2>&1
      "$RS" dev/rel079-ledger-diff.R > $L/ledger-diff.txt 2>&1
      "$RS" dev/rel079-bin2-diff.R > $L/bin2-diff.txt 2>&1;;
    plain|ref) ROUND_BASE=$LIBS/rellib-r17 bash dev/rel079-verify.sh "$s" 8;;
    sessions)
      NOT_CRAN=true bash dev/release/run-sessions.sh rel0.79.0 4 > $L/sessions-plain.txt 2>&1
      NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
        FRMTMB_SCALE_TESTS=true bash dev/release/run-sessions.sh rel0.79.0g 6 > $L/sessions-gated.txt 2>&1;;
    jobs-focus)
      git status --porcelain --untracked-files=all -- 'tests/testthat/test-*.R' 'extensions/*/tests/testthat/test-*.R' |
        cut -c4- | awk '{ if ($0 ~ /^extensions\//) { split($0,a,"/"); print a[2], $0 } else print "frmtmb", $0 }' > $L/focus.jobs;;
    focus-base) bash dev/rel079-verify.sh base 8;;
    obbuild)
      bash dev/ciharden-openblas.sh 0.3.26 lapack > $L/openblas-build.txt 2>&1
      bash dev/ciharden-openblas.sh 0.3.32 lapack >> $L/openblas-build.txt 2>&1;;
    ob) cp $L/focus.jobs $L/ob.jobs; bash dev/rel079-verify.sh ob 3;;
    check) bash dev/rel079-check.sh frmtmb frmtmb.sample frmtmb.latent frmtmb.eam frmtmb.learn > $L/check-3.txt 2>&1;;
    floor) bash dev/rel079-floor.sh 8 > $L/floor-2.txt 2>&1;;
    docs) bash dev/rel079-docs-check.sh > $L/docs-check2.txt 2>&1;;
    *) note "unknown step $s";;
  esac
  note "$s end, exit $?"
done
note "CHAIN2 END $*"
