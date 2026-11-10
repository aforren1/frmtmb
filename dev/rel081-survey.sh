#!/usr/bin/env bash
# The recorded survey of the 0.81.0 release. First the compat frames,
# which are gitignored: the summary's frame is 0.79.0's (rellib-r17), as
# at 0.80.0, so the corrected-weights line reads like for like; 0.80.0's
# (rellib-r19) and this build's (rellib-r20) are dumped beside it for
# lane compat's frame diff. The live row the flow now scores against is
# read by dev/equiv-live-row.R from FRMTMB_EQUIV_LIB, not from a dump.
# Then all 172 pairs into dev/rel081-equiv-out with RECORD=1, ONE pair at
# a time: each pair is 2 R processes (the measuring process and its
# live-row child), and the round's limit is three.
#
#   bash dev/rel081-survey.sh [frames|survey]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
G=dev/rel081-log
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
case "${1:-all}" in
  frames|all)
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r17" dev/equiv-frame-out > $G/frame-r17.txt 2>&1
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r19" $G/frame-r19 > $G/frame-r19.txt 2>&1
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r20" $G/frame-r20 > $G/frame-r20.txt 2>&1
    "$RS" dev/compat-frame-diff.R $G/frame-r19/compat.tsv $G/frame-r20/compat.tsv \
      > $G/frame-diff-r19-r20.txt 2>&1
    echo "FRAMES END";;
esac
case "${1:-all}" in
  survey|all)
    rm -rf dev/rel081-equiv-out
    FRMTMB_EQUIV_LIB="$L/rellib-r20;$L/rellib-r19" RECORD=1 bash dev/equiv-run.sh 1 dev/rel081-equiv-out \
      > $G/survey.txt 2>&1
    echo "SURVEY END $?" >> $G/survey.txt;;
esac
