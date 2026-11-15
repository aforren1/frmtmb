#!/usr/bin/env bash
# The recorded survey of the 0.83.0 release. First the compat frames of
# rellib-r21 (0.82.0) and of this build, and their diff (frm_compat() rows
# that moved this round), then every pair with RECORD=1, one pair at a
# time: each pair runs a live-row child, so it holds 2 R processes, which
# leaves one for the work beside it. Adapted from dev/rel081-survey.sh.
#
#   bash dev/rel083-survey.sh [frames|survey|all]
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel083 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
G=dev/rel083-log
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
case "${1:-all}" in
  frames|all)
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r21" $G/frame-r21 > $G/frame-r21.txt 2>&1
    "$RS" dev/equiv-dump-compat.R "$L/rellib-r22" $G/frame-r22 > $G/frame-r22.txt 2>&1
    "$RS" dev/compat-frame-diff.R $G/frame-r21/compat.tsv $G/frame-r22/compat.tsv \
      > $G/frame-diff-r21-r22.txt 2>&1
    echo "FRAMES END";;
esac
case "${1:-all}" in
  survey|all)
    rm -rf dev/rel083-equiv-out
    FRMTMB_EQUIV_LIB="$L/rellib-r22;$L/rellib-r21" RECORD=1 bash dev/equiv-run.sh 1 dev/rel083-equiv-out \
      > $G/survey.txt 2>&1
    echo "SURVEY END $?" >> $G/survey.txt;;
esac
