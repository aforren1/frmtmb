#!/usr/bin/env bash
# The recorded survey of the 0.79.0 release, in the brief's order. First the
# compat frame, which is gitignored and was lost: without it the corrected-weights
# line reads 83.4 for no reason. Two frames are dumped: 0.77.0's, which the sample
# was drawn from, in the place the summary reads it, and this build's, to compare.
# Then all 172 pairs into dev/rel079-equiv-out with RECORD=1.
#
#   bash dev/rel079-survey.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel079 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
export TMP=$L/../../AppData/Local/Temp/2 TEMP=$L/../../AppData/Local/Temp/2
"$RS" dev/equiv-dump-compat.R "$L/rellib-r15" dev/equiv-frame-out > dev/rel079-log/frame-r15.txt 2>&1
"$RS" dev/equiv-dump-compat.R "$L/rellib-r17" dev/rel079-log/frame-r17 > dev/rel079-log/frame-r17.txt 2>&1
rm -rf dev/rel079-equiv-out
FRMTMB_EQUIV_LIB="$L/rellib-r17;$L/rellib-r16" RECORD=1 bash dev/equiv-run.sh 8 dev/rel079-equiv-out \
  > dev/rel079-log/survey.txt 2>&1
echo "SURVEY END $?" >> dev/rel079-log/survey.txt
