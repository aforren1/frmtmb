#!/usr/bin/env bash
# The recorded survey of the 0.80.0 release, in the brief's order. First the
# compat frame, which is gitignored. The sample was drawn from 0.77.0's
# frame (rellib-r15), and that library is gone, so the frame the summary
# reads is 0.79.0's (rellib-r17), the nearest build that still exists: at
# 0.79.0 the corrected-weights line read 89.0 on 0.77.0's frame and 88.9 on
# 0.79.0's. This build's frame is dumped beside it to compare. Then all 172
# pairs into dev/rel080-equiv-out with RECORD=1, 3 at a time (the round's
# limit of three R processes; the brief's command says 8).
#
#   bash dev/rel080-survey.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel080/dev/rel080-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel080 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
"$RS" dev/equiv-dump-compat.R "$L/rellib-r17" dev/equiv-frame-out > dev/rel080-log/frame-r17.txt 2>&1
"$RS" dev/equiv-dump-compat.R "$L/rellib-r19" dev/rel080-log/frame-r19 > dev/rel080-log/frame-r19.txt 2>&1
rm -rf dev/rel080-equiv-out
FRMTMB_EQUIV_LIB="$L/rellib-r19;$L/rellib-r18" RECORD=1 bash dev/equiv-run.sh "${EQUIV_JOBS:-3}" dev/rel080-equiv-out \
  > dev/rel080-log/survey.txt 2>&1
echo "SURVEY END $?" >> dev/rel080-log/survey.txt
