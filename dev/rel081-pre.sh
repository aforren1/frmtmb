#!/usr/bin/env bash
# Before the recorded survey: remeasure a handful of each equivalence
# lane's pairs on the merged build and instrument (translator: p021,
# p047, p048, p131, p136; compat: p154, p010 and p146, two pairs the
# ordinal-only rule reaches, p167, p169), one pair at a time (each pair
# is 2 R processes with the live-row child), and compare each with the
# lane's own outdir (dev/rel081-pre-cmp.R).
#
#   bash dev/rel081-pre.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export NOT_CRAN=true
export FRMTMB_EQUIV_LIB="$L/rellib-r20;$L/rellib-r19"
O=dev/rel081-pre-equiv-out
mkdir -p $O/logs
T=$L/frmtmb-wt-translator/dev/translator-noreg3-out
C=$L/frmtmb-wt-compat/dev/compat-noreg2-edit-out
for id in p021 p047 p048 p131 p136 p154 p010 p146 p167 p169; do
  "$RS" dev/equiv-measure.R $id $O --force > $O/logs/$id.log 2>&1 ||
    echo "FAILED $id exit $?" >> $O/logs/$id.log
done
"$RS" dev/rel081-pre-cmp.R $O p021=$T p047=$T p048=$T p131=$T p136=$T \
  p154=$C p010=$C p146=$C p167=$C p169=$C
grep -h '^lib:' $O/logs/*.log | sort | uniq -c
echo "PRE END"
