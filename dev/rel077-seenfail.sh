#!/usr/bin/env bash
# The consolidation's new tests (dev/rel077-fix-*-test.R, appended by
# dev/rel077-release-fix.R), run from the release tree on the builds
# that lack the fix: each lane's own library (read only, in front of
# rellib-r15) and rellib-r14 (0.76.0). Every gate on, one file per
# process. A block that fails here and passes on rellib-r15 is seen to
# fail; the per-file DETAIL lines name the blocks.
#   bash dev/rel077-seenfail.sh
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
cd /c/Users/adf44/source/r/frmtmb-wt-rel077 || exit 1
L=dev/rel077-log
C=frmtmb; S=frmtmb.sample
T=tests/testthat; ST=extensions/frmtmb.sample/tests/testthat
printf '%s\n' "$C $T/test-autocor-sigma.R" > $L/sf-sigmainn.jobs
printf '%s\n' "$C $T/test-autocor-sigma.R" > $L/sf-mointerop.jobs
printf '%s\n' "$C $T/test-prior-one-slot.R" "$S $ST/test-flat-prior.R" > $L/sf-reslot.jobs
printf '%s\n' "$S $ST/test-prior-draws.R" > $L/sf-priordraws.jobs
cat $L/sf-sigmainn.jobs $L/sf-reslot.jobs $L/sf-priordraws.jobs > $L/sf-r14.jobs
for l in sigmainn mointerop reslot priordraws; do
  ROUND_BASE=C:/Users/adf44/source/r/rellib-r15 bash dev/ciharden-scan.sh \
    r077-sf-$l ref C:/Users/adf44/source/r/wt-$l-lib 2 $L/sf-$l.jobs > $L/sf-$l.txt 2>&1 &
done
ROUND_BASE=C:/Users/adf44/source/r/rellib-r14 bash dev/ciharden-scan.sh \
  r077-sf-r14 ref base 4 $L/sf-r14.jobs > $L/sf-r14.txt 2>&1 &
wait
for r in sigmainn mointerop reslot priordraws r14; do
  for f in dev/ciharden-log/r077-sf-$r-suite/*.txt; do
    echo "== $r $(basename "$f"): $(grep -h '^RESULT' "$f")"
    grep -h '^lib:' "$f" | head -1
    grep -h '^DETAIL	\(FAIL\|ERROR\)' "$f" | cut -f2-4 | sort | uniq -c | sed 's/^/   /'
  done
done
echo "SEENFAIL DONE"
