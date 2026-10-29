#!/usr/bin/env bash
# The consolidation's new tests, run from the release tree on the builds
# that lack them, every gate on, one file per process:
#   - test-prior-nlpar-intercept.R (release fix R1 appended its blocks,
#     dev/rel078-fix-r1-test.R): on meme's own library, which has the
#     order-dependent rule, in front of rellib-r15;
#   - test-rescor-me.R (the cross-lane test, dev/rel078-fix-rescorme-test.R):
#     on meme's library (no joint draw) and on rescorpp's (no me() prior
#     classes), each in front of rellib-r15;
#   - both on rellib-r15 alone (0.77.0).
# The lane libraries are read only. A block that fails here and passes on
# rellib-r16 is seen to fail; the per-file DETAIL lines name the blocks.
# (The first run, seenfail-1.txt, was the 0.77.0 script before this
# rewrite; its lane libraries no longer exist, and it is not evidence.)
#   bash dev/rel078-seenfail.sh
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
cd /c/Users/adf44/source/r/frmtmb-wt-rel078 || exit 1
L=dev/rel078-log
S=frmtmb.sample
ST=extensions/frmtmb.sample/tests/testthat
R15=C:/Users/adf44/source/r/rellib-r15
printf '%s\n' "$S $ST/test-prior-nlpar-intercept.R" "$S $ST/test-rescor-me.R" \
  > $L/sf-meme.jobs
printf '%s\n' "$S $ST/test-rescor-me.R" > $L/sf-rescorpp.jobs
cp $L/sf-meme.jobs $L/sf-r15.jobs
for l in meme rescorpp; do
  [ -d C:/Users/adf44/source/r/wt-$l-lib/frmtmb.sample ] ||
    { echo "no lane library for $l"; exit 1; }
  ROUND_BASE=$R15 bash dev/ciharden-scan.sh r078-sf-$l ref \
    C:/Users/adf44/source/r/wt-$l-lib 2 $L/sf-$l.jobs > $L/sf-$l.txt 2>&1 &
done
ROUND_BASE=$R15 bash dev/ciharden-scan.sh r078-sf-r15 ref base 2 \
  $L/sf-r15.jobs > $L/sf-r15.txt 2>&1 &
wait
for r in meme rescorpp r15; do
  for f in dev/ciharden-log/r078-sf-$r-suite/*.txt; do
    echo "== $r $(basename "$f"): $(grep -h '^RESULT' "$f")"
    grep -h '^lib:' "$f" | head -1
    grep -h '^DETAIL	\(FAIL\|ERROR\)' "$f" | cut -f2-4 | sort | uniq -c | sed 's/^/   /'
  done
done
echo "SEENFAIL DONE"
