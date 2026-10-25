#!/usr/bin/env bash
# The SESSION, LEAK and failing DETAIL lines of the lane's one-session
# runs (dev/release/run-sessions.sh, tiers brmsnames-splain and
# brmsnames-sgated), and the Status and problem lines of its R CMD
# checks, for the findings.
#
#   bash dev/brmsnames-1024-sesssum.sh
#
# Output: dev/brmsnames-1024-log/brmsnames-sessions-summary.txt and
# dev/brmsnames-1024-log/check-summary.txt.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
o=dev/brmsnames-1024-log/brmsnames-sessions-summary.txt
: > "$o"
for t in splain sgated; do
  l="dev/release/brmsnames-$t.log"
  [ -f "$l" ] || continue
  {
    echo "== $t"
    grep -a "^SESSION\|ran [0-9]* of" "$l"
    echo "LEAK lines: $(grep -ac 'LEAK' "$l"), untagged: $(grep -a 'LEAK' "$l" | grep -vc '\[with load of')"
    grep -a 'LEAK' "$l" | grep -v '\[with load of' | head -20
    echo "DETAIL FAIL/ERROR lines: $(grep -ac '^DETAIL \(FAIL\|ERROR\)' "$l")"
    grep -a '^DETAIL \(FAIL\|ERROR\)' "$l" | cut -c1-200 | head -20
    echo "lib lines of the sessions:"
    grep -ah "^lib" dev/release/brmsnames-$t-sessions/*.txt | sort | uniq -c
    echo
  } >> "$o"
done
c=dev/brmsnames-1024-log/check-summary.txt
: > "$c"
for p in frmtmb frmtmb.sample; do
  l="dev/brmsnames-1024-check/$p/check.log"
  [ -f "$l" ] || continue
  {
    echo "== $p"
    grep -a "^\* checking.*\(WARNING\|NOTE\|ERROR\)" "$l"
    grep -a -A4 "^\* checking.*\(WARNING\|NOTE\|ERROR\)$" "$l" | grep -v "^\*" | head -30
    echo "tests: $(grep -ah '\[ FAIL' dev/brmsnames-1024-check/$p/$p.Rcheck/tests/testthat.Rout* | tail -1)"
    grep -a "^Status" "$l"
    echo
  } >> "$c"
done
cat "$o" "$c"
