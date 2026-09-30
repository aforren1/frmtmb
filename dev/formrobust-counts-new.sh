#!/usr/bin/env bash
# The before and after arms of the lane's new test files, from their
# logs: bash dev/formrobust-counts-new.sh >> dev/formrobust-log/counts.md
cd /c/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-log
echo "### The new test files, before (rellib-r4) and after (lane library)"
echo
echo "| file | before | after |"
echo "|---|---|---|"
for a in after--*.txt; do
  k=${a#after--}
  b="before--$k"
  ra=$(grep -a '^RESULT' "$a" | tail -1 | cut -d' ' -f3-)
  rb=$(grep -a '^RESULT' "$b" 2>/dev/null | tail -1 | cut -d' ' -f3-)
  echo "| ${k%.txt} | ${rb:-not run} | $ra |"
done
echo
