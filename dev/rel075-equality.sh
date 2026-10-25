#!/usr/bin/env bash
# Each lane's equality and brms-comparison scripts, rerun on the merged
# build. Adapted from dev/rel074-equality.sh. The "before" arm reads
# rellib-r12 (frmtmb 0.74.0), the "after" arm rellib-r13 (0.75.0). A
# lane script that names its lane library, its worktree or its output
# directory runs from a copy in dev/rel075-log/eq/ with those strings
# replaced (sed prints each replaced line); the reviewer's sample of
# smallfix takes its libraries as an argument and runs from the lane's
# tree, unedited. Every arm is one process, so a difference is a
# difference of builds (dev/rtmb-pitfalls.md item 21); smallfix's arms
# are pinned to CPUs 0-1 as its review pinned them.
#
#   bash dev/rel075-equality.sh fast    # bhaz, sarfcor, brmsnames table,
#                                       # smallfix sample
#   bash dev/rel075-equality.sh suite   # brmsnames: names of every fit
#                                       # of the core and sample suites
#   (each > dev/rel075-log/equality-<stage>.txt 2>&1)
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel075
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R12=C:/Users/adf44/source/r/rellib-r12
R13=C:/Users/adf44/source/r/rellib-r13
LW=C:/Users/adf44/source/r/frmtmb-wt
E=dev/rel075-log/eq
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel075/dev/stan-cache
export NOT_CRAN=true
copy() {   # <lane script> <copy>: lane library -> rellib-r13, lane tree
           # -> this tree, lane output dirs -> dev/rel075-log/eq/
  sed -e "s#C:/Users/adf44/source/r/wt-[a-z]*-lib#$R13#g" \
      -e "s#C:/Users/adf44/source/r/frmtmb-wt-[a-z]*#C:/Users/adf44/source/r/frmtmb-wt-rel075#g" \
      -e "s#dev/sarfcor-out#$E/sarfcor#g" \
      -e "s#dev/brmsnames-1024-out#$E/bn-out#g" \
      -e "s#dev/brmsnames-1024-suite#$E/bn-suite#g" "dev/$1" > "$E/$2"
  diff "dev/$1" "$E/$2" | grep '^>' | sed "s/^/   $2: /"
}
case "${1:-fast}" in
fast)
  rm -rf "$E"; mkdir -p "$E/bhaz" "$E/sarfcor" "$E/bn-out" "$E/smallfix"
  echo "== bhaz: unchanged cox() fits, bitwise (dev/bhaz-bitwise.R)"
  copy bhaz-bitwise.R bhaz-bitwise.R
  for i in 1 2; do
    "$RS" "$E/bhaz-bitwise.R" base "$E/bhaz/base-$i.rds" > "$E/bhaz-base-$i.log" 2>&1 &
    "$RS" "$E/bhaz-bitwise.R" lane "$E/bhaz/lane-$i.rds" > "$E/bhaz-merged-$i.log" 2>&1 &
  done
  echo "== sarfcor: brms log densities, and newdata against brms"
  for s in sarfcor-lpcheck.R sarfcor-f3-epred.R sarfcor-fr2-rows.R; do
    copy "$s" "$s"
  done
  ( "$RS" "$E/sarfcor-lpcheck.R" > "$E/sarfcor/lpcheck.txt" 2>&1
    "$RS" "$E/sarfcor-f3-epred.R" > "$E/sarfcor/f3-epred.txt" 2>&1
    "$RS" "$E/sarfcor-fr2-rows.R" frmtmb > "$E/sarfcor/fr2-rows-after.txt" 2>&1 ) &
  wait
  echo "== brmsnames: the name table, both arms (dev/brmsnames-1024-frm.R)"
  cp "$LW-brmsnames"/dev/brmsnames-1024-out/brms-*.rds "$E/bn-out/"
  copy brmsnames-1024-frm.R bn-frm.R
  copy brmsnames-1024-table.R bn-table.R
  ids=$(seq 1 40)
  ( "$RS" "$E/bn-frm.R" base $ids > "$E/bn-frm-base.log" 2>&1 ) &
  ( "$RS" "$E/bn-frm.R" lane $ids > "$E/bn-frm-merged.log" 2>&1 ) &
  echo "== smallfix: intact fits through keep_missing(), pinned to CPUs 0-1"
  cmd //c "dev\\rel075-eq-smallfix.bat"
  wait
  for l in "$E"/*.log; do
    echo "-- $l: $(grep -E 'arm|frmtmb from|lib:' "$l" | head -2 | tr '\n' ' ')"
  done
  echo "== bhaz compare"
  "$RS" dev/bhaz-bitwise-cmp.R "$E/bhaz" 2>&1 | awk '{print}' | sort | uniq -c | sort -rn | head -40
  echo "== brmsnames table"
  "$RS" "$E/bn-table.R" > "$E/bn-table.md" 2>&1
  echo "table rows: $(grep -c '^|' "$E/bn-table.md")"
  echo "== smallfix compare"
  "$RS" dev/rel075-eq-smallfix.R 2>&1
  echo "EQUALITY FAST DONE";;
suite)
  rm -rf "$E/bn-suite" "$E/bn-suitelog"; mkdir -p "$E/bn-suitelog"
  copy brmsnames-1024-suitenames.R bn-suitenames.R
  copy brmsnames-1024-suitediff.R bn-suitediff.R
  N="${2:-8}"
  for pk in frmtmb frmtmb.sample; do
    dir=tests/testthat; [ "$pk" = frmtmb ] || dir=extensions/$pk/tests/testthat
    for f in "$dir"/test-*.R; do
      while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
      "$RS" "$E/bn-suitenames.R" lane "$pk" "$f" \
        > "$E/bn-suitelog/$pk--$(basename "$f" .R).txt" 2>&1 &
    done
    wait
  done
  echo "suite files recorded: $(ls "$E"/bn-suite/lane/*/*.rds | wc -l)"
  "$RS" dev/rel075-eq-suitenames.R 2>&1
  echo "EQUALITY SUITE DONE";;
esac
