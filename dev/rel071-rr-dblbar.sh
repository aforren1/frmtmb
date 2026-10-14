#!/usr/bin/env bash
# Lane dblbar's case-by-case brms comparison (dev/dblbar-cases.R through
# dev/dblbar-frm-probe.R, dev/dblbar-compare.R and dev/dblbar-table.R)
# rerun on a build of the release tree. The lane's scripts read and
# write its own worktree, so they run here as copies with three
# substitutions: the lane library becomes <lib>, the worktree the
# release tree, and dev/dblbar-log/ the rerun directory. brms's probe
# (brms 2.23.0, no frmtmb involved) and the base arm (rellib-r8) are
# the lane's, copied in; the merged build takes the lane arm's place.
# The table must equal the lane's line for line.
#
#   bash dev/rel071-rr-dblbar.sh <lib> <tag>
set -u
LIB="$1"; TAG="$2"
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LANE=/c/Users/adf44/source/r/frmtmb-wt-dblbar
W=$REL/dev/rel071-log/reruns/$TAG/dblbar
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rm -rf "$W"; mkdir -p "$W/dev/dblbar-log" "$W/scripts"
for f in brms-probe.rds frm-probe-base.rds compare-base.rds; do
  cp "$LANE/dev/dblbar-log/$f" "$W/dev/dblbar-log/"
done
cp "$REL/dev/dblbar-data.R" "$REL/dev/dblbar-cases.R" "$W/dev/"
for s in frm-probe compare table; do
  sed -e "s|C:/Users/adf44/source/r/wt-dblbar-lib|$LIB|g" \
      -e "s|C:/Users/adf44/source/r/frmtmb-wt-dblbar|${W/\/c\//C:/}|g" \
      "$REL/dev/dblbar-$s.R" > "$W/scripts/dblbar-$s.R"
done
cd "$REL"
"$RS" "$W/scripts/dblbar-frm-probe.R" lane > "$W/frm-probe.out" 2>&1
"$RS" "$W/scripts/dblbar-compare.R" lane > "$W/compare.out" 2>&1
"$RS" "$W/scripts/dblbar-table.R" > "$W/table.out" 2>&1
head -1 "$W/dev/dblbar-log/frm-probe-lane.txt"
for f in table.txt compare-lane.txt; do
  if diff -q "$LANE/dev/dblbar-log/$f" "$W/dev/dblbar-log/$f" > /dev/null; then
    echo "SAME as the lane: $f"
  else
    echo "DIFFERENT from the lane: $f"
    diff "$LANE/dev/dblbar-log/$f" "$W/dev/dblbar-log/$f" | head -40
  fi
done
# the frmtmb probe after its lib: line
if diff -q <(grep -v "^lib:" "$LANE/dev/dblbar-log/frm-probe-lane.txt") \
           <(grep -v "^lib:" "$W/dev/dblbar-log/frm-probe-lane.txt") > /dev/null; then
  echo "SAME as the lane: frm-probe-lane.txt (after the lib: lines)"
else
  echo "DIFFERENT from the lane: frm-probe-lane.txt"
  diff <(grep -v "^lib:" "$LANE/dev/dblbar-log/frm-probe-lane.txt") \
       <(grep -v "^lib:" "$W/dev/dblbar-log/frm-probe-lane.txt") | head -40
fi
