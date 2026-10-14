#!/usr/bin/env bash
# Lane postfix's nonlinear-start reruns on a build of the release tree,
# as its review asked (dev/reviews/2026-10-13-postfix.md, the overlap
# with optstart): dev/postfix-nlstart.R, dev/postfix-rev-nlscan.R and
# dev/postfix-fix2-nlnew.R, each compared line by line with the lane's
# log of its final build (fix round 2), after the line that names the
# library, with timings masked ("time 0.99s" reads "time Ts"). Three at
# once.
#
#   bash dev/rel071-rr-postfix.sh <lib> <tag>
set -u
LIB="$1"; TAG="$2"
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LANE=/c/Users/adf44/source/r/frmtmb-wt-postfix/dev/postfix-log
W=$REL/dev/rel071-log/reruns/$TAG/postfix
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rm -rf "$W"; mkdir -p "$W"
cd "$REL"
"$RS" dev/postfix-nlstart.R after "$LIB" > "$W/nlstart.txt" 2>&1 &
"$RS" dev/postfix-rev-nlscan.R "$LIB" > "$W/rev-nlscan.txt" 2>&1 &
"$RS" dev/postfix-fix2-nlnew.R "$LIB" > "$W/fix2-nlnew.txt" 2>&1 &
wait
mask() { sed -E 's/[0-9]+[.][0-9]+s([^a-z]|$)/Ts\1/g'; }
cmp1() {
  if diff -q <(tail -n +2 "$LANE/$1" | mask) <(tail -n +2 "$W/$2" | mask) \
       > /dev/null; then
    echo "SAME as the lane after line 1: $2 ($(grep -c '' "$W/$2") lines)"
  else
    echo "DIFFERENT from the lane: $2"
    diff <(tail -n +2 "$LANE/$1" | mask) <(tail -n +2 "$W/$2" | mask) |
      head -60
  fi
}
head -n 1 "$W/nlstart.txt" "$W/rev-nlscan.txt" "$W/fix2-nlnew.txt"
cmp1 nlstart-fix2.txt nlstart.txt
cmp1 rev-nlscan-fix2.txt rev-nlscan.txt
cmp1 fix2-nlnew.txt fix2-nlnew.txt
