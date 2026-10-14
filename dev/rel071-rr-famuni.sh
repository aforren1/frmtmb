#!/usr/bin/env bash
# Lane famuni's 720-fit check of gev_fit_check() (dev/famuni-r2-xiedge.R,
# both links, and dev/famuni-r2-xiedge-sum.R) rerun on a build of the
# release tree. The lane's scripts read and write its own worktree, so
# they run here as copies with the lane library replaced by <lib> and
# dev/famuni-out/ by the rerun directory. Every log is compared with
# the lane's, line by line.
#
#   bash dev/rel071-rr-famuni.sh <lib> <tag>
set -u
LIB="$1"; TAG="$2"
REL=/c/Users/adf44/source/r/frmtmb-wt-release
LANE=/c/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out
W=$REL/dev/rel071-log/reruns/$TAG/famuni
WW=${W/\/c\//C:/}
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rm -rf "$W"; mkdir -p "$W"
for s in r2-xiedge r2-xiedge-sum; do
  sed -e "s|C:/Users/adf44/source/r/wt-famuni-lib|$LIB|g" \
      -e "s|C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out|$WW|g" \
      "$REL/dev/famuni-$s.R" > "$W/famuni-$s.R"
done
grep -n "libPaths" -A1 "$W/famuni-r2-xiedge.R" | head -3
cd "$REL"
"$RS" "$W/famuni-r2-xiedge.R" log1p > "$W/r2-xiedge-log1p.txt" 2>&1 &
"$RS" "$W/famuni-r2-xiedge.R" identity > "$W/r2-xiedge-identity.txt" 2>&1 &
wait
"$RS" "$W/famuni-r2-xiedge-sum.R" > "$W/r2-xiedge-sum.txt" 2>&1
for f in r2-xiedge-log1p.txt r2-xiedge-identity.txt r2-xiedge-sum.txt; do
  if diff -q "$LANE/$f" "$W/$f" > /dev/null; then
    echo "SAME as the lane: $f ($(grep -c '' "$W/$f") lines)"
  else
    echo "DIFFERENT from the lane: $f"
    diff "$LANE/$f" "$W/$f" | head -40
  fi
done
