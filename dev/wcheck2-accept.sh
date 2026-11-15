#!/usr/bin/env bash
# Lane wcheck2: the acceptance runs of the brief, on rellib-r21 ("base")
# and the lane library ("lane"), 3 R processes at a time:
#   wcheck's 58-case table      dev/wcheck-cases.R <lib> <mv|cens|prior>
#   the 24 guard-absent models  dev/reviews/wcheck-review-common.R <lib>
#   the route table             dev/reviews/wcheck-short2-routes.R <lib>
# Each script takes a library path first; with a path its .libPaths()
# is that path, then rellib-r20, then the user library, so frmtmb comes
# from the path. Logs: dev/wcheck2-out/accept/<name>-<arm>.txt.
#   bash dev/wcheck2-accept.sh [arm ...]
cd "$(dirname "$0")/.."
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
OUT=dev/wcheck2-out/accept
mkdir -p $OUT
arms="${*:-base lane}"
jobs=()
for arm in $arms; do
  if [ "$arm" = base ]; then lib=C:/Users/adf44/source/r/rellib-r21; else
    lib=C:/Users/adf44/source/r/wt-wcheck2-lib; fi
  for p in mv cens prior; do
    jobs+=("dev/wcheck-cases.R $lib $p $OUT/cases-$p-$arm.txt")
  done
  jobs+=("dev/reviews/wcheck-review-common.R $lib - $OUT/common-$arm.txt")
  jobs+=("dev/reviews/wcheck-short2-routes.R $lib - $OUT/routes-$arm.txt")
done
printf '%s\n' "${jobs[@]}" | xargs -P 3 -L 1 bash -c '
  s=$0; lib=$1; part=$2; log=$3
  if [ "$part" = - ]; then part=""; fi
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" $s $lib $part > $log 2>&1
  echo "done $log"'
