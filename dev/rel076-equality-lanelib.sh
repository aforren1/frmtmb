#!/usr/bin/env bash
# dev/rel076-equality.sh found three lane logs that differ from the
# merged build's run: mocode's edge-lane.txt (06:57) and acprior's
# agree.txt and agree-up.txt (07:18, 07:22), each written before its
# lane's fix rounds (the libraries were last built at 10:19 and 11:03).
# This reruns those arms on each lane's own final library, read only,
# from the same copies with the library swapped, and compares them with
# the merged build's output: equal means the merge changed nothing.
#   bash dev/rel076-equality-lanelib.sh > dev/rel076-log/equality-lanelib.txt 2>&1
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel076
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
E=dev/rel076-log/eq
X=$E/lanelib; rm -rf "$X"; mkdir -p "$X/acprior"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache NOT_CRAN=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
sed "s#C:/Users/adf44/source/r/rellib-r14#C:/Users/adf44/source/r/wt-mocode-lib#" \
  "$E/mocode-edge.R" > "$X/mocode-edge.R"
for s in acprior-agree.R acprior-agree-up.R; do
  sed -e "s#C:/Users/adf44/source/r/rellib-r14#C:/Users/adf44/source/r/wt-acprior-lib#" \
      -e "s#\"dev/rel076-log/eq\", \"acprior\"#\"$X\", \"acprior\"#" \
      "$E/$s" > "$X/$s"
done
cp "$E"/acprior/*-brms.rds "$E"/acprior/agree-before.rds \
   "$E"/acprior/agree-up-before.rds "$X/acprior/"
"$RS" "$X/mocode-edge.R" lane > "$X/edge-lane.txt" 2>&1 &
"$RS" "$X/acprior-agree.R" after > "$X/agree-after.log" 2>&1 &
"$RS" "$X/acprior-agree-up.R" after > "$X/agree-up-after.log" 2>&1 &
wait
"$RS" "$X/acprior-agree.R" summary > "$X/agree.txt" 2>&1
"$RS" "$X/acprior-agree-up.R" summary > "$X/agree-up.txt" 2>&1
for p in "mocode/edge-lane.txt edge-lane.txt" "acprior/agree.txt agree.txt" \
         "acprior/agree-up.txt agree-up.txt"; do
  set -- $p
  grep -vE 'lib|== frmtmb' "$E/$1" | tr -d '\r' > "$X/a.tmp"
  grep -vE 'lib|== frmtmb' "$X/$2" | tr -d '\r' > "$X/b.tmp"
  if cmp -s "$X/a.tmp" "$X/b.tmp"; then
    echo "SAME    merged $E/$1 = lane library $X/$2 ($(wc -l < "$X/a.tmp") lines)"
  else
    echo "DIFFERS merged $E/$1 against lane library $X/$2"
    diff "$X/b.tmp" "$X/a.tmp" | head -20
  fi
done
grep -h "== frmtmb\|lib" "$X/edge-lane.txt" | head -2
echo "LANELIB DONE"
