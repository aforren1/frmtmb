#!/usr/bin/env bash
# Each lane's equality and brms-agreement scripts, rerun on the merged
# build. Adapted from dev/rel075-equality.sh. The "before" arm reads
# rellib-r13 (frmtmb 0.75.0), the "after" arm rellib-r14 (0.76.0). A
# lane script names its lane library, its worktree and its output
# directory, so it runs from a copy in dev/rel076-log/eq/ with those
# strings replaced (the log prints each replaced line). Every arm is one
# process, so a difference is a difference of builds (dev/rtmb-pitfalls.md
# item 21). Each output is then compared with the lane's own log of the
# same arm, line for line, past the lines that name a library.
#
# brms's arms of acprior's three sampler comparisons are not rerun: they
# do not read frmtmb, and their draws (dev/acprior-out/*-brms.rds in the
# lane's tree) are copied, as 0.75.0 copied brmsnames's stored names.
#
#   bash dev/rel076-equality.sh > dev/rel076-log/equality.txt 2>&1
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel076
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R14=C:/Users/adf44/source/r/rellib-r14
LW=/c/Users/adf44/source/r/frmtmb-wt
E=dev/rel076-log/eq
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel076/dev/stan-cache
export NOT_CRAN=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
rm -rf "$E"; mkdir -p "$E/mocode" "$E/acprior" "$E/loglik"
copy() {   # <lane script> -> $E/<lane script>
  sed -e "s#C:/Users/adf44/source/r/wt-[a-z]*-lib#$R14#g" \
      -e "s#C:/Users/adf44/source/r/frmtmb-wt-[a-z]*/#C:/Users/adf44/source/r/frmtmb-wt-rel076/#g" \
      -e "s#C:/Users/adf44/source/r/frmtmb-wt-[a-z]*\"#C:/Users/adf44/source/r/frmtmb-wt-rel076\"#g" \
      -e "s#\"dev\", \"acprior-out\"#\"$E\", \"acprior\"#g" \
      "dev/$1" > "$E/$1"
  diff "dev/$1" "$E/$1" | grep '^>' | sed "s/^/   $1: /"
}
# <ours> <lane log>: line-for-line comparison past library lines
same() {
  a=$(grep -vE '^lib|libPaths|rellib|wt-[a-z]*-lib|^\[1\] "C:' "$1" | tr -d '\r')
  b=$(grep -vE '^lib|libPaths|rellib|wt-[a-z]*-lib|^\[1\] "C:' "$2" | tr -d '\r')
  if [ "$a" = "$b" ]; then
    echo "SAME    $1 = $2 ($(echo "$a" | wc -l) lines)"
  else
    echo "DIFFERS $1 against $2:"
    echo "$b" > "$E/same-b.tmp"; echo "$a" > "$E/same-a.tmp"
    diff "$E/same-b.tmp" "$E/same-a.tmp" | head -30 | sed "s/^/      /"
  fi
}
for s in mocode-agree.R mocode-se3.R mocode-edge.R acprior-ml.R \
         acprior-frm-table.R acprior-roundtrip.R acprior-agree.R \
         acprior-agree-up.R acprior-innov.R loglik-brms-agree.R \
         loglik-check.R loglik-brms-lib.R; do
  copy "$s"
done
cp "$LW-acprior"/dev/acprior-out/agree-brms.rds \
   "$LW-acprior"/dev/acprior-out/agree-up-brms.rds \
   "$LW-acprior"/dev/acprior-out/innov-brms.rds "$E/acprior/"
# loglik's two scripts source helpers from the lane's tree; the copies
# read this tree's (the same files, copied by the merge), and the brms
# helper's copy compiles into this tree's Stan cache, not the lane's
sed -i "s#frmtmb-wt-rel076/dev/loglik-brms-lib.R#frmtmb-wt-rel076/$E/loglik-brms-lib.R#" "$E/loglik-brms-agree.R"
echo "== started $(date +%T)"
M=$E/mocode; A=$E/acprior; L=$E/loglik
for arm in base lane; do
  "$RS" "$E/mocode-agree.R" "$arm" > "$M/agree-$arm.txt" 2>&1 &
  "$RS" "$E/mocode-se3.R" "$arm" > "$M/se3-$arm.txt" 2>&1 &
  "$RS" "$E/mocode-edge.R" "$arm" > "$M/edge-$arm.txt" 2>&1 &
done
wait
for arm in before after; do
  "$RS" "$E/acprior-ml.R" "$arm" > "$A/ml-$arm.txt" 2>&1 &
  "$RS" "$E/acprior-frm-table.R" "$arm" > "$A/frm-table-$arm.txt" 2>&1 &
  "$RS" "$E/acprior-roundtrip.R" "$arm" > "$A/roundtrip-$arm.txt" 2>&1 &
done
"$RS" "$E/loglik-brms-agree.R" > "$L/loglik-brms-agree.txt" 2>&1 &
"$RS" "$E/loglik-check.R" > "$L/loglik-check.txt" 2>&1 &
wait
for arm in before after; do
  "$RS" "$E/acprior-agree.R" "$arm" > "$A/agree-$arm.log" 2>&1 &
  "$RS" "$E/acprior-agree-up.R" "$arm" > "$A/agree-up-$arm.log" 2>&1 &
done
"$RS" "$E/acprior-innov.R" frm before > "$A/innov-frm-before.log" 2>&1 &
"$RS" "$E/acprior-innov.R" frm > "$A/innov-frm-after.log" 2>&1 &
wait
"$RS" "$E/acprior-agree.R" summary > "$A/agree.txt" 2>&1
"$RS" "$E/acprior-agree-up.R" summary > "$A/agree-up.txt" 2>&1
"$RS" "$E/acprior-innov.R" summary > "$A/innov.txt" 2>&1
echo "== finished $(date +%T)"
for f in "$M"/*.txt "$A"/*.txt "$L"/*.txt; do
  echo "-- $f: $(grep -E 'lib:|libPaths|frmtmb from' "$f" | head -1)"
done
echo "== against the lanes' own logs"
MO="$LW-mocode/dev/mocode-out"; AO="$LW-acprior/dev/acprior-out"
for arm in base lane; do
  for s in agree se3 edge; do
    [ -f "$MO/$s-$arm.txt" ] && same "$M/$s-$arm.txt" "$MO/$s-$arm.txt"
  done
done
for arm in before after; do
  for s in ml frm-table roundtrip; do
    [ -f "$AO/$s-$arm.txt" ] && same "$A/$s-$arm.txt" "$AO/$s-$arm.txt"
  done
done
for s in agree agree-up innov; do
  [ -f "$AO/$s.txt" ] && same "$A/$s.txt" "$AO/$s.txt"
done
for s in loglik-brms-agree loglik-check; do
  [ -f "$LW-loglik/dev/$s.txt" ] && same "$L/$s.txt" "$LW-loglik/dev/$s.txt"
done
echo "EQUALITY DONE"
