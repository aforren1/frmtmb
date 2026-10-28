#!/usr/bin/env bash
# Each lane's equality and agreement scripts, rerun on the merged build.
# Adapted from dev/rel076-equality.sh, with its lesson
# (dev/round-handoff.md, 0.76.0: a lane's log can be older than its last
# fix round): every script runs from two copies, one with the lane
# library replaced by rellib-r15 (tag "m", the merged build) and one
# with the lane's own final library, read only (tag "l"); outputs go to
# dev/rel077-log/eq/<lane>/<tag>/ instead of the lane's tree, and the
# copies read this tree's dev/ files. The two are compared line for
# line past the lines that name a library, and each is also compared
# with the lane's own log. A difference between "m" and "l" is the
# merge's; it is explained in dev/round-20261028.md or it is a finding.
# The "before" arms read rellib-r14 (0.76.0) in both copies. Sampling
# scripts against brms run the frmtmb arm only; brms's draws are copied
# from the lane's output.
#
#   bash dev/rel077-equality.sh > dev/rel077-log/equality.txt 2>&1
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel077
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R15=C:/Users/adf44/source/r/rellib-r15
LW=C:/Users/adf44/source/r/frmtmb-wt
E=$ROOT/dev/rel077-log/eq
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export NOT_CRAN=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
rm -rf "$E"
for l in sigmainn mointerop priordraws reslot; do
  mkdir -p "$E/$l/m" "$E/$l/l"
done
# copy <lane> <script>: both copies, the lane's paths replaced
copy() {
  local l="$1" s="$2"
  for t in m l; do
    local d="$E/$l/$t"
    sed -e "s#file.path(wt, \"dev\", \"$l-out\", \"route-before.rds\")#\"$d/route-before.rds\"#g" \
        -e "s#file.path(wt, \"dev\", \"$l-out\")#\"$d\"#g" \
        -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l/dev/$l-log\"#$d\"#g" \
        -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l/dev/$l-out\"#$d\"#g" \
        -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l/dev/reslot-common.R#$d/reslot-common.R#g" \
        -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l/dev/stan-cache#$ROOT/dev/stan-cache#g" \
        -e "s#\"dev/priordraws-log/brms-names-%d.txt\"#\"$d/brms-names-%d.txt\"#g" \
        -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l\"#$ROOT\"#g" \
        "$LW-$l/dev/$s" > "$d/$s"
    if [ "$t" = m ]; then
      sed -i "s#C:/Users/adf44/source/r/wt-$l-lib#$R15#g" "$d/$s"
    fi
  done
  diff "$LW-$l/dev/$s" "$E/$l/m/$s" | grep '^>' | sed "s/^/   $l m $s: /"
}
for s in sigmainn-unchanged.R sigmainn-readers.R sigmainn-agree.R \
         sigmainn-same-point.R; do copy sigmainn "$s"; done
for s in mointerop-equal.R mointerop-mo12.R mointerop-sample.R; do
  copy mointerop "$s"; done
for s in priordraws-identity.R priordraws-names-frm.R; do
  copy priordraws "$s"; done
for s in reslot-common.R reslot-identical.R reslot-roundtrip.R \
         reslot-double.R reslot-flatpost.R; do copy reslot "$s"; done
for t in m l; do
  cp "$LW-sigmainn/dev/sigmainn-out/route-before.rds" \
     "$LW-sigmainn/dev/sigmainn-out/agree-before.rds" \
     "$LW-sigmainn/dev/sigmainn-out/agree-brms.rds" "$E/sigmainn/$t/"
  cp "$LW-priordraws"/dev/priordraws-log/brms-names-[123].txt "$E/priordraws/$t/"
done
echo "== started $(date +%T)"
for t in m l; do
  S=$E/sigmainn/$t; M=$E/mointerop/$t; P=$E/priordraws/$t; Q=$E/reslot/$t
  "$RS" "$S/sigmainn-unchanged.R" after > "$S/unchanged-after.txt" 2>&1 &
  "$RS" "$S/sigmainn-readers.R" after > "$S/readers-after.txt" 2>&1 &
  wait
  "$RS" "$S/sigmainn-same-point.R" > "$S/same-point.txt" 2>&1 &
  "$RS" "$M/mointerop-equal.R" lane > "$M/equal-lane.txt" 2>&1 &
  wait
  "$RS" "$M/mointerop-mo12.R" lane > "$M/mo12-lane.txt" 2>&1 &
  "$RS" "$P/priordraws-names-frm.R" \
    "$([ $t = m ] && echo $R15 || echo C:/Users/adf44/source/r/wt-priordraws-lib)" \
    > "$P/names-frm.txt" 2>&1 &
  wait
  "$RS" "$Q/reslot-identical.R" after > "$Q/identical-after.txt" 2>&1 &
  "$RS" "$Q/reslot-roundtrip.R" after > "$Q/roundtrip-after.txt" 2>&1 &
  wait
  "$RS" "$Q/reslot-double.R" after > "$Q/double-after.txt" 2>&1 &
  wait
  for r in 1 2 3; do
    "$RS" "$P/priordraws-identity.R" lane "$r" > "$P/identity-run-$r.log" 2>&1 &
  done
  wait
done
# the before arms (rellib-r14), from the "m" copies
S=$E/sigmainn/m; M=$E/mointerop/m; P=$E/priordraws/m; Q=$E/reslot/m
"$RS" "$S/sigmainn-unchanged.R" before > "$S/unchanged-before.txt" 2>&1 &
"$RS" "$S/sigmainn-readers.R" before > "$S/readers-before.txt" 2>&1 &
wait
"$RS" "$M/mointerop-equal.R" base > "$M/equal-base.txt" 2>&1 &
"$RS" "$Q/reslot-identical.R" before > "$Q/identical-before.txt" 2>&1 &
for r in 1 2 3; do
  "$RS" "$P/priordraws-identity.R" base "$r" > "$P/identity-runb-$r.log" 2>&1 &
done
wait
"$RS" "$Q/reslot-identical.R" compare > "$Q/identical-compare.txt" 2>&1
cp "$Q/identical-before.rds" "$E/reslot/l/" 2>/dev/null
"$RS" "$E/reslot/l/reslot-identical.R" compare > "$E/reslot/l/identical-compare.txt" 2>&1
# the sampling comparisons against brms, the merged build only
"$RS" "$S/sigmainn-agree.R" after > "$S/agree-after.txt" 2>&1 &
"$RS" "$M/mointerop-sample.R" lane > "$M/sample-lane.txt" 2>&1 &
"$RS" "$Q/reslot-flatpost.R" after > "$Q/flatpost-after.txt" 2>&1 &
wait
"$RS" "$S/sigmainn-agree.R" summary > "$S/agree-summary.txt" 2>&1
echo "== finished $(date +%T)"
# <a> <b>: line-for-line comparison past library lines
same() {
  a=$(grep -vE '^lib|libPaths|rellib|wt-[a-z]*-lib|^\[1\] "C:|^arm:' "$1" | tr -d '\r')
  b=$(grep -vE '^lib|libPaths|rellib|wt-[a-z]*-lib|^\[1\] "C:|^arm:' "$2" | tr -d '\r')
  if [ "$a" = "$b" ]; then
    echo "SAME    $1 = $2 ($(echo "$a" | wc -l) lines)"
  else
    echo "DIFFERS $1 against $2:"
    echo "$b" > "$E/same-b.tmp"; echo "$a" > "$E/same-a.tmp"
    diff "$E/same-b.tmp" "$E/same-a.tmp" | head -40 | sed "s/^/      /"
  fi
}
echo "== merged (m) against the lane's own library (l)"
for f in sigmainn/unchanged-after.txt sigmainn/readers-after.txt \
         sigmainn/same-point.txt mointerop/equal-lane.txt \
         mointerop/mo12-lane.txt priordraws/names-frm.txt \
         reslot/identical-after.txt reslot/roundtrip-after.txt \
         reslot/double-after.txt reslot/identical-compare.txt; do
  same "$E/${f%%/*}/m/${f#*/}" "$E/${f%%/*}/l/${f#*/}"
done
for r in 1 2 3; do
  same "$E/priordraws/m/identity-lane-$r.txt" "$E/priordraws/l/identity-lane-$r.txt"
done
echo "== the identity runs: distinct md5 lines per arm"
for t in m l; do
  for a in lane base; do
    ls "$E/priordraws/$t/identity-$a-"*.txt > /dev/null 2>&1 || continue
    echo "  $t $a: $(cat "$E/priordraws/$t/identity-$a-"*.txt | tr -d '\r' | grep -v '^lib' | sort -u | wc -l) distinct lines over $(ls "$E/priordraws/$t/identity-$a-"*.txt | wc -l) runs"
  done
done
echo "== against the lanes' own logs"
SO="$LW-sigmainn/dev/sigmainn-out"; MO="$LW-mointerop/dev/mointerop-out"
PO="$LW-priordraws/dev/priordraws-log"; QO="$LW-reslot/dev/reslot-out"
for s in unchanged-after unchanged-before readers-after readers-before \
         same-point agree-summary; do
  same "$E/sigmainn/m/$s.txt" "$SO/$s.txt"
done
for s in equal-lane equal-base mo12-lane sample-lane; do
  same "$E/mointerop/m/$s.txt" "$MO/$s.txt"
done
same "$E/priordraws/m/names-frm.txt" "$PO/names-frm.txt"
for s in identical-after identical-before identical-compare roundtrip-after \
         double-after flatpost-after; do
  same "$E/reslot/m/$s.txt" "$QO/$s.txt"
done
echo "EQUALITY DONE"
