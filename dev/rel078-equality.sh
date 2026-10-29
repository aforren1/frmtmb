#!/usr/bin/env bash
# Each code lane's equality and agreement scripts, rerun on the merged
# build. Adapted from dev/rel077-equality.sh, with its lesson
# (dev/round-handoff.md, 0.76.0: a lane's log can be older than its last
# fix round): every script runs from two copies, one with the lane
# library replaced by rellib-r16 (tag "m", the merged build) and one
# with the lane's own final library, read only (tag "l"). Outputs go to
# dev/rel078-log/eq/<lane>/<tag>/ instead of the lane's tree; the copies
# read the lane tree's inputs (brms's saved fits and summaries, test
# helpers) read only, and this tree's Stan cache. The two are compared
# line for line past the lines that name a library, and each is also
# compared with the lane's own log. A difference between "m" and "l" is
# the merge's; it is explained in dev/round-20261029.md or it is a
# finding. brms's arms are not rerun: their saved outputs are read.
# equiv is measurement only; its overlap pairs are rerun by
# dev/rel078-equiv-rerun.sh (release fix R4).
#
#   bash dev/rel078-equality.sh > dev/rel078-log/equality.txt 2>&1
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel078
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R16=C:/Users/adf44/source/r/rellib-r16
LW=C:/Users/adf44/source/r/frmtmb-wt
E=$ROOT/dev/rel078-log/eq
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export NOT_CRAN=true
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
rm -rf "$E"
for l in meme rescorpp; do mkdir -p "$E/$l/m" "$E/$l/l"; done
# copy <lane> <script>: both copies, the lane's output paths, its cache
# and its sourced helpers replaced; the "m" copy reads rellib-r16 where
# the lane read its own library
copy() {
  local l="$1" s="$2"
  for t in m l; do
    local d="$E/$l/$t"
    sed -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l/dev/$l-out\"#$d\"#g" \
        -e "s#C:/Users/adf44/source/r/frmtmb-wt-$l/dev/stan-cache#$ROOT/dev/stan-cache#g" \
        -e "s#source(\"C:/Users/adf44/source/r/frmtmb-wt-$l/dev/rescorpp-common.R\")#source(\"$d/rescorpp-common.R\")#g" \
        "$LW-$l/dev/$s" > "$d/$s"
    if [ "$t" = m ]; then
      sed -i "s#C:/Users/adf44/source/r/wt-$l-lib#$R16#g" "$d/$s"
    fi
  done
  diff "$LW-$l/dev/$s" "$E/$l/m/$s" | grep '^>' | sed "s/^/   $l m $s: /"
}
for s in meme-guard.R meme-before.R meme-l1.R meme-agree.R meme-agree-cmp.R; do
  copy meme "$s"
done
for s in rescorpp-common.R rescorpp-audit.R rescorpp-compare.R \
         rescorpp-repro.R; do
  copy rescorpp "$s"
done
for t in m l; do
  cp "$LW-meme/dev/meme-out/agree-brms.rds" "$E/meme/$t/"
done
echo "== started $(date +%T)"
# EQ_N processes at a time (default 3), so the reruns fit beside the
# one-session runs under the round's limit of about 8 R processes
EQ_N=${EQ_N:-3}
pace() { while [ "$(jobs -rp | wc -l)" -ge "$EQ_N" ]; do sleep 1; done; }
for t in m l; do
  M=$E/meme/$t; Q=$E/rescorpp/$t
  pace; "$RS" "$M/meme-guard.R" lane > "$M/guard-lane.txt" 2>&1 &
  pace; "$RS" "$M/meme-before.R" lane > "$M/before-lane.txt" 2>&1 &
  pace; "$RS" "$M/meme-l1.R" > "$M/l1.txt" 2>&1 &
  pace; "$RS" "$M/meme-agree.R" lane > "$M/agree-lane.txt" 2>&1 &
  pace; "$RS" "$Q/rescorpp-audit.R" lane > "$Q/audit-lane.txt" 2>&1 &
  pace; "$RS" "$Q/rescorpp-repro.R" lane > "$Q/repro-lane.txt" 2>&1 &
  for w in g3 gs st g2f; do
    pace; "$RS" "$Q/rescorpp-compare.R" lane "$w" > "$Q/compare-$w.txt" 2>&1 &
  done
done
wait
for t in m l; do
  "$RS" "$E/meme/$t/meme-agree-cmp.R" > "$E/meme/$t/agree-cmp.txt" 2>&1
done
echo "== finished $(date +%T)"
# drop the lines that name a library, a path or a time
strip() {
  grep -avE '^(lib|arm|frmtmb|frmtmb.sample)[: ]|rellib-r1[56]|wt-[a-z]+-lib|frmtmb-wt-|elapsed|seconds|Elapsed|sec\b' "$1"
}
cmp_pair() {   # <label> <a> <b>
  if [ ! -s "$2" ] || [ ! -s "$3" ]; then echo "MISSING $1"; return; fi
  if diff <(strip "$2") <(strip "$3") > /dev/null; then
    echo "SAME     $1 ($(strip "$2" | wc -l) lines)"
  else
    echo "DIFFERS  $1: $(diff <(strip "$2") <(strip "$3") | grep -c '^[<>]') lines"
    diff <(strip "$2") <(strip "$3") | head -16 | sed 's/^/      /'
  fi
}
echo "== merged (m) against the lane's own final library (l)"
for f in guard-lane before-lane l1 agree-lane agree-cmp; do
  cmp_pair "meme $f" "$E/meme/m/$f.txt" "$E/meme/l/$f.txt"
done
for f in audit-lane repro-lane compare-g3 compare-gs compare-st compare-g2f; do
  cmp_pair "rescorpp $f" "$E/rescorpp/m/$f.txt" "$E/rescorpp/l/$f.txt"
done
echo "== merged (m) against the lane's own log"
for f in guard-lane before-lane l1 agree-cmp; do
  cmp_pair "meme $f vs dev/meme-out/$f.txt" "$E/meme/m/$f.txt" \
    "$LW-meme/dev/meme-out/$f.txt"
done
cmp_pair "rescorpp audit vs rescorpp-audit-lane.log" "$E/rescorpp/m/audit-lane.txt" \
  "$LW-rescorpp/dev/rescorpp-audit-lane.log"
cmp_pair "rescorpp repro vs rescorpp-repro-lane.log" "$E/rescorpp/m/repro-lane.txt" \
  "$LW-rescorpp/dev/rescorpp-repro-lane.log"
for w in g3 gs st g2f; do
  cmp_pair "rescorpp compare-$w vs rescorpp-compare-$w.log" \
    "$E/rescorpp/m/compare-$w.txt" "$LW-rescorpp/dev/rescorpp-compare-$w.log"
done
echo "EQUALITY DONE"
