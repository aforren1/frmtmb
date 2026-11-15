#!/usr/bin/env bash
# Forward one-session runs of core and every extension on the lane
# library, three at a time, each with its own R_PARALLEL_PORT (the fix
# round's rerun; dev/release/run-session.R does the work).
#   bash dev/hygiene-sessions-fwd.sh <tag>
tag="$1"
cd "$(dirname "$0")/.."
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-hygiene-lib NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/hygiene-out/sess-$tag
mkdir -p "$out"
port=11900
for p in frmtmb frmtmb.sample frmtmb.eam frmtmb.coupling frmtmb.latent \
         frmtmb.learn frmtmb.ode frmtmb.spatial frmtmb.survival; do
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 5; done
  port=$((port + 7))
  ( R_PARALLEL_PORT=$port "$RS" dev/release/run-session.R "$p" forward \
      > "$out/$p.txt" 2>&1 ) &
done
wait
for p in frmtmb frmtmb.sample frmtmb.eam frmtmb.coupling frmtmb.latent \
         frmtmb.learn frmtmb.ode frmtmb.spatial frmtmb.survival; do
  grep -h "^SESSION" "$out/$p.txt" || echo "NO SESSION LINE $p"
done
echo "untagged LEAK lines: $(grep -h '^LEAK' "$out"/*.txt |
  grep -vc '\[with load of\|\[extension session state\]')"
grep -h "^DETAIL ERROR\|^DETAIL FAIL\|^DETAIL WARN" "$out"/*.txt
