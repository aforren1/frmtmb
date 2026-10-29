#!/usr/bin/env bash
# One-session runs of the 0.78.0 consolidation (dev/release/run-sessions.sh,
# tiers rel0.78.0 plain and rel0.78.0g gated): every SESSION line, the
# lib: lines, and the LEAK lines with and without a "[with load of ...]"
# tag, from the per-session logs. One table row per suite.
#   bash dev/rel078-sessions-sum.sh <tier> > dev/rel078-log/sessions-<tier>-sum.txt
cd /c/Users/adf44/source/r/frmtmb-wt-rel078 || exit 1
t="$1"; D=dev/release/$t-sessions
echo "tier $t: $(tail -1 dev/release/$t.log)"
echo "SESSION lines: $(grep -h '^SESSION' $D/*.txt | wc -l);" \
     "with fail, err or warn not 0: $(grep -h '^SESSION' $D/*.txt | grep -vc 'fail=0 err=0 skip=[0-9]* warn=0')"
echo "lib: lines: $(grep -h '^lib:' $D/*.txt | wc -l); naming rellib-r16: $(grep -h '^lib:' $D/*.txt | grep -c 'rellib-r16')"
echo "LEAK lines: $(grep -h '^LEAK' $D/*.txt | wc -l); without a tag: $(grep -h '^LEAK' $D/*.txt | grep -vc '\[with load of')" \
     "; tagged as a Stan model load (rstan): $(grep -h '^LEAK' $D/*.txt | grep -c 'with load of a Stan model')" \
     "; naming package:emmeans: $(grep -h '^LEAK' $D/*.txt | grep -c 'package:emmeans')"
grep -h '^LEAK' $D/*.txt | grep -v '\[with load of' | sed 's/^/   UNTAGGED /'
echo "| suite | files | pass forward, reverse, preload | skip forward, reverse, preload |"
echo "|---|---|---|---|"
for p in frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.survival frmtmb.spatial; do
  g() { grep -h '^SESSION' $D/$p--$1.txt | sed -E "s/.* $2=([0-9/]+).*/\1/"; }
  echo "| $p | $(g forward files) | $(g forward pass), $(g reverse pass), $(g preload pass) | $(g forward skip), $(g reverse skip), $(g preload skip) |"
done
