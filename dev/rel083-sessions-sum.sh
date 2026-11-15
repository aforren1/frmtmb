#!/usr/bin/env bash
# One-session runs of the 0.83.0 consolidation (dev/release/run-sessions.sh,
# tiers rel0.83.0 plain and rel0.83.0g gated): every SESSION line, the
# lib: lines, and the LEAK lines by tag, from the per-session logs. A LEAK
# line has one of two tags: "[with load of ...]" (what a package load
# brings) or "[extension session state]" (an extension's once-per-session
# cache, which has no restore seam; dev/release/run-session.R). A line
# with neither tag is a test's leak and is listed. One table row per suite.
#   bash dev/rel083-sessions-sum.sh <tier> > dev/rel083-log/sessions-<tier>-sum.txt
cd /c/Users/adf44/source/r/frmtmb-wt-rel083 || exit 1
t="$1"; D=dev/release/$t-sessions
echo "tier $t: $(tail -1 dev/release/$t.log)"
echo "SESSION lines: $(grep -h '^SESSION' $D/*.txt | wc -l);" \
     "with fail, err or warn not 0: $(grep -h '^SESSION' $D/*.txt | grep -vc 'fail=0 err=0 skip=[0-9]* warn=0')"
echo "lib: lines: $(grep -h '^lib:' $D/*.txt | wc -l); naming rellib-r22: $(grep -h '^lib:' $D/*.txt | grep -c 'rellib-r22')"
echo "LEAK lines: $(grep -h '^LEAK' $D/*.txt | wc -l); without a tag: $(grep -h '^LEAK' $D/*.txt | grep -v -e '\[with load of' -e '\[extension session state\]' | wc -l)" \
     "; tagged as a package load: $(grep -h '^LEAK' $D/*.txt | grep -c '\[with load of')" \
     "(a Stan model load (rstan): $(grep -h '^LEAK' $D/*.txt | grep -c 'with load of a Stan model'))" \
     "; tagged as extension session state: $(grep -h '^LEAK' $D/*.txt | grep -c '\[extension session state\]')" \
     "; naming package:emmeans: $(grep -h '^LEAK' $D/*.txt | grep -c 'package:emmeans')"
grep -h '^LEAK' $D/*.txt | grep '\[extension session state\]' | sed -E 's/^LEAK [^ ]+ //; s/ \[extension session state\]$//' |
  sort | uniq -c | sed 's/^/   EXTENSION STATE /'
grep -h '^LEAK' $D/*.txt | grep -v -e '\[with load of' -e '\[extension session state\]' | sed 's/^/   UNTAGGED /'
echo "| suite | files | pass forward, reverse, preload | skip forward, reverse, preload |"
echo "|---|---|---|---|"
for p in frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.survival frmtmb.spatial; do
  g() { grep -h '^SESSION' $D/$p--$1.txt | sed -E "s/.* $2=([0-9/]+).*/\1/"; }
  echo "| $p | $(g forward files) | $(g forward pass), $(g reverse pass), $(g preload pass) | $(g forward skip), $(g reverse skip), $(g preload skip) |"
done
