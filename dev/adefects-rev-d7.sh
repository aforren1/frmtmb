#!/bin/sh
D="C:/Users/adf44/source/r/frmtmb-wt-adefects/dev"
RS="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
ARM=${REV_ARM:-lane}
OUT="C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-log/d7-$ARM.txt"
: > "$OUT"
while IFS= read -r s; do
  [ -z "$s" ] && continue
  REV_ARM="$ARM" REV_SCEN="$s" "$RS" --vanilla "$D/adefects-rev-d7-probe.R" \
    >> "$OUT" 2>&1
  echo "P|$s|EXIT=$?" >> "$OUT"
done <<'EOF'
frmtmb
frmtmb|frmtmb.sample
frmtmb.sample
frmtmb.sample|frmtmb
brms|frmtmb
frmtmb|brms
brms|frmtmb|frmtmb.sample
frmtmb|frmtmb.sample|brms
frmtmb.sample|brms
brms|frmtmb.sample
ns:brms|frmtmb
ns:brms|frmtmb|frmtmb.sample
frmtmb|ns:brms
rstantools|frmtmb
frmtmb|rstantools
ns:rstantools|frmtmb|frmtmb.sample
frmtmb|frmtmb.sample|ns:rstantools
loo|frmtmb|frmtmb.sample
frmtmb|frmtmb.sample|loo
posterior|frmtmb|frmtmb.sample
frmtmb|frmtmb.sample|posterior
loo|posterior|brms|frmtmb|frmtmb.sample
frmtmb|frmtmb.sample|loo|posterior|brms
frmtmb|frmtmb.sample|ns:rstantools|un:rstantools
ns:rstantools|frmtmb|frmtmb.sample|un:rstantools
brms|frmtmb|frmtmb.sample|un:brms|un:rstantools
EOF
grep -c "DONE" "$OUT"
