#!/bin/bash
# Per-package totals of a run-par.sh tier, from its per-file logs, and a
# check that every log's lib line names the lane library for frmtmb.
# Usage: bash dev/acprior-sum.sh <tier>
D=dev/release/$1-files
for p in frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spatial frmtmb.survival; do
  ls $D/$p--*.txt >/dev/null 2>&1 || continue
  cat $D/$p--*.txt | grep -a "^RESULT" | awk -v P=$p '
    { for (i = 1; i <= NF; i++) { split($i, kv, "="); if (kv[1] in s || kv[1] ~ /^(pass|fail|err|skip|warn)$/) s[kv[1]] += kv[2] } n++ }
    END { printf "%-16s files=%d pass=%d fail=%d err=%d skip=%d warn=%d\n", P, n, s["pass"], s["fail"], s["err"], s["skip"], s["warn"] }'
done
echo "logs whose frmtmb is not the lane lib: $(grep -aL 'frmtmb 0.75.0 C:/Users/adf44/source/r/wt-acprior-lib' $D/*.txt | wc -l)"
echo "logs without a RESULT line: $(grep -aL '^RESULT' $D/*.txt | wc -l)"
