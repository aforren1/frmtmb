#!/usr/bin/env bash
# Compare every arm of dev/samplegen-brmsout.R against the brms-only
# control, call by call, on the value hash and the warnings hash.
cd /c/Users/adf44/source/r/frmtmb-wt-samplegen/dev/samplegen-out/brmsout || exit 1
rows() { awk '$2=="value"||$2=="ERROR"{print $1,$2,$3,$4}' "$1"; }
echo '```'
echo "== brms output with and without frmtmb.sample, dev/samplegen-brmsout.R =="
echo "cached 400-row lognormal brmsfit, seed 20260915 before every call"
echo "control: library(brms) alone, fix-none.txt"
n=$(rows fix-none.txt | wc -l)
echo "calls in the control: $n ($(rows fix-none.txt | awk '$2=="ERROR"' | wc -l) of them brms's own error)"
echo
printf "%-6s %-38s %s\n" build arm "calls differing from control"
for f in base-S base-U fix-S fix-T fix-U; do
  case $f in *-S) arm="library(brms); library(frmtmb.sample)";;
             *-T) arm="library(frmtmb.sample); library(brms)";;
             *-U) arm="library(frmtmb.sample); brms loaded";; esac
  k=$(diff <(rows fix-none.txt) <(rows $f.txt) | grep -c "^>")
  printf "%-6s %-38s %d of %d\n" "${f%%-*}" "$arm" "$k" "$(rows $f.txt | wc -l)"
done
echo
echo "base-S, the first three differing calls:"
grep ERROR base-S.txt | head -3 | awk '{$3="";$4="";print "  "$0}'
echo '```'
