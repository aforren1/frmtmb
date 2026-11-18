#!/usr/bin/env bash
# Lane drawdir: summarise dev/drawdir-v-<group>.tsv (dev/drawdir-verdict.R)
# by design: fits with a verdict, fired by the summed check (dir is NA),
# fired by the per-direction reading, the largest dir of an unfired fit,
# and each fit the reading fires on. No R process.
# Usage: bash dev/drawdir-v-sum.sh <group> [...]
cd "$(dirname "$0")/.."
for g in "$@"; do
  awk -F'\t' -v grp="$g" '
    NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    {
      id = $h["id"]; d = id; sub(/:[0-9]+$/, "", d)
      if (grp == "gglm" || grp == "gglmm" || grp == "wb" || grp == "kd") {
        n = split(d, a, ":"); d = a[1]; if (grp == "kd") d = a[1] ":" a[2]
      }
      if (!(d in fits)) order[++no] = d
      fits[d]++
      if ($h["verdict"] != "TRUE") { nov[d]++; next }
      dir = $h["dir"]; fire = $h["fire"]
      if (fire == "TRUE" && dir == "NA") summed[d]++
      else if (fire == "TRUE") { bydir[d]++; list[d] = list[d] " " id "(" sprintf("%.3g", dir) ")" }
      else if (dir != "NA" && (!(d in mx) || dir + 0 > mx[d])) mx[d] = dir + 0
    }
    END {
      for (k = 1; k <= no; k++) {
        d = order[k]
        printf "%-36s fits %3d noverdict %2d | summed %2d | dir %2d | max dir unfired %s%s\n",
          d, fits[d], nov[d], summed[d], bydir[d],
          (d in mx) ? sprintf("%.3g", mx[d]) : "-", list[d]
      }
    }' "dev/drawdir-v-$g.tsv"
done
