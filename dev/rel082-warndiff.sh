#!/usr/bin/env bash
# Files whose pass, fail, err or warn count moved between two per-file runs
# of dev/ciharden-scan.sh (dev/ciharden-log/<run>.sum), for the reruns
# after lane boundfam's OpenBLAS fix, whose new warning can reach any
# fitted value that reads a lost parameter with a NaN variance.
#
#   bash dev/rel082-warndiff.sh <run before> <run after>
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel082 || exit 1
key() { awk '/RESULT/ { k = $1 " " $3; v = ""; for (i = 4; i <= NF; i++) v = v " " $i; print k "\t" v }' "dev/ciharden-log/$1.sum" | sort; }
join -t $'\t' -a 1 -a 2 -e MISSING -o 0,1.2,2.2 <(key "$1") <(key "$2") |
  awk -F'\t' '$2 != $3 { print "MOVED", $1, "|", $2, "->", $3; n++ } END { print "files moved:", n + 0 }'
