#!/usr/bin/env bash
# R1's arms, at most 7 at once. Usage: bash dev/latres-r1run.sh
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
mkdir -p dev/latres-out/r1
for fam in bb beta zoib nb gamma weibull; do
  for a in "none 0" "nan 25" "box 25" "nan 15" "box 15"; do
    set -- $a
    while [ "$(jobs -rp | wc -l)" -ge 7 ]; do sleep 1; done
    "$RSCRIPT" dev/latres-r1.R $1 $2 $fam dev/latres-out/r1/$fam-$1-$2.tsv \
      > dev/latres-out/r1/$fam-$1-$2.log 2>&1 &
  done
done
wait
echo r1-done
