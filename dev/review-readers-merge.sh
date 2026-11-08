#!/usr/bin/env bash
# Review of lane readers: which top-level functions of each shared file each
# lane's diff touches, mapped through the base file 3cd2a0c8 (old-side line
# numbers of every hunk). Read-only git: diff and show only.
export GIT_OPTIONAL_LOCKS=0
R=C:/Users/adf44/source/r
for f in R/predict.R R/compat.R R/me.R R/conditional-effects.R \
         R/autocor-latent.R R/interop.R NEWS.md; do
  base=$(mktemp)
  git -C $R/frmtmb-wt-readers show 3cd2a0c8:$f > "$base" 2>/dev/null
  for L in readers latperf mfxse osa; do
    W=$R/frmtmb-wt-$L
    git -C $W diff -U0 -- $f 2>/dev/null | grep '^@@' |
      sed -E 's/^@@ -([0-9]+)(,([0-9]+))? .*/\1 \3/' |
      while read -r s c; do
        [ -z "$c" ] && c=1
        e=$((s + (c > 0 ? c - 1 : 0)))
        awk -v s="$s" -v e="$e" '
          /^[A-Za-z_.][A-Za-z0-9_.]* *(<-|=) *function/ {
            name=$1; start=NR }
          /^r\(/ && FILENAME ~ /compat/ { }
          NR>=s && NR<=e { if (name=="") name="(top)"; print name; exit }
          NR==e { print name; exit }' "$base"
      done | sort -u | sed "s|^|$f $L |"
  done
  rm -f "$base"
done
