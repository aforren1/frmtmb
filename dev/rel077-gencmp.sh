#!/usr/bin/env bash
# The files roxygen regenerated in the release tree (man/*.Rd and
# NAMESPACE of every package) against the lanes' own: a file one lane
# changed must equal that lane's copy; a file two or more lanes changed
# must equal the three-way merge of the lane copies over the 435e7ce6
# blob, taken lane by lane, both sides kept where they insert at one
# anchor (--union); a file no lane changed must equal the HEAD blob,
# unless a release fix changed its source (named in $FIXED). Line
# endings aside. Adapted from dev/rel076-gencmp.sh.
#   bash dev/rel077-gencmp.sh
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-rel077
cd "$REL" || exit 1
W=dev/rel077-log/merge-work
T=$(mktemp -d)
# the release fixes whose roxygen these pages carry: S1/S2 (mointerop)
# rewrite ?frmtmb-marginaleffects; the R1 refusal and as.brmsprior()'s
# code change no Rd text
FIXED="man/frmtmb-marginaleffects.Rd"
gen=$( (git ls-files 'man/*.Rd' 'NAMESPACE' 'extensions/*/man/*.Rd' \
          'extensions/*/NAMESPACE'; git ls-files --others --exclude-standard \
          'man/*.Rd' 'extensions/*/man/*.Rd') | sort -u)
same=0; diffs=0
for f in $gen; do
  lanes=""
  for l in sigmainn mointerop reslot priordraws; do
    cut -f2 "$W/$l.status" | grep -qxF "$f" && lanes="$lanes $l"
    [ -z "$(cd ../frmtmb-wt-$l && git ls-files --others --exclude-standard -- "$f")" ] ||
      lanes="$lanes $l"
  done
  tr -d '\r' < "$f" > "$T/rel"
  set -- $lanes
  case $# in
    0) git show "HEAD:$f" 2>/dev/null | tr -d '\r' > "$T/ref"; what=HEAD;;
    1) tr -d '\r' < "../frmtmb-wt-$1/$f" > "$T/ref"; what=$1;;
    *) git show "HEAD:$f" 2>/dev/null | tr -d '\r' > "$T/base"
       tr -d '\r' < "../frmtmb-wt-$1/$f" > "$T/ref"; what=$1; shift
       for l in "$@"; do
         tr -d '\r' < "../frmtmb-wt-$l/$f" > "$T/b"
         git merge-file -q --union "$T/ref" "$T/base" "$T/b"
         what="$what+$l"
       done;;
  esac
  if cmp -s "$T/rel" "$T/ref"; then
    same=$((same + 1))
    [ -n "$lanes" ] && echo "equal to $what: $f"
  else
    diffs=$((diffs + 1))
    echo "DIFFERS from $what: $f$(echo " $FIXED " | grep -q " $f " && echo ' (release fix)')"
    diff "$T/ref" "$T/rel" | head -20 | sed 's/^/      /'
  fi
done
echo "generated files: $(echo "$gen" | wc -w), equal: $same, differ: $diffs"
rm -rf "$T"
