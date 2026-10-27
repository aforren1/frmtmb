#!/usr/bin/env bash
# The files roxygen regenerated in the release tree (man/*.Rd and
# NAMESPACE of every package) against the lanes' own: a file one lane
# changed must equal that lane's copy; a file two lanes changed must
# equal the three-way merge of the two lane copies over the 70476dc3
# blob, both sides kept where they insert at one anchor (--union); a
# file no lane changed must equal the HEAD blob, unless a
# release fix changed its source (named in $FIXED). Line endings aside.
#   bash dev/rel076-gencmp.sh
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-rel076
cd "$REL" || exit 1
W=dev/rel076-log/merge-work
T=$(mktemp -d)
# the release fix to R/priors.R rewrites ?set_prior's Residual
# correlation section, which acprior also changed
FIXED="man/set_prior.Rd extensions/frmtmb.sample/man/frm_sample.Rd"
gen=$( (git ls-files 'man/*.Rd' 'NAMESPACE' 'extensions/*/man/*.Rd' \
          'extensions/*/NAMESPACE'; git ls-files --others --exclude-standard \
          'man/*.Rd' 'extensions/*/man/*.Rd') | sort -u)
same=0; diffs=0
for f in $gen; do
  lanes=""
  for l in mocode acprior loglik; do
    cut -f2 "$W/$l.status" | grep -qxF "$f" && lanes="$lanes $l"
  done
  tr -d '\r' < "$f" > "$T/rel"
  set -- $lanes
  case $# in
    0) git show "HEAD:$f" 2>/dev/null | tr -d '\r' > "$T/ref"; what=HEAD;;
    1) tr -d '\r' < "../frmtmb-wt-$1/$f" > "$T/ref"; what=$1;;
    *) git show "HEAD:$f" | tr -d '\r' > "$T/base"
       tr -d '\r' < "../frmtmb-wt-$1/$f" > "$T/ref"
       tr -d '\r' < "../frmtmb-wt-$2/$f" > "$T/b"
       git merge-file -q --union "$T/ref" "$T/base" "$T/b"
       what="$1+$2";;
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
