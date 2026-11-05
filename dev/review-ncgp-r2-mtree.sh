#!/usr/bin/env bash
# Fix round 1: a trial merge with gpgap that builds. A copy of this tree, gpgap's hunks by diff3 over the base, each
# conflict taken from this tree. Installed only into rev-ncgp-lib. Reads the two worktrees; changes neither.
export PATH="/c/rtools45/usr/bin:$PATH"
cd /c/Users/adf44/source/r
M=rev-ncgp-src/mtree; rm -rf $M; mkdir -p $M
for f in DESCRIPTION NAMESPACE NEWS.md R man src tests vignettes inst .Rbuildignore; do [ -e frmtmb-wt-ncgp/$f ] && cp -r frmtmb-wt-ncgp/$f $M/; done
base=frmtmb-wt-r16base; a=frmtmb-wt-ncgp; b=frmtmb-wt-gpgap; T=rev-ncgp-src/t3
(cd $b; git status --short | grep -v "dev/\|\.github\|_pkgdown\|extensions/\|^?? tests" | cut -c4-) > $T.list
c=0
while read -r f; do
  [ -f "$a/$f" ] || continue
  tr -d '\r' < $a/$f > $T.a; tr -d '\r' < $base/$f > $T.o; tr -d '\r' < $b/$f > $T.b
  diff3 -m -L ncgp -L base -L gpgap $T.a $T.o $T.b > $T.m
  n=$(grep -c '^<<<<<<< ncgp' $T.m); c=$((c + n))
  # keep this tree's side of a conflict
  # a "<<<<<<< base" region is one both changed the same way: keep the shared side
  awk '/^<<<<<<< ncgp/{s=1;next} /^<<<<<<< base/{s=4;next} /^\|\|\|\|\|\|\| base/{if(s==1){s=2;next}} /^=======$/{if(s==2){s=3;next}; if(s==4){s=5;next}} /^>>>>>>> gpgap/{if(s==3||s==5){s=0;next}} s==0||s==1||s==5' $T.m > $M/$f
  echo "$f conflicts $n"
done < $T.list
echo "conflicts taken from this tree: $c"
rm -f $T.a $T.o $T.b $T.m $T.list
grep -c "gp_nugget_message" $M/R/fit.R
grep -n "gp_block_nugget\|gp_sampler_nugget" $M/R/*.R | cut -c1-100 | head
mkdir -p rev-ncgp-lib
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-ncgp-lib --no-multiarch --no-docs $M 2>&1 | tail -3
