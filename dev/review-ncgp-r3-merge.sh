#!/usr/bin/env bash
# Short check: trial merge of ncgp, gpgap by the whole-file list, and refuse by diff3 with its mi_retape() scoped. Writes only under
# rev-ncgp-src/mtree3 and rev-ncgp-lib. No git state changes.
export PATH="/c/rtools45/usr/bin:$PATH"
cd /c/Users/adf44/source/r
M=rev-ncgp-src/mtree3; rm -rf $M; mkdir -p $M
for f in DESCRIPTION NAMESPACE NEWS.md R man src tests vignettes inst .Rbuildignore; do [ -e frmtmb-wt-ncgp/$f ] && cp -r frmtmb-wt-ncgp/$f $M/; done
base=frmtmb-wt-r16base; T=rev-ncgp-src/t4
take() { # lane file: three-way into the tree, conflicts keep the tree's side
  b=frmtmb-wt-$1; f=$2
  mkdir -p "$M/$(dirname $f)"
  if [ ! -f "$base/$f" ]; then [ -f "$M/$f" ] || cp "$b/$f" "$M/$f"; echo "$1 $f new"; return; fi
  [ -f "$M/$f" ] || cp "$base/$f" "$M/$f"
  tr -d '\r' < $M/$f > $T.a; tr -d '\r' < $base/$f > $T.o; tr -d '\r' < $b/$f > $T.b
  diff3 -m -L ncgp -L base -L lane $T.a $T.o $T.b > $T.m
  n=$(grep -c '^<<<<<<< ncgp' $T.m)
  awk '/^<<<<<<< ncgp/{s=1;next} /^<<<<<<< base/{s=4;next} /^\|\|\|\|\|\|\| base/{if(s==1){s=2;next}} /^=======$/{if(s==2){s=3;next}; if(s==4){s=5;next}} /^>>>>>>> lane/{if(s==3||s==5){s=0;next}} s==0||s==1||s==5' $T.m > $M/$f
  [ "$n" != 0 ] && echo "$1 $f conflicts $n kept this tree's side"
}
# gpgap: these stay this tree's whole; the rest merges
keep="R/covstruct.R R/fit.R R/predict.R R/conditions.R R/sampling-api.R man/frmtmb-conditions.Rd man/frmtmb-sampling-api.Rd NEWS.md tests/testthat/test-gp-nugget-floor.R"
(cd frmtmb-wt-gpgap; git status --short | grep -v "dev/\|\.github\|_pkgdown\|extensions/" | cut -c4-) > $T.list
while read -r f; do case " $keep " in *" $f "*) echo "gpgap $f whole from ncgp";; *) take gpgap "$f";; esac; done < $T.list
# refuse: everything of core
(cd frmtmb-wt-refuse; git status --short | grep -v "dev/\|\.github\|_pkgdown\|extensions/" | cut -c4-) > $T.list
while read -r f; do [ "$f" = NEWS.md ] && continue; take refuse "$f"; done < $T.list
rm -f $T.a $T.o $T.b $T.m $T.list
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" /c/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-r3-mipatch.R
mkdir -p rev-ncgp-lib
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-ncgp-lib --no-multiarch --no-docs $M 2>&1 | tail -2
grep -c "gp_block_nugget <-\|gp_sampler_nugget" $M/R/covstruct.R
