#!/usr/bin/env bash
# Trial merges of this lane with each other lane, by diff3 on file copies, outside
# every worktree. No git state is touched: base files come from frmtmb-wt-r16base.
R=/c/Users/adf44/source/r
W=$R/frmtmb-wt-nlfix
OUT=$R/rev-nlfix-src/merge
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
for t in latres refuse predread; do
  S=$OUT/src-$t
  rm -rf "$S"; mkdir -p "$S"
  cp -r $W/R $W/man $W/inst $W/data $W/DESCRIPTION $W/NAMESPACE "$S/"
  cd $R/frmtmb-wt-$t
  # every changed or new path of the other lane under R/, NAMESPACE, man/
  files=$(git status --short R NAMESPACE man inst data 2>/dev/null | awk '{print $NF}')
  for f in $files; do
    if [ -d "$f" ]; then continue; fi
    if [ -f "$R/frmtmb-wt-nlfix/$f" ] && [ -f "$R/frmtmb-wt-r16base/$f" ] && ! cmp -s "$R/frmtmb-wt-nlfix/$f" "$R/frmtmb-wt-r16base/$f"; then
      diff3 -m --strip-trailing-cr "$R/frmtmb-wt-nlfix/$f" "$R/frmtmb-wt-r16base/$f" "$R/frmtmb-wt-$t/$f" > "$S/$f.merged" 2>/dev/null
      n=$(grep -c '^<<<<<<<' "$S/$f.merged")
      echo "$t $f conflicts=$n"
      if [ "$n" != "0" ]; then
        # keep this lane's side of a conflict; the review reads the conflict separately
        awk '/^<<<<<<< /{m=1; next} /^\|\|\|\|\|\|\| /{m=2; next} /^=======$/{ if (m==2) {m=3; next} } /^>>>>>>> /{m=0; next} { if (m==0 || m==1) print }' "$S/$f.merged" > "$S/$f"
        rm "$S/$f.merged"
      else
        mv "$S/$f.merged" "$S/$f"
      fi
    else
      mkdir -p "$S/$(dirname $f)"
      cp "$f" "$S/$f"
    fi
  done
  L=$R/rev-nlfix-lib/merge-$t
  mkdir -p "$L"
  cd "$S/.."
  mv "src-$t" frmtmb
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library="$L" --no-multiarch --no-docs --no-test-load frmtmb > "$W/dev/review-nlfix-merge-$t.log" 2>&1
  tail -1 "$W/dev/review-nlfix-merge-$t.log" | cut -c1-160
  mv frmtmb "src-$t"
done
