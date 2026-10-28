#!/usr/bin/env bash
# The priordraws reviewer's checks the consolidation was asked to rerun
# on the merge (dev/reviews/2026-10-27-priordraws.md, findings 1 and 2):
# dev/review-priordraws-acdefault.R (default ar, ma and cosy under
# sample_prior) and O3 of dev/review-priordraws-only.R (frm arm, set 2)
# against the brms arm's saved draws, compared by
# dev/review-priordraws-compare.R. Each runs from a copy in
# dev/rel077-log/eq/priordraws-review/ with the lane libraries, the
# reviewer's Stan cache and log folder replaced by rellib-r15, this
# tree's cache and that folder (the lane tree is read only). brms's arm
# does not read frmtmb; its rds files are copied from the reviewer's log.
#   bash dev/rel077-priordraws-review.sh > dev/rel077-log/priordraws-review.txt 2>&1
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel077
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LT=C:/Users/adf44/source/r/frmtmb-wt-priordraws
D=$ROOT/dev/rel077-log/eq/priordraws-review
rm -rf "$D"; mkdir -p "$D"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
for s in review-priordraws-acdefault.R review-priordraws-only.R \
         review-priordraws-compare.R; do
  sed -e "s#C:/Users/adf44/source/r/wt-priordraws-rv-lib#C:/Users/adf44/source/r/rellib-r15#g" \
      -e "s#C:/Users/adf44/source/r/wt-priordraws-lib#C:/Users/adf44/source/r/rellib-r15#g" \
      -e "s#$LT/dev/review-priordraws-stancache#$ROOT/dev/stan-cache#g" \
      -e "s#$LT/dev/review-priordraws-log#$D#g" \
      "$LT/dev/$s" > "$D/$s"
  diff "$LT/dev/$s" "$D/$s" | grep '^>' | sed "s/^/   $s: /"
done
cp "$LT"/dev/review-priordraws-log/only-brms-O3.rds \
   "$LT"/dev/review-priordraws-log/only-brms-O4.rds "$D/"
"$RS" "$D/review-priordraws-acdefault.R" rv > "$D/acdefault-rv.txt" 2>&1 &
"$RS" "$D/review-priordraws-only.R" frm 2 > "$D/only-frm-2.txt" 2>&1 &
wait
"$RS" "$D/review-priordraws-compare.R" > "$D/compare.txt" 2>&1
echo "== acdefault (rellib-r15)"; cat "$D/acdefault-rv.txt"
echo "== only frm 2 (rellib-r15)"; cat "$D/only-frm-2.txt"
echo "== compare: the only rows"; sed -n '/#### only/,$p' "$D/compare.txt"
echo "== the reviewer's own O3 rows (lane alone)"
sed -n '/#### only/,$p' "$LT/dev/review-priordraws-log/compare.txt" 2>/dev/null | grep -E "O3|O4|only" | head -20
echo "PRIORDRAWS REVIEW DONE"
