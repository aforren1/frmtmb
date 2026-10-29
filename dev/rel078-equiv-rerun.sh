#!/usr/bin/env bash
# Release fix R4: the equivalence survey's pairs that the two code lanes
# reach, remeasured on the release build with the copied instrument
# (dev/equiv-measure.R, one pair per R process), into a scratch outdir,
# WITHOUT --record: dev/equiv-data/ and dev/equiv-findings.md are not
# written (equiv-measure.R writes only <outdir>/<id>.{rds,tsv}; the
# summary and the sweep are not run). The rescor pairs (rescorpp) and
# the me() pairs (meme) are the ones dev/equiv-findings.md, section 6,
# names. dev/rel078-equiv-cmp.R compares each with the lane's own
# measurement (../frmtmb-wt-equiv/dev/equiv-out/<id>.tsv, read only) and
# with the record (dev/equiv-data/results.tsv).
#
#   bash dev/rel078-equiv-rerun.sh [N] > dev/rel078-log/equiv-rerun.txt 2>&1
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/12 TEMP=C:/Users/adf44/AppData/Local/Temp/12
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel078
cd "$ROOT" || exit 1
N="${1:-4}"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export FRMTMB_EQUIV_LIB=C:/Users/adf44/source/r/rellib-r16
export NOT_CRAN=true
out=dev/rel078-log/equiv-r4
rm -rf "$out"; mkdir -p "$out/logs"
ids="p002 p030 p041 p048 p068 p151 p153 p159 p028 p088 p095 p115 p125"
md5sum dev/equiv-data/* dev/equiv-findings.md > dev/rel078-log/equiv-record-md5-before.txt
start=$(date +%s)
printf '%s\n' $ids | xargs -P "$N" -I{} bash -c \
  "timeout 2400 \"$R\" dev/equiv-measure.R {} \"$out\" \
     > \"$out/logs/{}.log\" 2>&1 ||
   echo \"FAILED {} exit \$?\" >> \"$out/logs/{}.log\""
echo "measured $(echo $ids | wc -w) pairs in $(( $(date +%s) - start )) s"
md5sum dev/equiv-data/* dev/equiv-findings.md > dev/rel078-log/equiv-record-md5-after.txt
if cmp -s dev/rel078-log/equiv-record-md5-before.txt dev/rel078-log/equiv-record-md5-after.txt; then
  echo "dev/equiv-data/ and dev/equiv-findings.md unchanged"
else
  echo "THE RECORD CHANGED"
fi
grep -h "^lib:\|^verdict\|FAILED" "$out"/logs/*.log
"$R" dev/rel078-equiv-cmp.R "$out"
