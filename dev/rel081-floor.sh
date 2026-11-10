#!/usr/bin/env bash
# The floors of the 0.81.0 release, measured, by dev/rel080-floor.sh's
# method: nothing is installed; each arm loads a COPY of an installed
# build (in a scratch folder outside the trees, deleted at the end) in
# front of the core it is measured on.
#   frmtmb.sample: rellib-r20's 0.28.0, floor lowered to 0.61.0, over
#     rellib-r19 (0.80.0); every test file of the suite from the release
#     tree, one per process, every gate on. Converse: rellib-r19's 0.27.0
#     over rellib-r20 (load).
#   frmtmb.survival: rellib-r20's 0.2.3 (floor 0.75.0 kept) over
#     rellib-r19; its test files from the release tree, every gate on.
# The six other extensions are unchanged and call no core function this
# release changed (grep of their R/ for row_lpdf, the aterm key helpers,
# the convergence and SE-check internals); their tests ran on rellib-r20
# in the plain and gated tiers.
#
#   bash dev/rel081-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/2}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/2}"
export GIT_OPTIONAL_LOCKS=0
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel081
SCR=C:/Users/adf44/AppData/Local/Temp/2/claude/ee0bfefdf562474c/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/floor
RBASE=C:/Users/adf44/source/r/rellib-r19
RNEW=C:/Users/adf44/source/r/rellib-r20
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
N="${1:-3}"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "== e0cce906 (0.80.0): the seam the new floor rests on"
printf "  row_lpdf() has exact = : e0cce906 %s, this tree %s\n" \
  "$(git show e0cce906:R/objective.R | tr -d '\r' | grep -c '^row_lpdf <- function.*exact')" \
  "$(tr -d '\r' < R/objective.R | grep -c '^row_lpdf <- function.*exact')"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(c(m\$Depends, m\$Imports), function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
jobs_of() {  # <pkg> -> dev/rel081-log/floor-<pkg>.jobs
  : > "dev/rel081-log/floor-$1.jobs"
  for f in extensions/$1/tests/testthat/test-*.R; do
    echo "$1 $f" >> "dev/rel081-log/floor-$1.jobs"
  done
}
p=frmtmb.sample
echo "== $p $(grep '^Version:' extensions/$p/DESCRIPTION | cut -d' ' -f2 | tr -d '\r'), floor $(grep -o 'frmtmb (>= [0-9.]*)' extensions/$p/DESCRIPTION | head -1)"
FL=$SCR/old-$p; rm -rf "$FL"; mkdir -p "$FL"
cp -r "$RNEW/$p" "$FL/"
echo "-- old: rellib-r20's build over rellib-r19 (floor lowered)"
lower_floor "$FL/$p"
out=$(load_line "$FL" "$RBASE" "$p"); echo "$out"
if echo "$out" | grep -q "LOADED"; then
  jobs_of $p
  ROUND_BASE=$RBASE /usr/bin/bash dev/ciharden-scan.sh \
    "r081-floor-$p" ref "$FL" "$N" "dev/rel081-log/floor-$p.jobs"
fi
echo "-- new: rellib-r20 as built"
load_line "$RNEW" "$RNEW" "$p"
FC=$SCR/conv-$p; rm -rf "$FC"; mkdir -p "$FC"
cp -r "$RBASE/$p" "$FC/"
echo "-- converse: $p $(grep '^Version:' "$FC/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r19 over rellib-r20"
load_line "$FC" "$RNEW" "$p"
p=frmtmb.survival
echo "== $p $(grep '^Version:' extensions/$p/DESCRIPTION | cut -d' ' -f2 | tr -d '\r'), floor $(grep -o 'frmtmb (>= [0-9.]*)' extensions/$p/DESCRIPTION | head -1)"
FS=$SCR/new-$p; rm -rf "$FS"; mkdir -p "$FS"
cp -r "$RNEW/$p" "$FS/"
echo "-- rellib-r20's build over rellib-r19, floor as built"
out=$(load_line "$FS" "$RBASE" "$p"); echo "$out"
if echo "$out" | grep -q "LOADED"; then
  jobs_of $p
  ROUND_BASE=$RBASE /usr/bin/bash dev/ciharden-scan.sh \
    "r081-floor-$p" ref "$FS" "$N" "dev/rel081-log/floor-$p.jobs"
fi
rm -rf "$SCR"
echo "FLOOR DONE"
