#!/usr/bin/env bash
# The floors of the 0.77.0 release, measured. Adapted from
# dev/rel076-floor.sh, with its method: nothing is installed. Each arm
# loads a COPY of an installed build (none of these packages has
# compiled code) in front of the core it is measured on.
#
# frmtmb.sample (floor raised to frmtmb (>= 0.77.0)):
#   old: rellib-r15's build, its floor lowered to 0.61.0 in DESCRIPTION
#        and Meta/package.rds, in front of rellib-r14 (frmtmb 0.76.0);
#        if it loads, every test file of its suite runs there (the
#        release tree's files, test-prior-draws.R among them), one per
#        process, every gate on. On 0.76.0 core, sample_prior = "only"
#        samples the posterior (priordraws' review, finding 3), which
#        test-prior-draws.R must catch.
#   new: rellib-r15 as built.
#   converse: rellib-r14's frmtmb.sample 0.23.0 in front of rellib-r15.
# frmtmb.spatial (no R file changed; one test helper did): its 8 test
# files from the release tree on rellib-r14 alone (frmtmb 0.76.0 with
# frmtmb.spatial 0.3.0). Equal counts to rellib-r15's tier mean its
# tests do not need the new core, and the version stays.
#
#   bash dev/rel077-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/12}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/12}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel077
SCR=C:/Users/adf44/source/r/rel077-floor-lib
R14=C:/Users/adf44/source/r/rellib-r14
R15=C:/Users/adf44/source/r/rellib-r15
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "== 435e7ce6 (0.76.0): the seam the new floor rests on"
printf '  prior_only in R/objective.R: 435e7ce6 %s, this tree %s\n' \
  "$(git show 435e7ce6:R/objective.R | grep -c 'prior_only')" \
  "$(grep -c 'prior_only' R/objective.R)"
printf '  "flat" density kind in R/priors.R: 435e7ce6 %s, this tree %s\n' \
  "$(git show 435e7ce6:R/priors.R | grep -c 'identical(dist$kind, "flat")')" \
  "$(grep -c 'identical(dist$kind, "flat")' R/priors.R)"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(c(m\$Depends, m\$Imports), function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
p=frmtmb.sample
v=$(grep '^Version:' "extensions/$p/DESCRIPTION" | cut -d' ' -f2)
fl=$(grep -o 'frmtmb (>= [0-9.]*)' "extensions/$p/DESCRIPTION" | head -1)
echo "== $p $v, floor $fl"
FL=$SCR/old-$p; rm -rf "${SCR:?}/old-$p"; mkdir -p "$FL"
cp -r "$R15/$p" "$FL/"
echo "-- old: rellib-r15's build over rellib-r14 (floor lowered)"
lower_floor "$FL/$p"
out=$(load_line "$FL" "$R14" "$p"); echo "$out"
if echo "$out" | grep -q "LOADED"; then
  : > "dev/rel077-log/floor-$p.jobs"
  for f in extensions/$p/tests/testthat/test-*.R; do
    echo "$p $f" >> "dev/rel077-log/floor-$p.jobs"
  done
  ROUND_BASE=$R14 /usr/bin/bash dev/ciharden-scan.sh \
    "r077-floor-$p" ref "$FL" "${1:-8}" "dev/rel077-log/floor-$p.jobs"
fi
echo "-- new: rellib-r15 as built"
load_line "$R15" "$R15" "$p"
FC=$SCR/conv-$p; rm -rf "${SCR:?}/conv-$p"; mkdir -p "$FC"
cp -r "$R14/$p" "$FC/"
echo "-- converse: $p $(grep '^Version:' "$FC/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r14 over rellib-r15"
load_line "$FC" "$R15" "$p"
echo "== frmtmb.spatial: the release tree's tests on rellib-r14"
: > dev/rel077-log/floor-spatial.jobs
for f in extensions/frmtmb.spatial/tests/testthat/test-*.R; do
  echo "frmtmb.spatial $f" >> dev/rel077-log/floor-spatial.jobs
done
ROUND_BASE=$R14 /usr/bin/bash dev/ciharden-scan.sh r077-floor-spatial ref \
  base "${1:-8}" dev/rel077-log/floor-spatial.jobs
rm -rf C:/Users/adf44/source/r/rel077-floor-lib
echo "FLOOR DONE"
