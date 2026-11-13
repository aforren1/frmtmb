#!/usr/bin/env bash
# The floors of the 0.82.0 release, measured as dev/rel081-floor.sh did:
# an extension as built into rellib-r21 is copied into a scratch library
# and loaded over an older core, with its floor lowered where it would
# refuse to load, and its own test files are run there one per process.
#   frmtmb.sample 0.29.0 (floor 0.82.0) over rellib-r20 (0.81.0)
#   frmtmb.spatial 0.4.0 (floor 0.82.0) over rellib-r20
#   frmtmb.survival 0.3.0 (floor 0.75.0, as built) over rellib-r20 and
#     rellib-r19 (0.80.0); no older core library is left to read
#   frmtmb.latent 0.7.4 (floor 0.80.0, as built) over rellib-r19
# First the seams: the core functions the new floors rest on, at
# e8abed0a and in this tree.
#
#   bash dev/rel082-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/2}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/2}"
export GIT_OPTIONAL_LOCKS=0
. C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/rel082-env.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel082
SCR=C:/Users/adf44/AppData/Local/Temp/2/claude/ee0bfefdf562474c/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/floor
R19=C:/Users/adf44/source/r/rellib-r19
R20=C:/Users/adf44/source/r/rellib-r20
RNEW=C:/Users/adf44/source/r/rellib-r21
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
N="${1:-1}"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "== the seams, e8abed0a against this tree"
for fn in integrated_design_zero brms_variable_sort autocor_loo_dense_max; do
  a=$(git grep -h "^$fn <- function" e8abed0a -- R | grep -c . || true)
  b=$(cat R/*.R | tr -d '\r' | grep -c "^$fn <- function")
  printf "  %s(): e8abed0a %s, this tree %s\n" "$fn" "$a" "$b"
done
for fn in integrated_design_zero brms_variable_sort; do
  printf "  frmtmb.sample calls %s(): %s lines\n" "$fn" \
    "$(cat extensions/frmtmb.sample/R/*.R | tr -d '\r' | grep -c "$fn(")"
done
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(c(m\$Depends, m\$Imports), function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
jobs_of() {  # <pkg> -> dev/rel082-log/floor-<pkg>.jobs
  : > "dev/rel082-log/floor-$1.jobs"
  for f in extensions/$1/tests/testthat/test-*.R; do
    echo "$1 $f" >> "dev/rel082-log/floor-$1.jobs"
  done
}
# <pkg> <old core lib> <tag> <lower the floor: yes|no>
arm() {
  local p=$1 base=$2 tag=$3 low=$4 d=$SCR/$3-$1
  echo "== $p $(grep '^Version:' extensions/$p/DESCRIPTION | cut -d' ' -f2 | tr -d '\r'), floor $(grep -o 'frmtmb (>= [0-9.]*)' extensions/$p/DESCRIPTION | head -1); over $base"
  rm -rf "$d"; mkdir -p "$d"
  cp -r "$RNEW/$p" "$d/"
  [ "$low" = yes ] && lower_floor "$d/$p"
  out=$(load_line "$d" "$base" "$p"); echo "$out"
  if echo "$out" | grep -q "LOADED"; then
    jobs_of "$p"
    ROUND_BASE=$base /usr/bin/bash dev/ciharden-scan.sh \
      "r082-floor-$tag-$p" ref "$d" "$N" "dev/rel082-log/floor-$p.jobs"
  fi
}
arm frmtmb.sample "$R20" r20 yes
arm frmtmb.spatial "$R20" r20 yes
arm frmtmb.survival "$R20" r20 no
arm frmtmb.survival "$R19" r19 no
arm frmtmb.latent "$R19" r19 no
rm -rf "$SCR"
echo "FLOOR DONE"
