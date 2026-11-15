#!/usr/bin/env bash
# The floors of the 0.83.0 release, measured as dev/rel082-floor.sh did:
# an extension as built into rellib-r22 is copied into a scratch library
# and loaded over the older core of rellib-r21 (0.82.0), with its floor
# lowered so that it loads, and its own test files are run there one per
# process. A floor is needed where files fail there and pass on rellib-r22.
#   frmtmb.sample 0.30.0 (floor 0.83.0)
#   frmtmb.spatial 0.5.0 (floor 0.83.0)
#   frmtmb.latent 0.7.5 (floor 0.83.0; tests only)
# First the seams: the core functions the new floors rest on, at ec8178d3
# (0.82.0) and in this tree.
#
#   bash dev/rel083-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/2}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/2}"
export GIT_OPTIONAL_LOCKS=0
. C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-env.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel083
SCR=C:/Users/adf44/AppData/Local/Temp/2/claude/ee0bfefdf562474c/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/floor
R21=C:/Users/adf44/source/r/rellib-r21
RNEW=C:/Users/adf44/source/r/rellib-r22
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
N="${1:-1}"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "== the seams, ec8178d3 against this tree"
for fn in ce_snl latres_err ce_grids_build; do
  a=$(git grep -h "^$fn <- function" ec8178d3 -- R | grep -c . || true)
  b=$(cat R/*.R | tr -d '\r' | grep -c "^$fn <- function")
  printf "  %s(): ec8178d3 %s, this tree %s\n" "$fn" "$a" "$b"
done
printf "  frmtmb_register_autocor(margin =): ec8178d3 %s, this tree %s\n" \
  "$(git show ec8178d3:R/autocor-seam.R | tr -d '\r' | grep -c 'margin')" \
  "$(tr -d '\r' < R/autocor-seam.R | grep -c 'margin')"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(c(m\$Depends, m\$Imports), function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
jobs_of() {  # <pkg> -> dev/rel083-log/floor-<pkg>.jobs
  : > "dev/rel083-log/floor-$1.jobs"
  for f in extensions/$1/tests/testthat/test-*.R; do
    echo "$1 $f" >> "dev/rel083-log/floor-$1.jobs"
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
      "r083-floor-$tag-$p" ref "$d" "$N" "dev/rel083-log/floor-$p.jobs"
  fi
}
arm frmtmb.sample "$R21" r21 yes
arm frmtmb.spatial "$R21" r21 yes
arm frmtmb.latent "$R21" r21 yes
rm -rf "$SCR"
echo "FLOOR DONE"
