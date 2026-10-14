#!/usr/bin/env bash
# The whole brms vignette re-measure against one installed build, as one
# command: the mechanical port (raw, spell, need and keep-prior passes),
# the hand translation, the family coverage, the estimate comparison
# with brms, and the reducers that fill the record.
#
#   bash dev/remeasure-run.sh <libs> <tag> [jobs]
#
# <libs>  ";"-separated libraries put before the user library (PORT_LIB
#         and BV_LIB). For a release build that is the release library
#         alone, for example C:/Users/adf44/source/r/rellib-r8.
# <tag>   names the output directory dev/remeasure-out/<tag>, and the
#         column the build gets in dev/brms-port/outcomes-by-build.tsv.
# [jobs]  R processes at once (default 8: the machine is shared).
#
# The brms reference fits do not depend on the frmtmb build, so they are
# fitted once into dev/remeasure-out/brms-ref and reused by every tag.
# Set REMEASURE_BRMS=1 to (re)fit them; the first run fits them when the
# directory is missing. They take most of an hour, on 4 processes.
#
# Wall time per stage and in total goes to dev/remeasure-out/<tag>/wall.txt.
set -u
wt="$(cd "$(dirname "$0")/.." && pwd)"
lib="$1"
tag="$2"
jobs="${3:-8}"
# [stages] "all" (default), "ports" (ports only), "brms" (reference fits
# only) or "reduce" (reducers only, on the records already there).
stages="${4:-all}"
rs="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out="$wt/dev/remeasure-out/$tag"
ref="$wt/dev/remeasure-out/brms-ref"
mkdir -p "$out/port" "$out/bv"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$wt/dev/stan-cache}"
export NOT_CRAN=true
# Rtools' bash drops TMP, and R then cannot make its temporary directory.
export TMP="${TMP:-${TEMP:-/tmp}}" TEMP="${TEMP:-${TMP}}"
export PORT_LIB="$lib" PORT_OUT="$out/port"
export BV_LIB="$lib" BV_OUT="$out/bv" BV_DIR="$wt/dev/brms-vignettes"
export PORT_PATCHES="brms_nonlinear.9.1,brms_multilevel.15.3"
vigs="brms_overview brms_multilevel brms_distreg brms_nonlinear
brms_phylogenetics brms_monotonic brms_multivariate brms_missings
brms_customfamilies"
wall="$out/wall.txt"
t_all=$(date +%s)
echo "## $(date '+%Y-%m-%d %H:%M:%S') stages=$stages jobs=$jobs lib=$lib" \
  >> "$wall"
stamp() { echo "$1 $(( $(date +%s) - $2 )) s" >> "$wall"; }
want() { [ "$stages" = all ] || [ "$stages" = "$1" ]; }

# One job per line: "port <vignette> <mode>" or "bv <script path>".
run_job() {
  case "$1" in
    port) "$rs" "$wt/dev/brms-port/run-vignette.R" "$2" 120 "$3" \
            > "$PORT_OUT/run-$2-$3.out" 2>&1 ;;
    bv)   "$rs" "$2" > "$BV_OUT/bv-$(basename "$2" .R).log" 2>&1 ;;
    brms) PORT_OUT="$ref" "$rs" "$wt/dev/brms-port/brms-fit.R" "$2" "$3" 2 \
            > "$ref/brms-$2-$3.out" 2>&1 ;;
  esac
}
export -f run_job
export rs wt ref

t0=$(date +%s)
"$rs" "$wt/dev/remeasure-provenance.R" > "$out/provenance.txt" 2>&1
stamp provenance "$t0"

if want ports; then
  t0=$(date +%s)
  {
    for m in raw spell need keepprior; do
      for v in $vigs; do echo "port $v $m"; done
    done
    for f in "$wt"/dev/brms-vignettes/brms_*.R \
             "$wt"/dev/brms-vignettes/brms-*.R; do
      echo "bv $f"
    done
  } | xargs -P "$jobs" -L 1 bash -c 'run_job "$@"' _
  for f in "$BV_OUT"/bv-*.log; do
    printf '%-36s %4s ok %4s err  %s\n' "$(basename "$f")" \
      "$(grep -c '^\[ok' "$f")" "$(grep -c '^\[err' "$f")" \
      "$(grep -m1 '^### done' "$f" || echo 'NO DONE LINE')"
  done > "$out/bv/done.txt"
  stamp "ports (4 mechanical passes and the hand translation)" "$t0"
fi

if [ "$stages" = brms ] || { [ "$stages" = all ] &&
     { [ "${REMEASURE_BRMS:-0}" = 1 ] || [ ! -d "$ref/brms" ]; }; }; then
  t0=$(date +%s)
  mkdir -p "$ref"
  {
    for v in $vigs; do echo "brms $v 1000"; done
    # brms_overview's fits do not mix at the cheap budget
    echo "brms brms_overview 0"
  } | xargs -P 3 -L 1 bash -c 'run_job "$@"' _
  stamp "brms reference fits" "$t0"
fi

if want reduce; then
  t0=$(date +%s)
  "$rs" "$wt/dev/brms-port/summarize.R" > "$out/port/summary.txt" 2>&1
  "$rs" "$wt/dev/brms-vignettes/_drift.R" > "$out/bv/drift.txt" 2>&1
  "$rs" "$wt/dev/brms-port/families-coverage.R" > "$out/families.txt" 2>&1
  PORT_BRMS="$ref" "$rs" "$wt/dev/brms-port/plausibility.R" \
    > "$out/plausibility.txt" 2>&1
  "$rs" "$wt/dev/remeasure-history.R" "$tag" > "$out/history.txt" 2>&1
  "$rs" "$wt/dev/remeasure-gaps.R" "$tag" > "$out/gaps.txt" 2>&1
  # the estimate comparison against every earlier build still on disk
  cmp=""
  for t in v067 r6 r7 r8; do
    if [ "$t" != "$tag" ] &&
       [ -f "$wt/dev/remeasure-out/$t/port/plausibility.rds" ]; then
      cmp="$cmp $t"
    fi
  done
  "$rs" "$wt/dev/remeasure-plauscmp.R" $cmp "$tag" > "$out/plauscmp.txt" 2>&1
  "$rs" "$wt/dev/remeasure-secmp.R" $cmp "$tag" > "$out/secmp.txt" 2>&1
  stamp reducers "$t0"
fi
stamp total "$t_all"
cat "$wall"
echo "remeasure $tag done"
