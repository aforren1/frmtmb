#!/usr/bin/env bash
# Lane optstart: the fragility scan (dev/ciharden-scan.sh), one test file
# per R process, with the fit logger (dev/optstart-profile.R) on.
#   bash dev/optstart-scan.sh <run> <config> <lib|base> [N] [list]
# Fit logs: dev/optstart-log/prof-<run>/<pid>.txt; per-file results:
# dev/ciharden-log/<run>-suite/.
set -u
ROOT=C:/Users/adf44/source/r/frmtmb-wt-optstart
export R_PROFILE_USER="$ROOT/dev/optstart-profile.R"
export OPTSTART_PROF_LOG="$ROOT/dev/optstart-log/prof-$1"
mkdir -p "$OPTSTART_PROF_LOG"
rm -f "$OPTSTART_PROF_LOG"/*.txt
bash "$ROOT/dev/ciharden-scan.sh" "$@"
