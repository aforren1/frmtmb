#!/usr/bin/env bash
# Lane optrap: the fragility scan (dev/ciharden-scan.sh) with the fit
# logger (dev/optrap-profile.R) on.
#   bash dev/optrap-scan.sh <run> <config> <lib|base> [N]
set -u
ROOT=C:/Users/adf44/source/r/frmtmb-wt-optrap
export R_PROFILE_USER="$ROOT/dev/optrap-profile.R"
export OPTRAP_PROF_LOG="$ROOT/dev/optrap-log/prof-$1"
mkdir -p "$OPTRAP_PROF_LOG"
rm -f "$OPTRAP_PROF_LOG"/*.txt
bash "$ROOT/dev/ciharden-scan.sh" "$@"
