#!/usr/bin/env bash
# dev/ciharden-scan.sh with lane optstart's fit logger
# (dev/optstart-profile.R) as R_PROFILE_USER, so that the release
# build's run of the suites can be compared fit by fit with the lane's
# (dev/optstart-fitdiff.R, which reads dev/optstart-log/prof-<run>/).
#   bash dev/rel071-scan-prof.sh <run> <config> <lib|base> [N] [list]
set -u
ROOT=C:/Users/adf44/source/r/frmtmb-wt-release
export R_PROFILE_USER="$ROOT/dev/optstart-profile.R"
export OPTSTART_PROF_LOG="$ROOT/dev/optstart-log/prof-$1"
mkdir -p "$OPTSTART_PROF_LOG"
rm -f "$OPTSTART_PROF_LOG"/*.txt
bash "$ROOT/dev/ciharden-scan.sh" "$@"
