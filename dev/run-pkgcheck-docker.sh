#!/usr/bin/env bash
# Run rOpenSci's pkgcheck container against this source tree, the same way
# .github/workflows/pkgcheck.yaml runs it on GitHub. See dev/pkgcheck-docker.md.
#
# Usage:
#   dev/run-pkgcheck-docker.sh                  # full pkgcheck, gates off (mirrors CI)
#   dev/run-pkgcheck-docker.sh --gates on       # full pkgcheck, chain-agreement gates on
#   dev/run-pkgcheck-docker.sh --tests 'v07|v15|sample-direct'   # tests only, no pkgcheck
#   dev/run-pkgcheck-docker.sh --tests 'v07' --gates on --repeat 2
#
# Run this from anywhere. It stages a clone of the committed HEAD, so
# uncommitted work is NOT checked: commit first, then run the gate.

set -euo pipefail

# MSYS rewrites container-side paths such as /github/workspace into
# C:/Program Files/Git/github/workspace, so docker needs the rewrite off. Only
# docker: git is a Windows binary too, and it needs the rewrite ON to read a
# /c/... argument, so the switch stays on the single call it belongs to.
dk() { MSYS_NO_PATHCONV=1 docker "$@"; }

IMAGE="ghcr.io/ropensci-review-tools/pkgcheck-action:latest"
DEPS_IMAGE="frmtmb-pkgcheck:deps"
CACHE_VOL="frmtmb-pkgcheck-cache"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d-%H%M%S)"

GATES="off"
TESTS=""
REPEAT=1
REUSE_DEPS=0

while [ $# -gt 0 ]; do
  case "$1" in
    --gates)  GATES="$2"; shift 2 ;;
    --tests)  TESTS="$2"; shift 2 ;;
    --repeat) REPEAT="$2"; shift 2 ;;
    --reuse-deps) REUSE_DEPS=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [ "$GATES" = "off" ]; then
  GATE_ENV="false"
else
  GATE_ENV="true"
fi

MODE="pkgcheck"
[ -n "$TESTS" ] && MODE="tests"
LOG="$REPO_ROOT/dev/pkgcheck-docker-$MODE-gates-$GATES-$STAMP.log"

# Stage a clone. The container writes .pkgcheck/, summary.md and full.md into
# its workspace and edits .Rbuildignore, so it must not get the real tree.
STAGE_ROOT="${TMPDIR:-/tmp}/dock-pkgcheck-$STAMP"
STAGE="$STAGE_ROOT/frmtmb"
mkdir -p "$STAGE_ROOT"
git clone --no-hardlinks -q "$REPO_ROOT" "$STAGE"
STAGE_WIN="$(cd "$STAGE" && pwd -W 2>/dev/null || echo "$STAGE")"

cleanup() { rm -rf "$STAGE_ROOT" >/dev/null 2>&1 || true; }
trap cleanup EXIT

# pkgcheck asks GitHub whether the package name is free and whether CI passes.
# Without a token both answers are rate-limited guesses, so hand it one.
GH_TOKEN_VAL="$(gh auth token 2>/dev/null || echo "")"

dk volume create "$CACHE_VOL" >/dev/null


# .github/pkgcheck.Rprofile ends with Sys.setenv(FRMTMB_SAMPLER_GATES="false"),
# which silently overrides -e. The profile exists to add the kaskr repository
# for RTMBode, and a --tests run needs no repositories because the dependency
# image already holds them, so a gates-on run drops the profile instead.
PROFILE_ARGS=(-e R_PROFILE_USER=.github/pkgcheck.Rprofile)
if [ "$GATES" = "on" ]; then
  if [ "$MODE" = "tests" ]; then
    # An empty value, not an absent one: docker commit bakes every -e of the
    # committed run into the image, so the dependency image can carry the
    # profile path already. Blank it rather than leave it out.
    PROFILE_ARGS=(-e R_PROFILE_USER=)
  else
    echo "--gates on with a full pkgcheck run needs a profile without the" >&2
    echo "Sys.setenv() line: the workflow profile forces the gates off." >&2
    exit 4
  fi
fi

COMMON_ARGS=(
  --rm
  -v "$STAGE_WIN:/github/workspace"
  -v "$CACHE_VOL:/root/.cache"
  -w /github/workspace
  "${PROFILE_ARGS[@]}"
  -e "FRMTMB_SAMPLER_GATES=$GATE_ENV"
  -e NOT_CRAN=true
  -e CI=true
  -e ROPENSCI=true
  # octolog turns a failed check into an annotation and carries on only when it
  # believes it is inside Actions. Without this the entrypoint aborts on the
  # first failed check and never prints the Check Results block.
  -e GITHUB_ACTIONS=true
  -e "GITHUB_TOKEN=$GH_TOKEN_VAL"
)

{
  echo "### frmtmb pkgcheck container run"
  echo "### mode=$MODE gates=$GATES (FRMTMB_SAMPLER_GATES=$GATE_ENV)"
  echo "### HEAD=$(git -C "$STAGE" rev-parse --short HEAD)"
  echo "### image=$IMAGE"
  echo "### started $(date -Iseconds)"
  echo
} | tee "$LOG"

if [ "$MODE" = "pkgcheck" ]; then
  RUN_IMAGE="$IMAGE"
  [ "$REUSE_DEPS" = "1" ] && RUN_IMAGE="$DEPS_IMAGE"
  # No --rm here: the finished container holds the installed dependency tree,
  # which --tests runs reuse instead of paying pak's install cost again.
  CNAME="frmtmb-pkgcheck-$STAMP"
  set +e
  dk run "${COMMON_ARGS[@]:1}" --name "$CNAME" "$RUN_IMAGE" 2>&1 | tee -a "$LOG"
  RC=${PIPESTATUS[0]}
  set -e
  # summary.md and full.md are written before the entrypoint reports failures,
  # so keep them even when the run ends badly.
  for f in summary.md full.md; do
    [ -f "$STAGE/.pkgcheck/$f" ] &&
      cp "$STAGE/.pkgcheck/$f" "$REPO_ROOT/dev/pkgcheck-docker-$STAMP-$f"
  done
  # docker commit copies the run's -e values into the image config, which would
  # store the GitHub token in a local image and pin the gate switch. Clear all
  # three on the way in.
  dk commit \
    --change 'ENV GITHUB_TOKEN=' \
    --change 'ENV R_PROFILE_USER=' \
    --change 'ENV FRMTMB_SAMPLER_GATES=' \
    "$CNAME" "$DEPS_IMAGE" >/dev/null &&
    echo "### committed deps image $DEPS_IMAGE" | tee -a "$LOG"
  dk rm -f "$CNAME" >/dev/null
else
  # Tests only. Needs the dependency image built by a prior pkgcheck run.
  if ! dk image inspect "$DEPS_IMAGE" >/dev/null 2>&1; then
    echo "no $DEPS_IMAGE yet: run without --tests first" >&2; exit 3
  fi
  RSCRIPT="
    library(testthat)
    for (i in seq_len($REPEAT)) {
      cat('##### repeat', i, 'of', $REPEAT, '\n')
      print(Sys.getenv('FRMTMB_SAMPLER_GATES'))
      res <- test_dir('tests/testthat', filter = '$TESTS',
                      package = 'frmtmb', load_package = 'installed',
                      stop_on_failure = FALSE, reporter = 'summary')
    }
  "
  set +e
  dk run "${COMMON_ARGS[@]}" --entrypoint Rscript "$DEPS_IMAGE" \
    -e "$RSCRIPT" 2>&1 | tee -a "$LOG"
  RC=${PIPESTATUS[0]}
  set -e
fi

{
  echo
  echo "### finished $(date -Iseconds) exit=$RC"
} | tee -a "$LOG"

echo
echo "log: $LOG"
exit "$RC"
