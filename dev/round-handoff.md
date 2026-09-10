# Handing a round to a new session

Written 2026-09-09, at frmtmb 0.55.2. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. `dev/machine-library.md` is new and is worth
reading before you run anything.

## Where the tree stands

Main is at the 0.55.2 release, `3ff5d8e`, plus its docs rebuild
`883dfc4`. Eleven commits sit ahead of `origin/main`. **The user
pushes; no session pushes for them.**

Versions: frmtmb 0.55.2, frmtmb.eam 0.7.0, frmtmb.ode 0.3.0,
frmtmb.sample 0.4.1, frmtmb.learn 0.3.0, frmtmb.spline 0.5.0,
frmtmb.coupling 0.3.0, frmtmb.latent 0.2.2. Every extension floors on
frmtmb 0.55.1; frmtmb.learn also floors on frmtmb.eam 0.6.0.

Verified at that commit, one test file per R process, all eight
packages installed into one library: 216 files, 12794 assertions, no
failures and no errors. Gated tiers with no skips, 23 files and 2475
assertions: BCM 363 over thirteen chapters, brms likelihood 404,
methods 950, agreement 187, priors 103, port 18,
reinforcement-learning identity 102, learn Stan identity 55, sample
loo 79 and sampling-ported 212, fuzz green. Scale tier, 7 files and 14
rows, no skips. `R CMD check --as-cran` OK on all eight with no NOTE
and no WARNING, though that run passed `--no-manual`, so the
environmental V8 note earlier releases carried was not exercised.

## What is next, in order

`dev/extension-gaps-plan.md` is the plan and it is edited in place
rather than appended to. Phase 1 is down to two items:

- **1.0b**, the per-group non-decision-time bound in `rlddm()`. It is
  unblocked: `ndt_group()` works as it stands. What it needs first is
  an export seam from frmtmb.eam for four functions that are `@noRd`
  today, whose shape depends on `rlddm()`, so 1.0a did not guess it.
  Half a day. The scale tier still records the failure it fixes:
  `learn-rlddm` at a maximum gradient of 6.36e+09 with a non positive
  definite Hessian and four bad standard errors.
- **1.0d**, the `frm_ode(n_ss =)` default, found by 1.0c. One day. It
  is a silent wrong answer during a fit and by Rule 3 it outranks the
  rest of the backlog. `?frm_ode` is already corrected; what is left
  is the default and the diagnostic. The review's unbuilt `"auto"`
  proposal is in the row and is FASTER than today's default on the
  common case.

Then Phase 2 opens. Item 2.1 is unblocked and its row carries three
warnings, of which the important one is: do not assert `sd(ndt)`. The
criterion is satisfiable without the model doing any work, because the
observed floors carry the spread. Assert the per-subject error and the
log-likelihood.

One unnumbered item is filed in the Phase 1 prose: `gddm()` reads every
dpar at the first row of its condition, so a dpar that varies within a
condition is silently ignored. That is wider than the `ndt` case that
exposed it.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes; a reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a shared reference build
of the base commit; punch rounds go back to the worker; the same
reviewer re-checks, because it holds the harness it already built.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Two rounds is the cap, with a third only for a blocker. Build ONE
reference library for the whole round rather than one per reviewer.
`dev/organizer-rules.md` has the rest, including the five cost rules.

Consolidation: commit each lane on its branch, merge, set versions,
roxygenise, install all eight into one private library, run the suites
and the gated tiers and the scale tier and the checks and the docs,
commit the release and the docs separately, remove the worktrees,
prune merged branches, regenerate `dev/suite-baseline.tsv`.

## The thing this round is actually evidence for

Three lanes and three reviewers found six defects in shipped code. The
consolidating session, running the release harness, produced FIVE
broken measurements of its own:

- a runner whose file paths were empty, so 216 files load-errored and
  the summary read `pass=0 fail=0 err=0` for all eight packages;
- a runner that never attached the packages, so every test reported
  `could not find function` and the damage looked like a regression in
  the one package the round had rewritten;
- an `R CMD build --no-build-vignettes` flag that manufactured two
  WARNINGs and a NOTE on all eight packages, including four the round
  never touched;
- a library restore that installed StanHeaders 2.39.1 against rstan
  2.32.7, which a populated Stan cache hid everywhere except the two
  blocks that compile something new;
- a PowerShell driver whose `$env:PATH` used forward slashes so `cmd`
  could not be resolved, which ran nothing and still wrote its
  completion marker.

Four of the five reported plausibly. None failed loudly. The evidence
standard in `dev/lane-rules.md` is written for lane code and every one
of these would have been caught by applying it to the harness instead:
read the file count and the skip count before the failure count, and
construct the case where the thing you are measuring is ABSENT.

The generalization worth carrying: **a release harness is a guard, and
every guard built in the last two rounds failed open on its first
try.**

## Open decisions that belong to the user, not to a session

- Whether to push. Eleven commits are waiting.
- The R user library still holds StanHeaders 2.39.1. A future session
  hits the same wall until it is pinned to 2.32.10 there; the release
  library has the pin, the shared one does not.
- Moving the R user library off `%LOCALAPPDATA%`, which
  `dev/machine-library.md` argues for. It was destroyed three times in
  nine days and the canary `ZZZ-canary.txt` is in place to identify the
  next one.
- Whether the tmbstan defect needs a disclosure for anyone who sampled
  before the guard shipped.
- Any change that alters what a shipped parameter MEANS. The standing
  policy is to break backward compatibility freely, but the user has
  wanted to hear about each one.

## Worktrees

`wt-ndt`, `wt-lincmt` and `wt-tmbstan` are merged and removed. Create
fresh ones off the current main.
