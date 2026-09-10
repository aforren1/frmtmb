# Handing a round to a new session

Written 2026-09-10, at the round 3 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: this machine has lost its R library four times in ten days,
and the fourth loss settled the cause.

## Where the tree stands

Main carries the round 3 release. Four lanes merged: `wt-nss`,
`wt-rlddm`, `wt-gddm` and `wt-latent`. The release commit and the docs
rebuild are separate commits, in that order. **The user pushes; no
session pushes for them.**

Versions: frmtmb **0.55.2, unchanged**, frmtmb.eam 0.8.0, frmtmb.ode
0.4.0, frmtmb.learn 0.4.0, frmtmb.latent 0.3.0, frmtmb.sample 0.4.1,
frmtmb.spline 0.5.0, frmtmb.coupling 0.3.0. Every extension floors on
frmtmb 0.55.1. frmtmb.learn now floors on frmtmb.eam 0.8.0, raised
this round because `rlddm()` takes its non-decision-time bound through
a seam that frmtmb.eam exports for the first time at 0.8.0.

**Core is not a lane right now.** The user is doing brms compatibility
in core separately, and a round that opens a core lane collides with
that work. Ask before opening one.

### Verified at the release commit

The harness is `dev/release/`, in the repo. Every tier ran against one
private library holding all eight packages.

- Suite, one R process per test file: 220 files, 13124 assertions, no
  failures and no errors. One file moved against 0.55.2 and the reason
  is known: `test-tmbstan-build-guard.R` goes 24 to 23, because the
  StanHeaders pin removes a conditional `succeed()`.
- Gated tiers, NO SKIPS: 23 of 23 files, 2475 assertions. Every row is
  identical to 0.55.2: BCM 363, brms likelihood 404, methods 950,
  agreement 187, priors 103, port 18, reinforcement-learning identity
  102, learn Stan identity 55, sample loo 79, sampling-ported 212.
- Scale tier: 7 of 7 files, 46 assertions, no skips. The rows this
  round moved: `learn-rlddm` at convergence code 0, maximum gradient
  0.00312, a positive definite Hessian, no bad standard errors, and
  100 of 100 learners below their own fastest response, where 0.55.2
  had code 1 at 6.36e+09 with four NaN standard errors; `eam` at
  -7003.01 with 30 of 30 below their own floor and 27.3 ms to spare.
  The `ode` row carries convergence code 1 at a maximum gradient of
  0.00224 with a positive definite Hessian, which it also did at
  0.55.2. That is a carried condition, not a regression, and it is the
  one row in the tier a future session should not read as clean.
- `R CMD check --as-cran`: 8 of 8. Five carry one NOTE, this machine's
  environmental V8 math-rendering note, and three are OK. This run did
  NOT pass `--no-manual`, so the note is exercised rather than skipped
  and the result is comparable to the lanes' own checks.
- Docs: 8 of 8 sites.

**Read the comment at the top of each `dev/release/` script before you
change it.** Each records a specific way a release was measured wrong
here, and the flags they avoid are avoided on purpose.

## What is next, in order

`dev/extension-gaps-plan.md` is the plan, and it is edited in place
rather than appended to.

**Phase 1 is complete.** It came in at about 13 days against an
estimate of 10, and the whole overrun is three items that were NOT in
the plan when the phase started: 1.0d, 1.0e, and the export seam that
1.0b turned out to need. All three were silent wrong answers, and all
three were found by measuring something next to them. The estimate was
not wrong about the work it knew about.

Phase 2 has four rows left, about 4.5 days, and they run as four
parallel lanes over four packages:

- **2.1**, eam hierarchy at 30 x 400 over 60 replicates. 1.5 days.
- **2.2**, learn `(1 | p | id)` at 100 x 200 over 60 replicates. 1.5
  days. Read its three warnings BEFORE you cost the lane: 60
  replicates of rlddm at 12.3 minutes each is 12 hours, and each of
  the three would otherwise spend that time on a wrong assertion.
- **2.5**, spline frailty against rstpm2. 1 day.
- **2.6**, coupling recovery at the realistic design. 0.5 days.

Two warnings repeat across 2.1 and 2.2 and are worth carrying into
either lane. Do NOT assert `sd(ndt)`: it is worse than switching the
component off, 0.0334 against 0.0380 with a truth of 0.0393, so assert
the per-subject error and the log-likelihood. And do not assert a
diagonal truth on a grouped arm: the design draws independent
deviations, but the block the model fits is not the block the design
drew, so the correlation in the fitted parameterization is about -0.78
before any fit and is computable per replicate with no fit at all.

Item 2.4 is the reason to distrust a plan row that predicts nothing
will be found. It said "none expected", and the measurement found a
starting rule that reached a local optimum 243 to 284 log-likelihood
units below poLCA on 8 of 200 replicates, seven of which passed the
convergence pair this plan's own scale tier asserts.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes. A reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a shared reference build
of the base commit. Punch rounds go back to the worker, and the SAME
reviewer re-checks, because it holds the harness it already built.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Two punch rounds is the cap, with a third only for a blocker. Build
ONE reference library for the whole round rather than one per
reviewer. `dev/organizer-rules.md` has the rest, including the five
cost rules.

Consolidation: commit each lane on its branch, merge, set versions,
roxygenise, THEN install all eight into one private library (that
order, because installing first once shipped a truncated NAMESPACE),
run the suite and the gated tiers and the scale tier and the checks
and the docs, commit the release and the docs separately, remove the
worktrees, prune merged branches, regenerate `dev/suite-baseline.tsv`.

## What the user has settled, so a session does not reopen it

- **The user pushes.** No session pushes for them, and no git
  operation at all without their say so.
- **StanHeaders is pinned OUTSIDE the user library**, in
  `C:/Users/adf44/source/r/pinlib` at 2.32.10, with the tarball beside
  it. The user library keeps 2.39.1 and is not to be downgraded.
  `dev/lane-rules.md` has the `.libPaths()` order and the two checks
  worth running before believing a Stan-backed tier. A later worker
  may prefer renv or the Posit package manager to install the earlier
  StanHeaders; the pin is what exists today.
- **The R user library stays under `%LOCALAPPDATA%`.** The user knows
  it has been destroyed four times and chose to keep it there for now.
  The canary is in place, and `dev/machine-library.md` has the restore
  recipe and the diagnosis.
- **The tmbstan defect gets no disclosure.** Nothing external depends
  on this package. The guard's own error message already tells a user
  to distrust draws from an affected installation, and that stays.
- **Backward compatibility breaks freely**, since nothing external
  depends on this yet, but the user wants to hear about each break,
  and each needs a NEWS bullet saying plainly what stops working. This
  round broke one thing: a constant `ndt` outside the link's range is
  now refused at parse rather than reaching the optimizer as a NaN.

## Open decisions that belong to the user, not to a session

- Any change that alters what a shipped parameter MEANS.
- Where a disclosure goes, on the rare occasion one is wanted. The
  0.6.0 `ndt` disclosure went in frmtmb.eam's NEWS only, at the user's
  direction, and is discharged as of 0.7.0.
- Whether a round opens a lane in core, while brms compatibility is
  the user's own work.

## What round 3 is evidence for

Four lanes and four reviewers found six defects in shipped code. The
consolidating session, running the release harness, produced three
more broken measurements of its own, which makes eight across the two
consolidations. None of the eight failed loudly. Every one reported a
plausible number.

Round 3 added these three: a runner that called `test_file()` from the
global environment, where a symbol brought in by `importFrom` is
invisible; a driver whose scripts lived in a session temp directory
that was cleaned mid-release, after which `powershell -File <missing>`
exited 0 and the job reported success with no log; and
`$ErrorActionPreference = "Stop"` with `2>&1`, which turned an
ordinary rstan stderr line into a terminating error and killed a gated
tier thirteen files in.

The generalization from round 2 still holds and got stronger: **a
release harness is a guard, and every guard built in the last three
rounds failed open on its first try.** Apply the lane evidence
standard to the harness and not only to lane code: read the file count
and the skip count before the failure count, and construct the case
where the thing you are measuring is ABSENT.

Round 3 adds a second one, about review. Measurement moved the
headline in BOTH directions, on all four lanes:

- The rlddm reviewer proved the lane's headline defect did not exist.
  `bf(ndt = 0.2)` fitted 0.2 exactly, and the lane's script had
  reconstructed an eta that nothing uses. The lane retracted.
- The latent reviewer found `hmm_starts()` false-alarming on 6 of 6
  unimodal chains, then measured its own recommended remedy and proved
  it does not work, 87.3 percent against 87.3 percent.
- The gddm lane built the predicate the organizer pushed for, measured
  it, and DECLINED it, on the ground that Rule 3 does not license
  trading a reachable refusal for an unreachable acceptance.
- The nss lane rejected its reviewer's `"auto"` proposal by
  measurement, and separately smoothed a discontinuity that the
  organizer ranked above the reviewer's own ordering.

A review that only confirms is not doing the job. Neither is one that
only rejects. Both directions were load-bearing this round.

## Worktrees

`wt-nss`, `wt-rlddm`, `wt-gddm` and `wt-latent` are merged and
removed, and their evidence scripts, logs and tables are committed on
main under `dev/`. Create fresh worktrees off the current main.
