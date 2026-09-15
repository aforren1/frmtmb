# Handing a round to a new session

Written 2026-09-14, at the Phase 2 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: this machine has lost its R library five times, and the cause
is OPEN rather than settled.

## Where the tree stands

Main carries the Phase 2 release. Four lanes merged: `wt-frailty`,
`wt-eamhier`, `wt-learnhier` and `wt-coh`. The release commit and the
docs rebuild are separate commits, in that order. **The user pushes; no
session pushes for them.**

Versions: frmtmb **0.55.2, unchanged for the second round running**,
frmtmb.eam 0.8.1, frmtmb.learn 0.4.1, frmtmb.spline 0.5.1,
frmtmb.coupling 0.3.1, frmtmb.ode 0.4.0, frmtmb.latent 0.3.0,
frmtmb.sample 0.4.1. All four bumps are PATCH because no exported
function, argument, likelihood, estimate or fitted value moved in any
lane: this round was documentation, tests and measurement.

**Core is not a lane by default.** The user does brms compatibility
there separately. Phase 2.5 is the exception and they asked for it
directly; see below.

### Verified at the release commit

- Suite: 221 files, 13159 assertions, no failures and no errors.
- Gated, NO SKIPS: 23 of 23 files, 2491 assertions. Every row identical
  to the previous release except `test-stan-identity.R` at 71 from 55,
  which is the correlated-block Stan programs item 2.2 added, exactly
  the +16 the total moved by.
- Scale tier and `R CMD check --as-cran` and pkgdown: see the release
  commit message for the numbers.
- One baseline count fell, `frmtmb.spline/test-surface.R` 48 to 47, and
  it was reconciled in review before release: three assertions removed,
  two added.

**Read the comment at the top of each `dev/release/` script before you
change it.** Two flags are forbidden there for reasons this round paid
for again.

## What is next, in order

**Phase 2.5 comes before Phase 3**, and it is the only thing in the plan
that is a defect in shipped code rather than a gap. Attaching frmtmb
after brms breaks 15 of 17 generics for brms's OWN objects: frmtmb
defines rival S3 generics for 27 exported names, brms has methods on 24
of them, and frmtmb's generic wins the search path while its method
table has no `brmsfit` entry. The user reported it from their own
session and has approved 2.5a, 2.5b and 2.5d.

**Item 2.5c is a decision waiting on the user and it gates 2.5a.**
`fixef`, `ranef` and `VarCorr` have one owner, lme4; nlme does not
export them. The options are priced in the row. Do not build 2.5a
before it is answered, because the answer changes its shape.

Then Phase 3, about 10.5 days, which now includes item 3.6, the
`rp_floored()` widening filed by item 2.5.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes. A reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a shared reference build
of the base commit. Punch rounds go back to the worker, and the SAME
reviewer re-checks.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Two punch rounds is the cap, a third only for a blocker. **Build ONE
reference library for the round.** This round did not build one at all:
the previous release's library was already the base commit at the right
versions, so it served as the reference and four lanes read it
read-only. Check whether that is true again before paying for a build.

Consolidation: commit each lane on its branch, merge, set versions,
roxygenise, THEN install (that order), verify each installed NAMESPACE
against the tree, run the tiers, commit the release and the docs
separately, remove the worktrees, prune branches, regenerate
`dev/suite-baseline.tsv`.

**Run `R CMD check` on a QUIET machine.** Its examples-timing NOTE
measures load here, established with a control: a fixed arithmetic
control swings a factor of 4.7 on identical work, and the base build
itself crossed the five second threshold at 6.75 s under load.

## What the user has settled

- **The user pushes.** No session pushes, and no git operation without
  their say so.
- **Nobody uses this package yet.** Break backward compatibility
  freely: ship the refusal, bump, and say plainly in NEWS what stops
  working. No disclosure is owed beyond that, which is why item 2.2's
  correction to what `?rlddm` says needs a NEWS bullet and nothing more.
- **StanHeaders is pinned OUTSIDE `%LOCALAPPDATA%`**, in `pinlib` at
  2.32.10. It survived all five library losses. The user library keeps
  2.39.1 and is not to be downgraded.
- **The R user library stays under `%LOCALAPPDATA%`** by the user's
  decision, knowing it will recur.
- Core is the user's lane except where they ask otherwise.

## What this round is evidence for

Four lanes, four reviewers, six punch rounds. Every row's claim held,
and three of the six rows found something the row was not sent for.

**The thing worth carrying is how often an early number dissolved.**
Two signals that looked like findings did not survive their own
replicates: a coverage of 3 of 6 that became 0.85 at 20, and an `se/sd`
of 0.694 at 14 that became 0.984 at 60. The second is the instructive
one, because the lane first explained it as "a ratio of two spreads is
noisy" and the review showed the real reason: 0.694 sat at the 0.3rd
PERCENTILE of 20,000 random 14-subsets, so it was an unlucky prefix,
and the coverage count had agreed with it at 14. Both instruments were
too few. A dissolved signal is a result about the count and belongs on
the page, not in the bin.

**And a target can be the biased quantity.** Item 2.2 spent most of a
long mechanism hunt asking why an estimate sat 22 percent below a
target, refuting three candidate explanations including two of its own,
before finding that `qlogis(ndt/floor) = log(ndt) - log(m)` identically
and the target was 94 percent nuisance. Nobody had asked whether the
target was right. When an estimator looks biased against one quantity
and unbiased against eight others, ask what is different about the
quantity before asking what is wrong with the estimator.

**Guards improved.** Two failed CLOSED on their first spelling, after
three rounds in which every guard failed open. Both were written with
their inverse case at the same time as the assertion rather than after
it, and one was built deliberately to FAIL when a filed defect is
fixed, with a comment saying to flip it rather than delete it.

**The provenance failures were all one failure.** A count taken from
launches, a guard whose condition had never been observed true, and a
count taken from the plan rather than from the files. In each the thing
being measured was not the thing being reported. The third is the one
that generalizes: a lane's verifier checks files, and that was a
SENTENCE. It was closed structurally, by generating counts into the
document, and the fix found a real error on its first run.

## Worktrees

`wt-frailty`, `wt-eamhier`, `wt-learnhier` and `wt-coh` are merged and
removed, with their evidence committed on main under `dev/`. Create
fresh worktrees off the current main.
