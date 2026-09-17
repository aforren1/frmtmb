# Handing a round to a new session

Written 2026-09-17, at the 2.6c release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: this machine has lost its R library five times, and the cause
is OPEN rather than settled.

## Where the tree stands

Main carries item 2.6c, the defects brms's own test suite found. Three
lanes merged: `wt-famlink` (family objects and links), `wt-priorform`
(`bf()`, formula grammar, priors) and `wt-brmsnames` (brms's names,
`VarCorr()`, `hypothesis()`, frmtmb.sample output). A fourth worktree,
`wt-vectorize`, is NOT merged and waits on the user (see below). The
release commit and the docs rebuild are separate commits, in that order.
**The user pushes; no session pushes for them.**

Versions: frmtmb **0.59.0**, frmtmb.sample **0.7.0**, frmtmb.learn 0.4.2,
frmtmb.latent 0.3.1, frmtmb.eam 0.8.2, frmtmb.spline 0.5.2,
frmtmb.coupling 0.3.2, frmtmb.ode 0.4.1. Every extension floors on core
0.59.0, because its tests read brms's names or `varcorr_matrices()`.

### Verified at the release commit

- Suite: 231 files, 15278 assertions, no failures and no errors. Three
  counts fell, and each lane recorded the same count before merging;
  `dev/suite-baseline.md` has the details.
- The merge broke one design that each lane passed alone:
  `(1 | h:g) + (1 | g/h)`. reformulas expands `g/h` into `h:g` where
  brms writes `g:h`, so brmsnames gave two different blocks one name.
  A block from a slash now carries the flag, and brms's names read it.
  Fixed in commit 4e179c0, before the release suite.
- Gated, no skips: 23 of 23 files, 2496 assertions, and three files
  gained assertions. Scale: 7 of 7, with every logLik identical to
  0.58.0. `R CMD check --as-cran`: 8 of 8, with no warnings or errors.
  Five packages carry the environmental "V8 unavailable" NOTE.

**Read the comment at the top of each `dev/release/` script before you
change it.** Two flags are forbidden there.

## PAUSED 2026-09-17 mid-round (read this first)

The user shut the machine down with agents running. They were stopped
cleanly; nothing is half-merged. Worktrees and their state:

- `wt-conditions` (2.6e, classed conditions): worker DONE, all suites
  green, uncommitted. Its reviewer was STOPPED before reporting; rerun
  the review from scratch (brief: message fidelity on 60+ sites, the
  runtime subclass rule, census gaps, the 29 unclassed sweep cases,
  interop code checking `simpleError`). The punch round must also carry
  the user's decisions: keep warnings/messages classed; turn the
  user-reachable `stopifnot()`/`match.arg()` refusals into named
  classed refusals (`set_prior("normal(0)")`, `student_t(3, 0)`,
  `exponential(-1)`, `gamma(1)`, `lkj()`, six match.arg sites); apply
  the one-line `getME` S4 fix at `R/interop.R:421`.
- `wt-brmsport` (2.6b, bin-1 port): reviewer said MERGEABLE. The worker
  was STOPPED partway through a small last round (R5 emmeans rows back
  to cannot transfer, target totals 192/112/34/12/144; R1 caveats and a
  pattern specificity guard; R2 stale check in `brms_port_own()`; R6
  helper-drift test into the gated tier; cheap R3 hollow forms; the
  run-gated.ps1 header NIT). Check the worktree state, finish or redo
  that round, then merge. It adds tests only.
- `wt-tmbstan121`: findings DONE (`dev/tmbstan121-findings.md`), one test
  fix, uncommitted. Its DESCRIPTION floor, CI and advice-text changes
  wait for 2.6e. The user must add `CXX17FLAGS += -std=gnu++17` to their
  Makevars and reinstall tmbstan from CRAN before the pin can drop.
- `wt-vectorize`, `wt-vecshape`: deferred by the user, committed on
  their branches, worktrees removed.

Queue after 2.6e merges: 2.6d and 2.6f together (predict and return
shapes), the port's defects (eight silent ones ranked in
`dev/brmsport-findings.md`), `hmm_starts(1)`, `brmshypothesis` class,
the tmbstan floor and CI.

## What is next, in order

**Waiting on the user:** whether `frmtmb_control(vectorize = FALSE)`
merges. `wt-vectorize` measured `TapeConfig(vectorize = "enable")` on
17 models and found it SLOWER on every random-effects model, from 1.36x
(`s(x)`) to 1,116x (`ar1()`), and faster on none. It breaks the
sparsity of the Laplace sparse Hessian tape, which grows as rows times
groups. The option is built, off by default, with tests seen failing,
in that worktree, uncommitted. Merging it needs a core bump, because
frmtmb.sample calls the new `make_adfun()` export. The record is
`dev/vectorize-findings.md` there.

- **2.6d**, `predict()` becomes brms's predictive summary (user
  decision, 2026-09-16).
- **2.6e**, classed conditions, now that the three lanes are merged.
- Filed during 2.6c and not fixed. The details are in each lane's
  findings file.
  - Two prior-scope differences from brms: a `sd` prior with a group
    and no dpar reaches `phi`; mv `sd group = g` without `resp` is
    accepted.
  - Two frmtmb.sample default priors: rescor is flat where brms uses
    `lkj(1)`; the offset intercept location.
  - `mi()` is still `b_` where brms uses `bsp_`.
  - A written `theta` formula's reference component.
- Then **2.6b**, the port of bin 1, and tmbstan 1.2.1's verification in
  a private library.

**Settled, do not reopen:**
- The non-generic name collisions with brms stay, because `::` is
  sufficient in both load orders.
- A gratia older than 0.9.0, loaded before frmtmb.sample, stops it
  loading; the user accepted that.
- `inverse.gaussian` defaults to brms's `1/mu^2`, knowing 135 of 240
  designs fit cleanly on it against 234 on `log`.
- Duplicate priors on one slot, and duplicate group-level effects
  including the animal model, are refused as brms refuses them.
- `frm_simulate(newparams =)` takes brms names only.

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
