# Handing a round to a new session

Written 2026-09-17, at the 0.60.0 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has now been lost SIX times, and the sixth has a
dated trigger.

Main carries items 2.6b (brms's own bin-1 suite ported as a gated tier),
2.6c (the defects that suite found) and 2.6e (classed conditions), plus
the tmbstan 1.2.1 verification. Four lanes merged since the last
release: `wt-brmsport`, `wt-conditions`, `wt-tmbstan121`, and before
them the three 2.6c lanes. No worktree is left. `wt-vectorize` and
`wt-vecshape` are branches, deferred by the user and never merged.

Versions: frmtmb **0.60.0**, frmtmb.sample **0.8.0**, frmtmb.eam 0.9.0,
frmtmb.learn 0.5.0, frmtmb.latent 0.4.0, frmtmb.spline 0.6.0,
frmtmb.coupling 0.4.0, frmtmb.ode 0.5.0. Every extension floors on
frmtmb 0.60.0, because each imports `frm_stop()` and the other condition
helpers.

**The StanHeaders 2.32.10 pin is GONE.** tmbstan 1.2.1 fixes the build
that sampled a standard normal, and `CXX17FLAGS += -std=gnu++17` in the
user Makevars fixes rstan's compile against StanHeaders 2.39.1. The
release scripts no longer put `pinlib` on the path; `run-tests.R`
asserts tmbstan >= 1.2.1 and that the Makevars carries the flag, and
names `R_MAKEVARS_USER` because HOME depends on the launcher.

### Verified at the release commit

- Suite: 262 files, 15567 assertions, 0 fail, 0 error, nothing below its
  baseline. 31 files are new; see `dev/suite-baseline.md`.
- Gated: 38 of 38 files, 3127 assertions, 0 skips, and 109 Stan programs
  compiled from an EMPTY cache against StanHeaders 2.39.1.
- Scale: 7 of 7, every logLik identical to 0.59.0.
- `R CMD check --as-cran`: 8 of 8, no warnings or errors, 5 packages
  carrying the environmental "V8 unavailable" NOTE and no timing NOTE.
- Two harness defects found and fixed on the way. `run-gated.ps1` had a
  loop variable `$s` that clobbered `$S`, the runner path, so the tier
  reported 0 of 38 while running nothing: PowerShell names are
  case-insensitive. And the sixth library loss, recorded in
  `dev/machine-library.md`, happened in the minute this session
  hard-killed R processes to pause a release. Do not do that.

## What is next, in order

**2.6d and 2.6f together, next.** Both reshape the post-fit methods in
the same files, so one lane.
- 2.6d: `predict()` becomes brms's predictive summary (user decision,
  2026-09-16). Audit every caller of `predict(type = )` first, including
  `conditional_effects()` and the interop packages.
- 2.6f: the return shapes of `fitted`, `residuals`, `fixef`, `ngrps`,
  `vcov` and `summary` follow brms, and frmtmb.sample offers
  `nsamples()` and `posterior_samples()` (user decision, 2026-09-17).

**Then the defects the ported suite found**, ranked with constructions
in `dev/brmsport-findings.md`. Eight are SILENT: `ar()`/`ma()` taking an
expression as the time index; `gr = g1/g2` on numeric codes grouping by
the quotient; `variables()` omitting an ordinal fit's thresholds; a
hypothesis with no relation answered as `= 0`; `fit$data`
partial-matching `fit$data2`; `point_estimate` ignored on draws;
`fitted()`/`residuals()` returning NULL on draws; ordinal residuals
answered where brms refuses. Beside them: `hmm_starts(1)`, the
`brmshypothesis` class on `hypothesis()` output (it breaks the user's
own rule), `predict(allow_new_levels = TRUE)` without the grouping
column, and `log_lik()` on a fit giving "no applicable method".

**Filed during 2.6c and not fixed**, details in each lane's findings:
two prior-scope differences from brms (a `sd` prior with a group and no
dpar reaches `phi`; mv `sd group = g` without `resp` is accepted); two
frmtmb.sample default priors (rescor flat where brms uses `lkj(1)`, and
the offset intercept location); `mi()` still `b_` where brms uses
`bsp_`; and a written `theta` formula's reference component.

**Deferred by the user, 2026-09-17:** `frmtmb_control(vectorize = FALSE)`
on branch `wt-vectorize`, with the census on `wt-vecshape`. Vectorizing
is slower on every random-effects model, 1.36x to 1,116x, because it
breaks the sparsity of the Laplace sparse Hessian tape. Only the sampler
tape could gain, and only for gaussian-like models, after three
reshapings whose first one can break custom likelihoods. Revisit only if
sampling speed becomes a priority.

**Bin 2 of brms's suite is NOT ported and the audit recommends against
it** (`dev/brms-suite-audit.md` section 9). Quote any "passes brms's
suite" fraction against 823, bins 1 and 2, never 2,011. Today bin 1 is
192 of 494.

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
