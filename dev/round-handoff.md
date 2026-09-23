# Handing a round to a new session

Written 2026-09-22, at the 0.61.0 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has been lost SIX times, and the sixth has a
dated trigger (a hard kill of R processes).

Main carries items 2.6d (`predict()` is brms's predictive summary; the
old predictor is `frm_linpred()`) and 2.6f (brms's post-fit shapes and
names), lane `wt-shapes`; the silent wrong answers of lane `wt-adefects`;
and the RTMB 2.0 fallout. `wt-vectorize` and `wt-vecshape` stay
branches, deferred by the user and never merged.

Versions: frmtmb **0.61.0**, frmtmb.sample **0.9.0**, frmtmb.eam 0.10.0,
frmtmb.learn 0.6.0, frmtmb.latent 0.5.0, frmtmb.spline 0.7.0,
frmtmb.coupling 0.5.0, frmtmb.ode 0.6.0. Every extension floors on
frmtmb 0.61.0, because each reads `frm_linpred()` or the new exports.

**RTMB 2.0 is in the user library** (the user's install, 2026-09-21).
It changes no answer: the objective agrees with 1.9 to 2 ulp. It
exports `atan2`, now in the `nl()` shadow set. RTMB 1.9 remains the
Windows R 4.6 binary on CRAN, so the transparency block skips there.

### Verified at the release commit

- Suite: 269 files, 15,919 assertions, 0 fail, 0 error; every drop
  against 0.60.0 explained in `dev/suite-baseline.md`.
- Gated: 38 of 38 files, 3,130 assertions, 0 fail.
- Scale: 7 of 7.
- `R CMD check --as-cran`: 8 of 8, no warnings or errors after one
  merge-artifact Rd fix in frmtmb.sample; 5 packages carry the
  environmental "V8 unavailable" NOTE.
- Ported brms bin 1: 243 of 494 (was 231), ledger rebuilt on the merged
  build.

**Consolidation found four defects the lane's three review rounds did
not**, all now fixed and pinned (`dev/shapes-findings.md` section 13):
`summary()` stopped on every REML and profile fit; `predict()` under
REML dropped the fixed-effect uncertainty (1.0147 against 1.1611);
the adefects forward of `allow_new_levels` on draws was bit-identical
to `re_formula = NA`; and two `skip_if_not()` guards in
`test-unpinned-seams.R` skipped a measurement in silence after the
`vcov()` rename. The lesson for a review: no harness in the lane or
review touched a REML fit, and a gate the lane never set
(`FRMTMB_BRMS_FIT_TESTS`) hid 25 errors until round 3. Ask every lane
to run the release tiers, not its own subset.

## What is next, in order

**One lane, next: `predict()` and `fitted()` carry the group effects'
uncertainty at a level the fit saw** (user decision, 2026-09-22: match
brms), with the partial-`re_formula` defect below, because a partial
`re_formula` must add only its own terms' variance. Draw each
replicate's group effects jointly from their conditional law, not a
per-row variance, or `summary = FALSE` is wrong across rows of one group
(the shared-new-level-draw lesson). The lane states which coverage the
interval claims and measures that one (`dev/test-backlog.md`).

**Silent first.**
- A partial `re_formula` is accepted and not honored, with nothing said
  (`dev/adefects-findings.md` section 11, item 7), and
  `re_formula = ~(1 | nosuch)` is treated as `NULL` silently
  (`dev/test-backlog.md`).
- `residuals()` on an ordinal fit still answers where brms refuses.

**Then loud ones.** `pp_check(type = "*_grouped")` is broken for every
grouped type, and `type = "violin"`; `fitted()` on a multivariate fit
refuses while `predict()` answers; on draws, an unseen level without
the flag gets a hint that leads to a refusal (`dev/test-backlog.md`).

**Open for the user:** new levels on draws are refused rather than
drawn per posterior draw.

**Filed during 2.6c and not fixed**, details in each lane's findings:
two prior-scope differences from brms (a `sd` prior with a group and no
dpar reaches `phi`; mv `sd group = g` without `resp` is accepted); two
frmtmb.sample default priors (rescor flat where brms uses `lkj(1)`, and
the offset intercept location); `mi()` still `b_` where brms uses
`bsp_`; and a written `theta` formula's reference component.

**Small follow-ups:** a sentence in `?frm` that frmtmb's REML and
mgcv's differ on a location-scale smooth, and a test measuring the gap
(`dev/test-backlog.md`).

**Deferred by the user:** `frmtmb_control(vectorize = FALSE)`
(2026-09-17; branches `wt-vectorize`, `wt-vecshape`): slower on every
random-effects model, 1.36x to 1,116x. EM-seeded starts for the
latent-discrete families (2026-09-21; `dev/extension-gaps-plan.md`).

**Bin 2 of brms's suite is NOT ported and the audit recommends against
it** (`dev/brms-suite-audit.md` section 9). Quote any "passes brms's
suite" fraction against 823, bins 1 and 2, never 2,011. Today bin 1 is
243 of 494.

## Decisions the user made on 2026-09-23, for the four-lane round

- **The tiebreaker, refined:** match brms, UNLESS it is clearly obvious
  that brms should be doing it the other way.
- **`re_formula = ~(1 | nosuch)` stays REFUSED**, against brms, which drops
  an unmatched term silently: a user who names a term meant something by
  it. The same reading governs `~(1 | a/b)` on a fit that has `a` but not
  `a:b`, where brms keeps the half it recognizes.
- **A partial `re_formula` stays REFUSED on a fit with a factor-smooth
  term** (user, 2026-09-23). `NA` drops the smooth and `NULL` keeps it, and
  a partial formula cannot say which, so frmtmb refuses rather than guess;
  brms would keep the smooth. Same principle as the `~(1 | nosuch)`
  decision: naming some terms and not others must not silently answer a
  different question. This is the second departure from brms in the same
  argument, and both go the same way.
- **Prefix-less `b`, `Intercept` and `sigma` priors in a multivariate model
  are refused**, as brms refuses them. This is wider than the `sd` half the
  correctness lane took, and it breaks more call sites; it is a lane of its
  own.
- **`param_uncertainty` is renamed `propagate_error`** (user, 2026-09-23),
  defaulting to `TRUE`, so `predict(fit, propagate_error = FALSE)` holds
  the parameters AND the group effects at their estimates. The name states
  the axis: `re_formula` chooses WHICH terms are in the prediction, and
  this chooses whether the error in the estimates is propagated into the
  interval. Neither can express the other, which is why both exist: only
  this one can say "include the group effect but treat it as known", and
  only `re_formula` can say "predict for an average group". brms needs no
  such argument because its draws always carry both. Spelled out, not
  `propagate_err`: 191 documented argument names in core contain no
  clipped word. `plug_in` and `incl_uncertainty` were considered and
  rejected, the first as jargon that reads backwards, the second because
  it reads as "which things are included". The old name goes outright
  rather than deprecated. It is a `predict()` argument only; `fitted()`
  and `frm_linpred()` do not take it.
- **Build `simulate(newdata = )`** and let `pp_check()` pass it through, so
  `pp_check(newdata = )` answers as brms does instead of being refused. Two
  ported rows go to passes. That lane owns `simulate()`'s argument surface,
  and it BRINGS `simulate()` INTO LINE with `predict()` (user decision,
  2026-09-23): `~1` means no group effects, and a partial formula is honored.
- **Versions at this consolidation:** frmtmb 0.62.0, frmtmb.sample 0.10.0,
  a bump on every extension that changed, and every extension floor moved
  to frmtmb 0.62.0.

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
- `variables()` keeps frmtmb's order (2026-09-17); frmtmb objects carry
  no brms class, so `hypothesis()` output is `frmtmb_hypothesis` alone.
- `REML = TRUE` integrates the `mu` coefficients only; distributional
  coefficients stay outer (the double-GLM REML, `cd5bb83`). Integrating
  them would make a gaussian sigma the ML estimate again, and the inner
  problem can be unbounded.

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

`wt-shapes` and `wt-adefects` are merged and removed, with their evidence
committed on main under `dev/`. Create fresh worktrees off the current
main.
