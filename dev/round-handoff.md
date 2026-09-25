# Handing a round to a new session

Written 2026-09-24, at the 0.63.0 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has been lost EIGHT times, and the last three
losses each followed processes being killed, not low disk.

Five lanes merged:

- `wt-simnewdata`: `simulate(newdata = )` and `pp_check(newdata = )`,
  and `re_formula` in `simulate()` read as in `predict()`.
- `wt-mvprior`: a prior whose target brms refuses is refused, not
  broadcast, in multivariate, nonlinear, categorical, mixture and
  several-location extension families.
- `wt-predfix`: the filed `predict()`/`fitted()` defects, `cs()` in
  `predict()`, and the `autoscale = NULL` default.
- `wt-phase3a`: items 3.1, 3.2 and 3.6 (`frm_ode_records()`, coupling
  ingestion, `rp_floored()` on groups with no events).
- `wt-phase3b`: items 3.3, 3.4 and 3.5 (learn sessions, `wiener()`
  censoring, the `wiener()` contaminant).

Versions: frmtmb **0.63.0**; frmtmb.sample **0.11.0**, frmtmb.eam
**0.11.0**, frmtmb.latent **0.6.0** and frmtmb.learn **0.7.0**, which
floor on frmtmb 0.63.0 because each needs something only it has;
frmtmb.coupling **0.6.0**, frmtmb.spline **0.8.0** and frmtmb.ode
**0.7.0**, whose floors stay at 0.61.0 because they call nothing new.

### Verified at the release commit

- Suite: 292 files, 17,167 assertions, 0 fail, 0 error. 13 files are
  new; one fell by one assertion, by design (`dev/suite-baseline.md`).
- Gated: 39 of 39 files, 3,288 assertions, 0 fail, 0 files with a skip.
- Scale: 7 of 7.
- `R CMD check --as-cran`: 8 of 8. Five with 1 NOTE, the environmental
  V8 one on the HTML manual; three OK. No examples-timing NOTE, on a
  quiet machine. In-check tests FAIL 0 everywhere; core PASS 11,115,
  kept as `dev/release/frmtmb-testthat.Rout`.
- Ported brms bin 1: 252 of 494, from 250.

**Silent wrong answers found and fixed this round**, most by adversarial
review rather than by any tier:

- `simulate(re_formula = NA)` redrew population smooths, so
  `pp_check()` was wrong on every smooth model by default.
- `predict()` and draws `posterior_predict()` dropped `cs()` entirely.
- Priors: `b` without `nlpar` on `a ~ 1 + z` went to `a_z` only;
  categorical and mixture models broadcast `b` and `Intercept` across
  every `mu`; `class = "b"` left `cs()` coefficients out.
- `wiener()`'s contaminant gradient was a staircase (in-lane, never
  released).
- The new autoscale default hid separation and turned returned fits into
  errors (in-lane, caught by review, never released).

## What is next, in order

**Silent wrong answers filed and not fixed**, highest first:

1. `cs()` on a factor is fitted on the factor's integer codes, and a
   newdata factor is re-coded from its own levels. brms builds treatment
   dummies. `dev/test-backlog.md`, wt-predfix section.
2. `predict(re_formula = NA)`, and so `simulate(NA)`, drops
   `s(g, bs = "re")`, `fs` smooths and `t2` smooths with an `re` margin;
   brms keeps every smooth. User decision 2026-09-24: match brms.

**Filed defects and gaps**, in `dev/test-backlog.md`: autoscale on
covariance structures other than `us`, `diag` and Student-t; `frm_sample()`
on a one-parameter model; mixture sampling defaults for `sigma` and
`theta`; six prior spellings where frmtmb and brms disagree on accept or
refuse; missing `default_prior()` rows; `bf + bf + bf`, `cumulative` in a
multivariate model, `me()`, `0 + Intercept`, and `student` with
`rescor`; a negative-hazard penalty for `royston_parmar()`;
`frm_curve()` refusing on `gamma1 ~ x`; two core seams for frmtmb.eam
(a log-difference `lcdf` slot, a refusal of NA `dec()` on censored
rows); `lba()` and `gddm()` censoring and `lba()`/`rdm()` contaminant,
not built; base under-coverage of `log sd(log bs | s)`.

**Upstream reports, drafted and NOT filed**: RTMB `log_pnorm_both`'s
derivative (`dev/phase3b-rtmb-report-pnorm.md`); TMB `TanhOp::reverse`
(`dev/phase3b-rtmb-report-tanh.md`); drmTMB's REML and `beta_sigma`
(`dev/drmtmb-findings.md`); TMB's macOS binary and OpenMP.

**From the drmTMB comparison** (`dev/drmtmb-findings.md`):
`hurdle_negbinomial`, `zero_one_inflated_beta`, a negbinomial CDF for
`trunc()`, `gr(g, by = f)`, `fcor()`, boundary-corrected
variance-component tests, heritability and ICC accessors.

**Phases 4 and 5 of `dev/extension-gaps-plan.md` are untouched.**

## Decisions the user made on 2026-09-24

- **Prior spellings brms refuses:** refuse `class = "b"` where the
  predictor has no slope, `rescor` with `resp`, and `resp = "y"` on a
  univariate model. Keep `cor` with `resp` and `Intercept` with `nlpar`
  as frmtmb extensions.
- **Several-location extension families** (`lca`, `hmm`, `lba`, `rdm`,
  `mixture_mvn`) need `dpar` on `b` and `Intercept` priors, frmtmb's own
  rule where brms has no such family (organizer's call, not objected to).
- **`frm_bootstrap()` keeps its whole-model default**, redrawing group
  effects and smooths; `re_formula = NULL` conditions on them. No new
  argument.
- **Residual correlation at newdata counts fitted time levels**, a
  documented departure from brms, whose position-based reading is not
  consistent under marginalization.
- **Random-effect smooths under `re_formula = NA` should match brms**
  (keep them). Filed, not built.
- **Keep reporting nonzero optimizer convergence codes.**
- **`autoscale` defaults to `NULL`** (BREAKING), engaging below sd 1e-3,
  or 0.05 for a random-slope column, and falling back to the plain fit
  when its pre-fit fails.
- **eam contaminant window:** `contaminant_range =`, or with a `trunc()`
  upper bound, the fastest response to the deadline; otherwise refused.
- **eam left and interval censoring** use the boundary `dec()` names.
- **`rp_floored()`** warns when every flagged row is censored and refuses
  on a flagged event row.

**Settled earlier, do not reopen:**
- The tiebreaker: match brms, UNLESS it is clearly obvious that brms
  should be doing it the other way.
- `re_formula = ~(1 | nosuch)` stays refused, and a partial `re_formula`
  stays refused on a fit with a factor-smooth term.
- `propagate_error` is the name, spelled out.
- The non-generic name collisions with brms stay, because `::` is
  sufficient in both load orders.
- A gratia older than 0.9.0, loaded before frmtmb.sample, stops it
  loading; the user accepted that.
- `inverse.gaussian` defaults to brms's `1/mu^2`.
- Duplicate priors on one slot, and duplicate group-level effects, are
  refused as brms refuses them.
- `frm_simulate(newparams =)` takes brms names only.
- `variables()` keeps frmtmb's order; frmtmb objects carry no brms class.
- `REML = TRUE` integrates the `mu` coefficients only.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes. A reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a shared reference build
of the base commit. Punch rounds go back to the worker, and a reviewer
re-checks.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Two punch rounds is the cap, a third only for a blocker. The previous
release's library served as the reference build again this round.

Consolidation: commit each lane on its branch, merge, set versions,
roxygenise, THEN install (that order), verify the tree is unchanged by
roxygenise, run the tiers, commit the release and the docs separately,
remove the worktrees, prune branches, regenerate
`dev/suite-baseline.tsv`. Do not install into the release library while
any lane still uses it as its base build.

**Run `R CMD check` on a QUIET machine.** Its examples-timing NOTE
measures load here; it did not appear at this release, run alone.

**Memory is shared.** At most 3 R fitting processes per lane, each
started with 5 GB free (`dev/lane-rules.md`). One lane running 11 at
once crashed the machine this round.

## What the user has settled

- **The user pushes.** No session pushes, and no git operation without
  their say so.
- **Nobody uses this package yet.** Break backward compatibility
  freely: ship the refusal, bump, and say plainly in NEWS what stops
  working.
- **StanHeaders is pinned OUTSIDE `%LOCALAPPDATA%`**, in `pinlib` at
  2.32.10. It survived every library loss. The user library keeps
  2.39.1.
- **The R user library stays under `%LOCALAPPDATA%`** by the user's
  decision, knowing it will recur. RTMB 2.0 comes from r-universe after
  a restore, because CRAN's Windows binary is 1.9.
- Core is the user's lane except where they ask otherwise.

## What this round is evidence for

Five lanes, each with at least one adversarial review and one or two
punch rounds; one machine crash and one power loss.

**Review still finds silent wrong answers at about one per lane.** Most
of those listed above were found by a reviewer, not by a tier, and
several predate this round. Until that rate falls, a result should not
be trusted without a check against brms or an exact reference.

**A shared seed is a hidden replicate count.** Phase 3b reported a
drift-intercept interval covering 84.8 percent over 231 fits. The arms
drew their random effects right after the same `set.seed()`, so the 231
fits were 80 seeds, and an interval that knew every true drift did no
better. Fresh seeds on the released build cover 65 of 70. Coverage
claims now use seeds independent across arms, or say they are shared.

**A fix can be worse than the defect it fixes.** The autoscale default
closed a 62-unit silent stall and opened a path where a separated fit
reported success with no warning. The rule that closed it is general:
a default must never be LESS diagnostic than the setting it replaces.

**Interrupted runs are a standing condition, not an accident.** A crash
and a battery shutdown each cut every lane mid-tier this round. What
held: treat every log after the cut as void, check each saved result
reads back, verify each lane library against its source before trusting
it.

## Worktrees

`wt-simnewdata`, `wt-mvprior`, `wt-predfix`, `wt-phase3a` and
`wt-phase3b` are merged and removed, with their evidence committed on
main under `dev/`. Create fresh worktrees off the current main.
