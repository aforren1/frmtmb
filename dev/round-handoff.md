# Handing a round to a new session

Written 2026-09-24 at the 0.63.0 release, rewritten 2026-09-28 at the
0.64.0 brms-parity release, 2026-09-29 at the 0.65.0 and 0.66.0
releases and 2026-09-30 at the 0.67.0 release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has been lost EIGHT times, and the last three
losses each followed processes being killed, not low disk.

## The 0.67.0 round, 2026-09-30

`dev/round-20260930.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-09-30-<lane>.md`. Three lanes of brms parity merged:

- `ordinal`: brms's `disc` on every ordinal family,
  `threshold = "equidistant"` and `"sum_to_zero"`, and every brms link
  on `acat()`, with brms's draw names in frmtmb.sample.
- `ceplot`: `plot()` of conditional effects and of a hypothesis
  returns plot objects and takes brms's arguments; brms's effect check;
  new levels on crossed and `mm(by = )` terms in
  `conditional_effects()`; `fitted(sample_new_levels = "old_levels")`
  with one choice per grouping factor and per call; `parnames()` on a
  fit; `nsamples(incl_warmup = TRUE)`; brms's order in
  `posterior_samples(pars = )`.
- `formrobust`: addition-term expressions, `update()` that keeps the
  parameter formulas, bernoulli on any two values, `ar()` and the other
  autocorrelation functions, `autocor()`, `drop_unused_levels`, offsets
  in the grids of `conditional_effects()` and `emmeans()`, and the
  ARMA fill of a missing response on newdata.

Versions: frmtmb **0.67.0**; frmtmb.sample **0.15.0**, floor frmtmb
0.67.0 (it imports `parnames()` and calls `ord_delta_info()`,
`arma_cond_fill_dpars()`, `arma_cond_fill_epred()` and
`response_codes_newdata()`). The other six extensions are unchanged.
The release library is `C:/Users/adf44/source/r/rellib-r5`; `rellib-r4`
keeps the 0.66.0 build. In brms's own ported suite, bin 1 passes
387 of 494 (368 at 0.66.0).

### Decisions the user made on 2026-09-30, for this round

- **The ordinal default of `conditional_effects()` stays the
  per-category display**, a deliberate divergence from brms. brms's
  default is the expected category number, and its own warning calls
  that display likely invalid for ordinal families and asks for
  `categorical = TRUE`, which is frmtmb's default. brms's warning comes
  with `categorical = FALSE`. Ledger row `brmsfit-methods:217` is a
  divergence.
- **A family-default fixed `disc` is hidden, as brms hides it.** A
  `disc` held at 1 shows nowhere: not on the `Links:` line, not as a
  `Fixed dpar:`, not in `fixef(flatten = TRUE)` or `coef()`. A `disc`
  the user fixes at another value still shows.
- **Addition terms read data variables only**, as brms requires. A
  variable found only in the formula environment is refused by name,
  because predictions on newdata and refits read the term again, where
  an outside value can be missing or changed. Predictor formulas and
  offsets keep R's rule of looking in the formula environment.
- **Filed as an idea, not built: snapshot outside constants at fit
  time**, so that a predictor or offset that reads a value outside the
  data reads the same value in a later prediction or refit.
- **bernoulli codes any two values by level order, as brms does**, and
  warns when both values lie strictly between 0 and 1, since that looks
  like a proportion (lme4#682). Effect coding such as -0.5 and 0.5 does
  not warn.
- **The ARMA fill of a missing response**: frmtmb.sample's draws sample
  the fill, as brms does, and the maximum likelihood `fitted()` fills
  with the expected value, a deliberate divergence, since a fit has no
  draws to carry the fill through.

### What the 0.67.0 round is evidence for

**A worker's report can fail to arrive; read its findings file.** The
findings file is the record the lane is judged on. When a report does
not reach the coordinator, the lane's state is still in
`dev/<lane>-findings.md`, and the consolidation read each lane's
intent there: for every conflict, and for the ledger rows each lane
said it flipped.

**A small punch round can reverse the direction of a change.** Lane
formrobust's first `emmeans()` offset change went one way, its review
found brms does the other, and the punch round turned it around. Run
brms on the case, rather than only reading its code: the reversal came
from a run (`dev/formrobust-rev-log/brms-emm.txt`), not from the
source.

**A merge of two clean lanes can make a wrong answer that neither
lane has.** formrobust put offset variables in the model frame for the
grids; ceplot's effect check read every variable of the terms. Together
they accepted `effects = "time"` for a variable only an offset reads,
which brms refuses, while each lane alone stopped or answered for
another reason. Run the interplay of two lanes that touch one function
on the merged tree, with brms beside it (`dev/rel067-effoffset.R`).

**A lane's record of a ledger row can be wrong in the direction it
expects.** formrobust recorded `brmsfit-methods:955` as a class-name
divergence because `update()` "refits". On the lane's own library the
refit stops at its starting values. Read the message the merged build
records for every row whose message moved
(`dev/rel067-msgdiff.R`), not only the rows a builder stops on.

## The 0.66.0 round, 2026-09-29

`dev/round-20260929b.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-09-29-<lane>.md`. Six lanes of brms parity merged:

- `sampfix`: frmtmb.sample reads `frm_sample(laplace = TRUE)` draws in
  their own layout (they returned NaN or finite wrong numbers) and
  refuses by name what reads the integrated values; `pp_check()`'s
  `loo_*` types; `stanfit = NULL` draws; brms's names for ordinal
  draws; one-parameter models.
- `fams2`: `xbeta()`, `zero_inflated_beta_binomial()`,
  `hurdle_cumulative()` and `cse()`.
- `aterms2`: `subset()`, `index()`, `mi(x, idx = )`, `rate()` and
  `cat()`; `nobs()` is brms's.
- `formula2`: dpar equations (`sigma1 = "sigma2"`), `cmc`, `y ~ .`, a
  list of families, and `update()` with a `bf()` delta.
- `defects`: 38 rows of the ported brms suite, among them two silent
  wrong answers (`fitted(scale = "linear")` with `cs()`, residuals of a
  missing `mi()` response).
- `postfit2`: `conditional_smooths()`, `make_conditions()`,
  `update_adterms()`, `posterior_average()`, `conditional_effects()`'s
  `surface`, `too_far`, `select_points` and `spaghetti`, a Wald band
  for a nonlinear predictor, and brms's per-row rule for which group a
  display row reads under `re_formula = NULL` (a silent wrong answer
  found by lane sampfix). `frm_bootstrap()` no longer redraws smooths.

Versions: frmtmb **0.66.0**; frmtmb.sample **0.14.0**, floor frmtmb
0.66.0 (new sampling-API exports, the ordinal inverse maps, equated
parameters); frmtmb.eam **0.11.2**, floor 0.66.0 (its `test-family.R`
reads `frmtmb:::row_aterms`; no code change); frmtmb.learn **0.7.1**
and frmtmb.ode **0.7.1** (tests only: `test-stan-identity.R` and the
scale tier's `test-scale.R` let no warning escape). frmtmb.coupling,
frmtmb.latent and frmtmb.spline are unchanged. In brms's own ported
suite, bin 1
passes 368 of 494 (306 on lane defects alone, 275 at 0.65.0).

### Decisions the user made on 2026-09-29, for this round

- **`default_prior()` and `validate_prior()` skip the response check**,
  as brms's do: a `Beta()` model on a response outside (0, 1) gets its
  table, and `frm()` still refuses the fit. An ordinal response the
  family cannot read is still refused, as brms's `extract_nthres()`
  refuses it.
- **`update(fit, data = )` is kept** and refits on the new data, as
  `stats::update()`, lme4 and glmmTMB read it. brms refuses that
  spelling. Recorded as a divergence (`brmsfit-methods:927`).
- **Laplace conditional draws stay refused.** `frm_sample(laplace =
  TRUE)` samples an approximate posterior whose use is checking the
  approximation; a call that reads the integrated values (group-level
  coefficients, a smooth's coefficients, `mi()` values) is refused by
  name, and the values are not filled in from `N(b_hat, H^-1)`.
- **`frm_bootstrap()` never redraws a smooth, `gp()` or `hsgp()`
  term**, one indexed by a grouping factor included, and
  `conditional_effects(band = "boot")` inherits the rule. This
  withdraws the whole-model default of 2026-09-24 (below).
- **Upstream bugs are listed in `dev/upstream-bugs.md`** for the user
  to file. No session files them.
- **Nobody uses the package, so break compatibility freely.** This
  round did: `nobs()`, `frm_bootstrap()`, `fitted(scale = "linear")`
  on `cs()` fits, the prior-table spellings and the ordinal draw names
  all changed, and NEWS says so under "Breaking changes".

### What the 0.66.0 round is evidence for

**Run every extension suite that can reach a changed path.** Lane
aterms2 added `subset()` and `index()` to core's compatibility table
and ran core and frmtmb.sample twice. The change broke 5 assertions in
frmtmb.eam's `test-family.R`, which the reviewer found only by running
all eight suites. `dev/lane-rules.md` now says so.

**A `test_file()` runner must attach the package.** Lane runners that
called `test_file(package = p)` without `library(p)` reported failures
in `test-conditions.R`, `test-data2.R` and `test-id-kron.R` that are
the runner's. Measured on the release build
(`dev/round-20260929b.md`): unattached, `test-conditions.R` has 3
failures and 3 errors and the other two 1 error each; attached, all
three pass. `dev/lane-rules.md` now says so.

**A small punch round can add large code, so give it a final check.**
fams2's second punch round added a "some rows" rule to the `kappa`
warning, and the final check found it warned on well-identified fits
(K1). postfit2's punch rounds replaced the whole new-level
construction of `conditional_effects()` with a per-row, per-term plan
(`R/ce-levels.R`, new), and each re-check found a new blocker in it.
Neither would have been seen without a review after the last punch
round.

**A clean textual merge can call a function another lane removed.**
sampfix's laplace probe in frmtmb.sample's `conditional_effects()`
called `ce_draw_new_levels()`, which postfit2 removed, and the two
lanes' hunks merged without a conflict. The probe is now written on
postfit2's `ce_plan_eval()`. Grep for every removed export after a
merge; the suite finds the rest.

## The 0.65.0 round, 2026-09-29

`dev/round-20260929.md` is this round's record, with the verification
of the release tree. Each lane has `dev/<lane>-findings.md` and
`dev/reviews/2026-09-29-<lane>.md`. Seven lanes merged:

- `gradcheck`: the convergence warning no longer fires on correct fits
  (289 of 720 before, 0 after); bound-aware Newton decrement.
- `csfactor`: `cs()` on a factor builds brms's treatment dummies;
  `x + cs(x)` and `mo(m) + cs(m)` refused by a rank test.
- `resmooth`: `re_formula = NA` keeps every smooth; ordinal `Est.Error`
  now carries the smooth's coefficients.
- `thresrefit`: `influence()` refits keep the fitted threshold layout
  and report failed refits.
- `arcovsample`: frmtmb.sample `log_lik()` and `loo()` for
  `cov = FALSE` ARMA; the Student-t rescor row density at large `nu`.
- `records`: the pkgcheck failures and this round's records.
- `docsci`: `.github/workflows/pkgdown.yaml` and
  `dev/release/build-docs.R`. Its deploy job skips until the owner
  switches Pages to the workflow builder; the cutover runbook is in
  `dev/docsci-findings.md`.

One more lane was opened DURING consolidation, because the release
check found frmtmb.spline's `curve-inference` vignette failing against
the new `re_formula = NA` rule:

- `splinecurve`: `allow_new_levels` through `frm_curve()`,
  `frm_curve_deriv()` and `frm_curve_feature()`, so the population
  curve of a factor-smooth model is read at an unseen level; a refusal
  of a difference between unseen `(1 | g)` levels, whose standard error
  would otherwise be exactly 0.

Versions: frmtmb **0.65.0**; frmtmb.sample **0.13.0**, floor frmtmb
0.65.0 (it calls `arma_cond_resp()` and `arma_cond_dpars()`);
frmtmb.spline **0.9.0**, floor 0.65.0 (the new `re_formula` rule and
`allow_new_levels`); frmtmb.eam **0.11.1** (a Suggests line). The
other four extensions are unchanged.

## The 0.64.0 round, 2026-09-25

`dev/parity-round-20260925.md` is the 0.64.0 round's own record. Each
lane's findings file has its validation, its numbers and its scripts.

Nine lanes merged:

- `me`: `me()` noise-free predictors and `set_mecor()`.
  `dev/me-findings.md`
- `thres`: `thres()` for the ordinal families.
  `dev/thres-findings.md`
- `grby`: `gr(g, by = )` and `mm(g1, g2, by = )`.
  `dev/grby-findings.md`
- `arcov`: `ar()`, `ma()` and `arma()` without `cov = TRUE`.
  `dev/arcov-findings.md`
- `mv`: three or more responses summed with `+`; `student()` with
  `rescor`; ordinal families in a multivariate model.
  `dev/mv-findings.md`
- `icpt0`: `0 + Intercept` and `bf(center = FALSE)`.
  `dev/icpt0-findings.md`
- `fams`: `hurdle_negbinomial()` and `zero_one_inflated_beta()`.
  `dev/fams-findings.md`
- `emm`: emmeans on nonlinear and multivariate fits.
  `dev/emm-findings.md`
- `sratio`: `sratio()` thresholds unordered, found by lane thres.
  `dev/sratio-findings.md`

Versions: frmtmb **0.64.0**; frmtmb.sample **0.12.0**, which floors on
frmtmb 0.64.0 because it needs `rescor_row_loglik()` and the new
`fit_extras(resp = )`. Every other extension is unchanged from the
0.63.0 release: frmtmb.eam **0.11.0**, frmtmb.latent **0.6.0** and
frmtmb.learn **0.7.0** keep their floor at frmtmb 0.63.0;
frmtmb.coupling **0.6.0**, frmtmb.spline **0.8.0** and frmtmb.ode
**0.7.0** keep theirs at 0.61.0.

### Verified at the release commit

Verification ran on a LINUX CONTAINER, not on the Windows box, which is
why `R CMD check` was `--no-manual` there: the container has no LaTeX.
The Windows run with the manual was done on 2026-09-28 on the fresh
install (`dev/release/check.log`; core and eam rerun in
`dev/release/check-rerun.log` after drmTMB and makeindex were
installed): core, frmtmb.eam, frmtmb.latent, frmtmb.ode and
frmtmb.spline at 1 NOTE, the environmental V8 one; frmtmb.coupling,
frmtmb.sample and frmtmb.learn OK. The ungated suite on the same tree
was 305 files, 17,893 assertions, 0 fail, 0 error
(`dev/suite-baseline.md`). The reference log-likelihoods of
`dev/machine-library.md` reproduce exactly on the container:
-295.602189818 and -332.137876329.

- Suite, core and all seven extensions, one file per process: 305
  files, 17,815 assertions. Two files fail, and both fail identically
  on the unmodified base build there: `test-pp-check-types.R` (4), from
  conda's bayesplot 1.15.0, and frmtmb.sample `test-reparam.R` (1).
- Gated: 54 files, 3,767 assertions, 0 fail, 0 error, with Stan
  compiling. `test-drmtmb-agreement.R` (131) and `test-fuzz.R` ran
  separately behind their own gates, 0 skip.
- Ported brms suite: the regenerated verdicts file and every generated
  test file are byte-identical to the committed ones. Bin 1 passes 275
  of 494, up from 252.
- `R CMD build`: OK, vignettes included.
- `R CMD check --as-cran --no-manual`: in-check tests 9,098 pass, 4
  fail, all four `test-pp-check-types.R`. The other warnings and notes
  are the container's: the locale, a missing `qpdf`, no CRAN index and
  the clock check. One warning was real, an undeclared `ordinal::` in
  `test-thres.R`; `ordinal` is now suggested.

The Windows machine got a FRESH R INSTALL on 2026-09-28, R 4.6.1, with
the user library restored from `dev/restore-cran-2026-09-28.txt` and
the manifest `dev/win-library-4.6-manifest-2026-09-28.csv` (412
packages, 402 of them from CRAN). Several Windows-side records were
stale after that and are corrected in this round: the StanHeaders pin
is gone (`dev/lane-rules.md`, `dev/machine-library.md`), and
`codemeta.json` is regenerated at 0.64.0.

**Defects the merge created, found and fixed on the branch.** A lane
can only test its own feature, so these appeared where two lanes met or
in test files a lane did not run: `ordinal_ncat()` read the first
response's family; `cs()` in a multivariate ordinal response had zero
coefficients; `lp$center` was read with `$`; two refusals shared one
message; a nonlinear body with `m[, 1]` stopped with 'argument "ei" is
missing'; nine manual port verdicts were stale. Both
`test-parity-integration.R` pins were seen failing on a build without
the fixes. The table with each cause is in
`dev/parity-round-20260925.md`.

## What is next, in order

**The docs cutover is DONE, 2026-09-29.** GitHub Pages builds from
`.github/workflows/pkgdown.yaml`: every push to main that touches the
documented sources rebuilds and deploys the eight sites, and `docs/` is
no longer tracked (it is in `.gitignore`; a local build still writes
there by default). The first deploy was of `c9d9126b`, and all 22 checks
of runbook step 5 passed on the live site. To roll back, re-track a
built `docs/` and switch Pages to `legacy` from `main:/docs` (runbook
step 7). A released site beside a development site is recorded in
`dev/docsci-findings.md` as later work, for when a tagged release
exists that people install.

**Filed at 0.66.0, in `dev/test-backlog.md`** under "Filed at the
0.66.0 release": what the six lanes left, led by a NaN the laplace
probe can pass, `predict(newdata = )` with NA responses under
`cov = FALSE` ARMA, and the bernoulli recoding of a two-valued
response.

**Filed at 0.65.0, in `dev/test-backlog.md`:** closed at 0.66.0 by
lane sampfix: `pp_check()`'s four `loo_*` types,
`posterior_predict()` and `posterior_epred()` on `laplace = TRUE`
draws, `frm_sample(laplace = TRUE)` with no random effect, and
`variables()` on the draws of an ordinal fit. Still open:
`fitted(newdata = )` on a multivariate ordinal fit errors;
`mo(m) + m` is unidentified and nothing reports it;
the finite-difference `Est.Error` at `newdata` off a `gp()`'s fitted
positions costs one evaluation pair per coefficient; standard errors
are not bound-aware at a constrained optimum.

**Then, in order:**

- **Priors on `me()` hyperparameters.** Classes `meanme`, `sdme` and
  `corme` are refused, not implemented (lane me).
- **`gr(g, by = f, cov = A)`.** Refused by name today, because brms
  correlates the by-levels through `A`, which is not a by-split
  (lane grby).
- **Latent-residual AR for non-gaussian families**, which is brms's
  other autocorrelation form (lane arcov).
- **Phases 4 and 5 of `dev/extension-gaps-plan.md` are untouched.**

**Also open, and filed in `dev/test-backlog.md`:** REML with `me()` in
`mu` is approximate and registered as conditional, and the same
argument applies to `mi()`, which is registered as working; autoscale
on covariance structures other than `us`, `diag` and Student-t;
`frm_sample()` on a one-parameter model; mixture sampling defaults for
`sigma` and `theta`; six prior spellings where frmtmb and brms disagree
on accept or refuse; missing `default_prior()` rows; a negative-hazard
penalty for `royston_parmar()`; `frm_curve()` refusing on `gamma1 ~ x`;
two core seams for frmtmb.eam (a log-difference `lcdf` slot, a refusal
of NA `dec()` on censored rows); `lba()` and `gddm()` censoring and
`lba()`/`rdm()` contaminant, not built; base under-coverage of
`log sd(log bs | s)`.

**Upstream reports, drafted and NOT filed**: RTMB `log_pnorm_both`'s
derivative (`dev/phase3b-rtmb-report-pnorm.md`); TMB `TanhOp::reverse`
(`dev/phase3b-rtmb-report-tanh.md`); drmTMB's REML and `beta_sigma`
(`dev/drmtmb-findings.md`); TMB's macOS binary and OpenMP. One more is
ready to draft: brms's `posterior_predict_hurdle_negbinomial()` does
not draw the zero-truncated negative binomial (lane fams).

**From the drmTMB comparison** (`dev/drmtmb-findings.md`), what is
left after this round: `fcor()`, boundary-corrected variance-component
tests, heritability and ICC accessors. `hurdle_negbinomial`,
`zero_one_inflated_beta` and `gr(g, by = f)` shipped at 0.64.0; a
negbinomial CDF for `trunc()` has not.

## Decisions the user made on 2026-09-29

- **Redirects** for the 104 retired core reference pages live as a
  `redirects:` block in `_pkgdown.yml`, not as static files.
- **The docs workflow has no `release:` trigger.** A push to main
  deploys; a tagged release does not redeploy an older tree. The Pages
  cutover's one-build gap is acceptable.
- **`y ~ x + cs(x)` stays refused**, in `frm()` and in `frm_sample()`
  alike. brms samples it and the prior sets the split between the global
  and category-specific parts, which the data cannot inform; `cs(x)`
  alone gives the same predictive distribution. A deliberate departure
  from brms.
- **A partial `re_formula` beside a factor smooth is ACCEPTED** and
  keeps the smooth, as brms does. This lifts the 0.63.0 refusal, whose
  reason (`NA` drops the smooth) stopped being true at 0.65.0. The
  refusal is gone entirely: its other documented case, a `car()` or
  `spde()` field, was never reachable, because a field carries the bar
  `(1 | loc)` and a formula names it like any group-level term.
- **A family added after several multivariate responses fills only the
  responses that have none.** brms gives it to every response, which
  silently overwrites an ordinal response's family; frmtmb keeps its
  own rule, a documented departure.

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
  argument. WITHDRAWN 2026-09-29 for smooths: `frm_bootstrap()` never
  redraws a smooth, `gp()` or `hsgp()` term (see the 0.66.0 round).
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
- `re_formula = ~(1 | nosuch)` stays refused. (The companion rule, that
  a partial `re_formula` stays refused on a fit with a factor-smooth
  term, was REVERSED by the user on 2026-09-29; see above.)
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
measures load here; it did not appear at the 0.63.0 release, run alone.

**Memory.** On the machine of 2026-09-28, which has 63 GB of RAM and a
fast disk, test runs need no process cap. Keep ONE TEST FILE PER R
PROCESS and run as many at once as the work needs. The history stays
because the failure mode is real on a smaller box: in the 0.63.0 round
one lane running 11 R fitting processes at once crashed the machine of
that time, and 34 of its fits had already failed with `std::bad_alloc`
(`dev/lane-rules.md`).

## What the user has settled

- **The user pushes.** No session pushes, and no git operation without
  their say so.
- **Nobody uses this package yet.** Break backward compatibility
  freely: ship the refusal, bump, and say plainly in NEWS what stops
  working.
- **Nothing pins StanHeaders.** The `pinlib` at 2.32.10 was retired on
  2026-09-17. rstan 2.32.7 compiles against the user library's
  StanHeaders 2.39.1 because
  `C:/Users/adf44/Documents/.R/Makevars.win` carries
  `CXX17FLAGS += -std=gnu++17`, and tmbstan 1.2.1 no longer samples a
  standard normal. That file is outside `%LOCALAPPDATA%`, so a library
  loss does not touch it (`dev/tmbstan121-findings.md`).
- **The R user library stays under `%LOCALAPPDATA%`** by the user's
  decision, knowing it will recur. RTMB 2.0 comes from r-universe after
  a restore, because CRAN's Windows binary is 1.9.
- Core is the user's lane except where they ask otherwise.

## What the 0.65.0 round is evidence for

Seven lanes, each with an adversarial review and one or two punch
rounds.

**The reviews found false sentences more often than false numbers.** In
four of the five code lanes the first review blocked on documentation:
a manual page that promised every short fit would warn, a help page
that said `cs()` is not available in a multivariate fit, a help page
that described which rows of an ARMA group are shifted, and a findings
file whose reasoning was wrong in both directions while its outcome was
right. Each was measured false by construction. Review the prose
against the code with the same care as the code.

**A guard failed closed again.** The withholding of the gradient
verdict on importance-corrected fits also withheld the clean line on
every such fit, well converged or not, where 0.64.0 printed it. The
worker had argued in a code comment that withholding was the safe side.
The rule stands: construct the case where the guarded thing is ABSENT,
and a default must never be less diagnostic than what it replaces.

**A count settles what a reason cannot.** A reviewer attributed a
28.9x slowdown to a batching routine returning NULL. The worker counted
batches and model evaluations and showed the routine never returned
NULL: it split per column position, which on a smooth is one batch per
coefficient. The reviewer withdrew its cause on the page.

**The `tail -f` trap bit the consolidation itself, twice.** A watcher
tailing the four tier logs made the PowerShell drivers' `Add-Content`
throw, exactly as `dev/lane-rules.md` records. Stopping the watcher did
NOT release the logs: a monitor that is stopped or expires leaves its
`tail` process running, and ten such orphans from the session's
earlier watchers were still holding files open, so the relaunched
gated tier ran all 39 files, logged none of them, and printed "GATED
ran 39 of 39". Kill `tail.exe` before trusting a log, and wait on the
driver PROCESS, never on the log.

**A lane that changes core behavior must run every extension's check,
not only their suites.** The resmooth lane updated frmtmb.spline's help
pages and `test-curve.R` for the new `re_formula = NA` rule and ran the
extension suites green, but never built frmtmb.spline, whose
`curve-inference` vignette drew a population curve at `re_formula = NA`
on a grid with no grouping column. The release check found it; the fix
needed a new `allow_new_levels` argument through the `frm_curve()`
family, a lane of its own (`dev/splinecurve-findings.md`).

**A refusal can outlive its reachable cases.** Lifting the factor-smooth
case of the partial-`re_formula` refusal left its other documented
case, a `car()` or `spde()` field, which a new guard test showed was
never reachable: both carry a `(1 | loc)` bar. The help page had said
so since 0.63.0. A guard test for a surviving branch is what found it.

**A fresh machine has no git identity.** Commits could not be made
until the user configured one, so the consolidation built and verified
the integrated tree in `../frmtmb-wt-release` first and committed
afterwards.

## What the 0.64.0 round is evidence for

Nine lanes, consolidated and verified on a Linux container.

**A lane cannot test the seam it shares with another lane.** Six
defects came out of the merge itself, not out of any lane, and two of
them were wrong answers rather than refusals. The consolidating session
found them by running the whole suite on the merged tree, which is the
only place they exist.

**Verification on one platform leaves a hole.** `--no-manual` on the
container skipped the manual sections, which is where an unescaped `%`
in Rd has bitten this project before, and the pkgcheck failures this
lane fixed were on a Windows-side check that the container run never
reached. A record round after a platform move is not optional.

## What the 0.63.0 round was evidence for

Five lanes, each with at least one adversarial review and one or two
punch rounds; one machine crash and one power loss.

**Review still finds silent wrong answers at about one per lane.** Most
of the defects that round listed were found by a reviewer, not by a
tier, and several predated the round. Until that rate falls, a result
should not be trusted without a check against brms or an exact
reference.

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
and a battery shutdown each cut every lane mid-tier that round. What
held: treat every log after the cut as void, check each saved result
reads back, verify each lane library against its source before trusting
it.

## Worktrees

None should be live once 0.67.0 is committed. The 0.67.0 round used
three lane worktrees off `57c25589` (`wt-ordinal`, `wt-ceplot`,
`wt-formrobust`) and one integration worktree, `wt-release` off
`57c25589`, where the three were combined and verified. Their private
libraries are `C:/Users/adf44/source/r/wt-<lane>-lib`; the release
library is `C:/Users/adf44/source/r/rellib-r5`, and `rellib-r4` holds
the 0.66.0 reference build. Remove a worktree and prune its branch only after its
work is merged and its evidence is committed on main under `dev/`.
