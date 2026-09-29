# Handing a round to a new session

Written 2026-09-24 at the 0.63.0 release, rewritten 2026-09-28 at the
0.64.0 brms-parity release and 2026-09-29 at the 0.65.0 release. Read
this, then `dev/extension-gaps-plan.md`, then `dev/organizer-rules.md`
and `dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has been lost EIGHT times, and the last three
losses each followed processes being killed, not low disk.

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

**The docs cutover is the owner's next step.** The workflow and the
104 redirects (as a `redirects:` block in `_pkgdown.yml`) shipped with
0.65.0; switching Pages to the workflow builder and untracking `docs/`
follow the runbook in `dev/docsci-findings.md`, steps 3 to 7. A
released site beside a development site is recorded there as later
work, for when a tagged release exists that people install.

**Filed at 0.65.0, in `dev/test-backlog.md`:** `pp_check()`'s four
`loo_*` types fail on every model (no PSIS object is built);
`posterior_predict()` and `posterior_epred()` on `laplace = TRUE` draws
return NaN silently; `frm_sample(laplace = TRUE)` with no random effect
dies with an internal error; `fitted(newdata = )` on a multivariate
ordinal fit errors; `variables()` on a draws object gives internal names
for ordinal fits; `mo(m) + m` is unidentified and nothing reports it;
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

None should be live once 0.65.0 is committed. The round used seven
lane worktrees off `9b4bb650` (`wt-gradcheck`, `wt-csfactor`,
`wt-resmooth`, `wt-thresrefit`, `wt-arcovsample`, `wt-records`,
`wt-docsci`) and one integration worktree, `wt-release`, where the
seven were combined and verified. Remove a worktree and prune its
branch only after its work is merged and its evidence is committed on
main under `dev/`.
