# Handing a round to a new session

Written 2026-09-24 at the 0.63.0 release, rewritten 2026-09-28 at the
0.64.0 brms-parity release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: the library has been lost EIGHT times, and the last three
losses each followed processes being killed, not low disk.

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

**Six lanes are running now**, each in its own worktree:

1. `wt-gradcheck`: the convergence check fires on correct fits. "Large
   maximum absolute gradient" appears on ordinal and multivariate fits
   whose likelihood identities hold to 1e-12, and on fits with an
   active bound. The 1e-3 threshold is absolute.
2. `wt-csfactor`: `cs()` on a factor is fitted and predicted on the
   factor's INTEGER CODES, a silent wrong answer; and `y ~ x + cs(x)`
   is not identified and is not refused.
3. `wt-resmooth`: `predict(re_formula = NA)`, and so `simulate(NA)`,
   drops re-indexed smooths. User decision 2026-09-24: match brms.
4. `wt-thresrefit`: refits recount the ordinal thresholds, so
   `frm_bootstrap()` can lose one; and `draw_prior_entry()` on an
   unordered threshold vector may copy one draw into every threshold.
5. `wt-arcovsample`: frmtmb.sample `log_lik()` and `loo()` for
   `cov = FALSE` ARMA.
6. `wt-records`: the documentation and repository records of this
   round, including the pkgcheck failures on `R/me.R` and
   `R/ad-env.R`.

**Then, in order:**

- **Priors on `me()` hyperparameters.** Classes `meanme`, `sdme` and
  `corme` are refused, not implemented (lane me).
- **`gr(g, by = f, cov = A)`.** Refused by name today, because brms
  correlates the by-levels through `A`, which is not a by-split
  (lane grby).
- **Latent-residual AR for non-gaussian families**, which is brms's
  other autocorrelation form (lane arcov).
- **The multivariate family-fill decision.** In
  `bf(o ~ x) + cumulative() + bf(y ~ x) + gaussian()`, frmtmb fills
  only the responses that have no family, so `o` stays ordinal; brms
  gives the last family to every response. frmtmb behaved this way
  before the round and the behavior is now documented. The user has to
  decide whether it should follow brms.
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

Six are live, all off the 0.64.0 release commit `ec0e6d5f`:
`wt-gradcheck`, `wt-csfactor`, `wt-resmooth`, `wt-thresrefit`,
`wt-arcovsample` and `wt-records`, each on the branch of its own name
in `../frmtmb-wt-<name>`. What each one carries is listed under "What
is next". Remove a worktree and prune its branch only after its work
is merged and its evidence is committed on main under `dev/`.
