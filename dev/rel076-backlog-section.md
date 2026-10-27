## Filed at the 0.76.0 release (2026-10-27)

Three lanes off 70476dc3 (frmtmb 0.75.0), merged onto 70476dc3
(`dev/round-20261027.md`). What each lane, its reviewer or the
consolidation filed and did not fix, with the record to read. Nothing
was posted upstream.

### The next round's first lane: brms's `sigma` under `ar/ma/arma(cov = TRUE)`

- **MAJOR (a reporting divergence). Under a first-order `cov = TRUE`
  `ar()`, `ma()` or `arma()`, frmtmb reports the MARGINAL sd as
  `sigma`, where brms's `sigma` is the innovation sd.** The model is the
  same; only the name differs, by a factor `c = sqrt(1 - ar^2)`,
  `1 / sqrt(1 + ma^2)` or `sqrt((1 - ar^2) / (1 + ma^2 + 2 ar ma))`.
  Measured by the acprior reviewer (`dev/review-acprior-r1-report.R` in
  the lane's tree, `r1-report.txt`; the AR(1) 0.9 data of
  `dev/acprior-agree.R`, seed 20261027): brms's posterior mean of
  `sigma` is 0.432; frmtmb gives `summary(fit)` 1.08,
  `VarCorr(fit)$residual__$sd` 1.078, `sigma(fit)` 1.078,
  `posterior_summary(ds)` and `summary(ds)` 1.158, the draws' `sigma`
  column 1.158, and `variables(ds)` lists `sigma`. brms's draws carried
  to the marginal sd give 1.149, and frmtmb's carried to the innovation
  sd give 0.4324. 0.76.0 reads every prior on class `sigma` at the
  innovation sd, as brms does, so the priors already agree; the
  reporters do not. A fix touches `summary()`, `VarCorr()`, `sigma()`,
  the draws and `variables()`, the predictions and `posterior_predict()`
  (which read sigma), the ML standard errors (a delta-method row in the
  new coordinates), `log_lik()` (whose conditional reads the marginal
  covariance; the 0.76.0 consolidation test "ar(cov = TRUE) under the
  default priors" in frmtmb.sample's `test-loglik-autocor.R` pins that
  reading), `hypothesis()` names, and the "Divergence from brms" in
  `?frmtmb-autocor`. A larger API change; frmtmb has no users, so it is
  allowed, in a lane of its own. Record: `dev/acprior-findings.md`, Fix
  round 2, item 3; `dev/reviews/2026-10-27-acprior.md`, Fix round 1
  check, item 2. The coordinator asked for this to be the next round's
  FIRST lane.
- **`sigma ~ x` under a first-order `cov = TRUE` term: brms's
  `Intercept_sigma` is the log innovation sd; frmtmb's default and user
  priors on the sigma intercept read the log marginal sd.** brms passes
  `sigma = exp(Intercept_sigma + ...)` to `normal_time_het_lpdf()` with
  the AR correlation divided by `1 - ar^2`
  (`dev/review-acprior-r2-brmssigx.R`, `r2-brmssigx.txt`). So a prior on
  `class = "Intercept", dpar = "sigma"` differs from brms's by a shift of
  `log c`; the slopes (`class = "b", dpar = "sigma"`) are the same in
  both. The default and a user prior read the same scale, so the
  sampling table still round-trips. The fix is small: read the density
  at `x + log c`, with no Jacobian term (a shift). Covers the DEFAULT
  and the USER prior both. `?set_prior` (Residual correlation) states
  the difference since 0.76.0. Record: `dev/reviews/2026-10-27-acprior.md`,
  Fix round 2 check, consolidation note 1.

### mo()

- **Two-category `mo()`.** brms 2.23.0 fits `mo(x)` on two values
  (`Jmo = 1`, codes 0 and 1, a simplex of one weight fixed at 1, a
  `simo mox1` row with `dirichlet(1)` in `default_prior()`); one value
  gives `Jmo = 0`. frmtmb refuses both ("mo() needs at least 3 ordered
  categories; mo(x) has 2"). A build needs a zero-length simplex through
  frame assembly, `mo_simplex()` and `mo_coords()` (`mo_sphere_basis(1)`
  calls `contr.helmert(1)`, which stops), `mo_search()`, the Dirichlet
  prior, `se_mo_gaps()`, `lp_delta_mo()`, and frmtmb.sample's softmax
  chart, its `dirichlet(1)` density and `draws_mo_cols()`. The values 1
  and 2, which 0.75.0 fitted with an unidentified extra cell, are
  refused since 0.76.0. Record: `dev/mocode-findings.md`, Fix round 1,
  item 7; `dev/reviews/2026-10-27-mocode.md`, finding 7.
- **`fixef()` and `variables()` of saved 0.75.0 `frm_sample()` draws**
  of an integer `mo()` variable that does not start at 0 are not
  refused; every computed quantity on such draws is. frmtmb.sample may
  not call `frmtmb:::mo_check_fit()` (guardrail 3), so the refusal needs
  a new exported seam, its docs and a test. Record:
  `dev/mocode-findings.md`, Fix round 2, Nits.
- **Gap-aware marginaleffects standard errors.** On a fit with a `mo()`
  category that no row has, every marginaleffects standard error is
  `NaN`, with one warning per fit, also at the observed categories;
  frmtmb's own routes give `NaN` only at the gap. marginaleffects forms
  `J V J'` itself, so no covariance can be `NaN` only where a gap code
  loads. A later improvement: pass the identified reparametrization (the
  observed categories' cumulative steps) through `get_coef()` on a gap
  fit. Record: `dev/reviews/2026-10-27-mocode.md`, Fix round 2 check,
  item 2.
- **`predictions(re_formula =)` through marginaleffects** is refused by
  marginaleffects ("not known to be supported ... re_formula"); not a
  `mo()` matter, pre-existing. Record: the same review, item 1.
- emmeans refuses a `mo()` variable held at its mean, as brms does.
  Holding it at `round(median())` instead is a one-line change in
  `emm_check_mo_grid()` if the user prefers convenience (a user
  decision). Record: `dev/mocode-findings.md`, "What needs the user
  after fix round 1".
- `attr(x, "min")` is not read; brms reads it for new data (brms-40,
  `dev/upstream-bugs.md`, is brms's own inconsistency there). Record:
  `dev/mocode-findings.md`, "Where frmtmb still differs from brms".
- `variables()` of an ML fit has no `simo_` entries (filed at 0.75.0,
  still open).

### Residual correlation priors

- **The random-effect and `rescor` double prior slot.** Class `sd` or
  `cor` and `class = "theta", coef = "theta_k"`, and `class = "rescor"`
  and `theta thetar_1`, are two spellings of one parameter. Measured by
  the acprior reviewer (`dev/review-acprior-re.R`, block `(1 + x + z |
  g)`): in `frm()`, `cor lkj(2)` then `theta theta_4 normal(0, 1)`
  leaves only the theta entry on position 4, and the LKJ on theta_5
  and theta_6 is gone; in `frm_sample()` one `theta theta_4` prior
  retires the default `lkj(1)` on all three correlations (a changed
  sampled density), and `sd` with `coef = "x"` plus `theta theta_2`
  drops the sd spec with no message. A filter in `prior_table()` and a
  refusal map in the theta branch of `resolve_priorlist()`, about 40 to
  60 lines in `R/priors.R`, plus tests; the registered covariance
  structures' parameters, which only class theta reaches, keep the
  spelling. It removes the box on a single RE correlation, which brms
  does not have, and flips the ported divergences `priors:5` and
  `priors:123`: a USER DECISION. Record: `dev/acprior-findings.md`,
  "Not done, and why"; `dev/reviews/2026-10-27-acprior.md`, "The same
  double slot elsewhere".
- The sampling table prints `uniform(-1, 1)` for `ar` and `ma` where
  brms prints `(flat)` with `lb`/`ub` (as `sar()`'s row did at 0.75.0).
- frmtmb.sample's `prior = "flat"` still means flat on the internal
  parameter, so for `ar` it is the old `(1 - rho^2)^(-3/2)`: the
  documented meaning of the opt-out.

### log_lik() of residual covariance terms

- **`frm()` gives no message for an interior time gap** under `ar()`,
  `ma()` or `arma()`: `autocor_warn_gaps()` warns only on
  non-consecutive GLOBAL levels, so a ported brms script fits the
  by-level model where brms fits the by-position one. Pre-existing.
  Record: `dev/loglik-findings.md`, Fix round 1, "Not done".
- **brms accepts `se()` and `mi()` beside `ar(cov = TRUE)`**, where
  frmtmb refuses them. Pre-existing. Record: the same.
- **The Jacobian check is literal.** `frmtmb_register_autocor(jacobian
  =)` is checked against the map entry by entry, so a slot that returns
  `O K` for an orthogonal `O`, with the same `K'K`, is refused although
  the precision it gives is right. Record: `dev/reviews/2026-10-27-loglik.md`,
  Fix round 2 check, "Still open".
- **Optional: a documented option for the dense-fallback threshold**
  (`getOption("frmtmb.autocor_dense_max", 1000L)`, named in "The
  Jacobian"), which would also let frmtmb.spatial's `test-sar.R` set it
  instead of reading `frmtmb:::autocor_loo_dense_max`. Record: the same
  review, Fix round 2 check, item 2.
- `fcor()` allocates two `n x n` matrices per draw. Record:
  `dev/loglik-findings.md`, "Open items and defects".

### Found by the consolidation

- **marginaleffects' `vcov =` matrix is not read for a frmtmb fit.**
  `get_vcov.frmtmb_fit(model, ...)` returns `interop_vcov(model)` and
  ignores the argument, so a user's robust or clustered covariance is
  dropped without a message. Measured (`dev/rel076-me-vcov.R`, seed 1,
  `y ~ x` on 80 rows): `predictions(fit, vcov = 2 * get_vcov(fit))`
  gives the standard errors of `predictions(fit)`, ratio 1, 1, 1, where
  the `lm()` control gives 1.414214; the same on rellib-r13 and
  rellib-r14 (`dev/rel076-log/me-vcov-rellib-r1[34].txt`), so it is
  pre-existing. Either read a matrix given, or refuse it.
- RESOLVED, not a defect (the 0.76.0 release review, finding 2):
  **draws under the sampling table passed back as `prior =` are not
  the defaults' draws.** The two prior stacks sum their terms in
  another order: evaluated right after each run, the log densities
  differ at 9 of 500 points by one ulp (1.4e-14), the gradients at none,
  and the warmup paths part at iteration 14 (`dev/review-rel076-trip.R`).
  The consolidation's first comparison, bitwise equal at 800 points,
  was invalid: tmbstan holds the objective in a process-global pointer
  (`dev/rtmb-pitfalls.md` item 23). Kept here as a record; nothing to
  do.
- **The ported block "update has reasonable outputs" now stops at its
  refits.** Its newdata hold `Exp = 4` on every row, so since 0.76.0
  `mo(Exp)` has one category and both refits (setup lines 921 and 935)
  are refused. `brmsfit-methods:924` is a divergence by policy for that
  reason; `:923`, `:927` and `:937` keep their outcome and class
  (divergence) but no longer reach what their reasons describe
  (frmtmb's class name; the `data =` spelling). A later pass could port
  the block with a non-constant `Exp` as an adapted input, so the three
  reasons are observed again. Record: `dev/round-20261027.md`, "The
  ported brms suite".
- `dev/brmsport-bin2-survey.tsv` keeps its hand-written gap note
  ("frmtmb has no sar term; frm() refuses it") on the 13 `sar()` and
  `fcor()` rows the 0.76.0 round ported; the ledger's reason overrides
  it, but the survey reads wrong.
