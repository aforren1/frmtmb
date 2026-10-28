# Lane sigmainn: brms's `sigma` is the innovation sd under `ar/ma/arma(cov = TRUE)`

Worktree `frmtmb-wt-sigmainn` (branch `wt-sigmainn` at 435e7ce6,
frmtmb 0.76.0), library `wt-sigmainn-lib`, base `rellib-r14`. Nothing
is committed. Every number below is from a script in `dev/sigmainn-*`
and its log in `dev/sigmainn-out/` (gitignored; the numbers are pasted
here).

## 1. What brms 2.23.0 does, run (not recalled)

`dev/sigmainn-brms-code.R` (`brms-code.txt`): `stancode()`,
`standata()`, `default_prior()` for every residual term frmtmb has.

| term | brms's `sigma` | frmtmb before | frmtmb after |
|---|---|---|---|
| `ar(cov = TRUE)`, p = 1 | innovation sd: `cholesky_cor_ar1()` is `chol(T / (1 - ar^2))`, scaled by `sigma` in `normal_time_hom_lpdf()` | marginal sd | innovation sd |
| `ma(cov = TRUE)`, q = 1 | innovation sd: `cholesky_cor_ma1()`, diagonal `1 + ma^2` | marginal | innovation |
| `arma(cov = TRUE)` | innovation sd: `cholesky_cor_arma1()`, diagonal `(1 + ma^2 + 2 ar ma) / (1 - ar^2)` | marginal | innovation |
| `ar/ma(cov = TRUE)`, order 2 | brms refuses: "Covariance formulation of ARMA structures is only possible for effects of maximal order one" | marginal | innovation (one meaning at every order) |
| `ar/arma(cov = FALSE)`, order 1 and 2 | sd of a row around the one-step mean `mu*` (`normal_lpdf(Y | mu, sigma)` after the `Err` recursion), the innovation sd | same likelihood, same `sigma` | unchanged |
| `cosy()` | marginal: `cholesky_cor_cosy()` is a correlation | marginal | unchanged |
| `unstr()` | marginal: `Lcortime` is `cholesky_factor_corr` | marginal | unchanged |
| `fcor(M)` | scale of `M`: `multi_normal_cholesky(y | mu, sigma * chol(M))` | `sigma^2 M` | unchanged |
| `sar(lag)`, `sar(error)` | innovation sd: `normal_lagsar_lpdf()` and `normal_errorsar_lpdf()` divide `(I - rho W) y - mu` by `sigma` | innovation (frmtmb.spatial's whitening map divides by `sigma` the same way) | unchanged |
| `car()` | observation residual sd: `normal_lpdf(Y | mu + rcar, sigma)`; the field has `sdcar` | residual sd; field sd is the covstruct's | unchanged |
| `ar(cov = TRUE)` + `sigma ~ x` | `sigma = exp(Intercept_sigma + Xc_sigma * b_sigma)` into `normal_time_het_lpdf()`, row by row: the log innovation sd | log marginal sd | log innovation sd |
| `ar(cov = TRUE)` + `student()` | `student_t_time_hom_lpdf(Y | nu, mu, sigma, Lcortime)`: innovation | marginal | innovation |
| `sar()` + `sigma ~ x` | refused ("SAR models are not implemented when predicting 'sigma'") | refused | unchanged |

So the only terms where frmtmb differed are the covariance form of
`ar()`, `ma()` and `arma()`, with constant `sigma` and under `sigma ~
x`, gaussian and student. Every other term already agreed.

## 2. The route, decided by measuring

Two routes:

- A: the model's own `sigma` parameter (and the `sigma` linear
  predictor) is the innovation sd; the likelihood, the simulator, the
  pointwise conditional and the pearson residual read the marginal sd
  `sigma * exp(autocor_log_marg(theta))`.
- B: keep the marginal sd internally and convert at reporting time.

`dev/sigmainn-route.R` (seed 20261027; arms `before` on rellib-r14,
`after` on the lane build, `ref` an independent dense multivariate
normal written from brms's covariance, `optim()` then Newton steps on a
`numDeriv` Hessian; `route-summary.txt`). Data: ar, the AR(1) 0.9 data
of `dev/acprior-agree.R`; ma, arma as in `dev/acprior-innov.R`; ar2 an
AR(2) (0.5, 0.3) in 15 groups of 8; sigx the ar model with
`sigma ~ x2`, 20 groups of 6.

| model | log sigma_innov: A (after) | B (before, marginal + log c) | ref | SE: A | B (delta) | ref | cor(log sigma, theta_1): marginal param | innovation param |
|---|---|---|---|---|---|---|---|---|
| ar | -0.902022061 | -0.902021793 | -0.902021912 | 0.108668 | 0.108668 | 0.108668 | 0.877 | -0.343 |
| ma | -0.172216023 | -0.172216856 | -0.172216852 | 0.095300 | 0.095300 | 0.095300 | 0.454 | -0.287 |
| arma | -0.096200939 | -0.096200514 | -0.096200762 | 0.093972 | 0.093972 | 0.093972 | 0.690 | -0.083 |
| ar2 | -0.066383898 | -0.066383058 | -0.066383282 | 0.066255 | 0.066255 | 0.066255 | 0.795 | -0.186 |
| sigx (intercept) | -0.864216597 | -0.864214384 | -0.864216383 | 0.065829 | 0.065829 | 0.065829 | 0.695 | -0.189 |

- ML: A, B and the reference agree to optimizer precision (at most
  1.3e-5 relative on the log scale, the size of nlminb's own stopping
  differences; the betas likewise, at most 4.1e-6). The SEs agree to 6
  digits: B's delta method in (log sigma_m, theta) and A's Hessian in
  (log sigma_i, theta) are the same quantity at the optimum.
- logLik: the same model at the same point is an identity.
  `dev/sigmainn-same-point.R` (`same-point.txt`) evaluates the lane
  build's objective at rellib-r14's optimum carried to the innovation
  sd: relative differences 2.18e-16, 0, 1.33e-16, 1.68e-16 and 0 for
  ar, ma, arma, ar2, sigx. The two optima themselves differ by the
  optimizer's path: logLik relative differences 5.3e-14, 3.7e-12,
  6.6e-14, 1.3e-12, 1.6e-11.
- Geometry: the correlation of log sigma with the first correlation
  parameter falls from 0.45 to 0.88 (marginal) to -0.08 to -0.34
  (innovation). That is the parameterization the optimizer and the
  sampler see: route B would keep the 0.88.
- Readers: under B every reader of `sigma` and of the `sigma`
  coefficients must convert (sigma(), VarCorr(), summary(), confint(),
  vcov(), fixef()/coef() of `sigma ~ x`, frm_linpred/fitted(dpar),
  conditional_effects/emmeans/marginaleffects on dpar sigma,
  hypothesis(), the draws, variables(), prior_summary, the
  default-prior tables), the intercept of `sigma ~ x` needs a shift by
  log c with a delta-method SE across components, and the priors keep
  the `lc` layer and need it for the sigma intercept too. Under A four
  internal consumers convert (objective, simulator, pointwise
  conditional, pearson residual, plus the SE check's grid scale), the
  prior layer is deleted, and every reporter is right with no code.

Decision: route A. The 0.76.0 prior adjustment (`autocor_innov()`,
`innov_logc()`, `autocor_innov_sd()`, `prior_logdens(lc =)`, the
`innov` entry field, `target_coefs()`'s `resp`/`dpar` fields, the
refusal of `lb > 0` / `ub` on sigma, frmtmb.sample's `innov_sd` flag
and its message tag) is removed, not left as a dead layer.

## 3. What changed

Core:

- `R/autocor.R`: `autocor_arma_acv()` (the unnormalized
  autocovariance at unit innovation variance; `autocor_arma_acf()` is
  it over its first entry, the same operations as before);
  `autocor_is_innov()`, `autocor_log_marg()` (closed forms at order
  one, `0.5 log gamma(0)` above), `autocor_marg_sigma()`,
  `autocor_marg_dpars()`, `autocor_check_block()`,
  `autocor_check_fit()`; `autocor_block()` marks an innovation block
  `sigma_is = "innovation"`. Docs: the header comment, the
  `?frmtmb-autocor` intro, a new section "The scale of `sigma`", the
  "Divergence from brms" section (only the order limit is left),
  `?autocor_matrix`.
- `R/objective.R`: the autocor density reads the marginal sd.
- `R/families.R` (`sim_autocor_rows()`), `R/autocor-loglik.R`
  (`autocor_row_loglik()`), `R/predict.R` (`residuals()` pearson and
  deviance, and `?residuals.frmtmb_fit`), `R/se-check.R`
  (`se_sd_scale()`): read the marginal sd.
- `R/priors.R`: the innovation prior layer removed (see above);
  `?set_prior` Residual correlation paragraph rewritten.
- `R/sampling-api.R`, `NAMESPACE`: `autocor_innov_sd` no longer
  exported or documented.
- Saved fits: `autocor_check_fit()` beside each `mo_check_fit()` call
  (`require_fit_packages()`, `summary()`, `vcov()`, `coef()`,
  `confint()`) and in `sigma()` and `VarCorr()`; not in `logLik()`,
  whose number is the same model's. `autocor_marg_sigma()` refuses at
  the point of use, so the objective, simulation, `log_lik()` and
  residuals of an old fit or of old draws refuse too.
- `R/compat.R`: the two `group:autocor` residual notes.
- `vignettes/brms-migration.Rmd`: the `sigma` paragraph.
- `NEWS.md`: a development section, Breaking changes.

frmtmb.sample: `R/sample.R` (the `innov_sd` flag and its message line
removed, `?frm_sample` paragraph rewritten), `NEWS.md`.

Tests (every one that pins the change was seen to fail on rellib-r14;
section 5): core `test-autocor.R`, `test-autocor-sigma.R` (new),
`test-autocor-loglik.R`, `test-priors-autocor-classes.R`,
`test-brms-likelihood.R`, `helper-brms.R` (the translator's sigma
scale rule removed); frmtmb.spatial `helper-brms.R` (the same rule);
frmtmb.sample `test-default-priors-autocor.R`,
`test-loglik-autocor.R`, `test-loo.R`, `test-default-priors-fit.R`.

Ported suite: `dev/brmsport-bin2-adapters.R`,
`dev/brmsport-bin2-adapters-spec.R`, and the two regenerated bin-2
files `dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R` and
`-posterior_predict.R` (section 7).

## 4. Every reader of sigma, before and after

`dev/sigmainn-readers.R` (`readers-before.txt`, `readers-after.txt`;
the ar data, seed 20261027):

| reader | before (rellib-r14) | after |
|---|---|---|
| `exp(sigma_(Intercept))`, the parameter | 1.07848 (marginal) | 0.405748 (innovation) |
| `sigma(fit)` | 1.07848 | 0.405748 |
| `summary(fit)` sigma: est, SE, CI | 1.07848 0.229332 [0.710898, 1.63612] | 0.405748 0.0440918 [0.327913, 0.502059] |
| `VarCorr(fit)$residual__$sd` | 1.07848 0.229332 [0.628993, 1.52796] | 0.405748 0.0440918 [0.31933, 0.492167] |
| `confint(fit)` sigma row, exp | [0.710898, 1.63612] | [0.327913, 0.502059] |
| `vcov(full = TRUE)` SE of log sigma | 0.212645 | 0.108668 |
| `frm_linpred(dpar = "sigma")`, `fitted(dpar = "sigma")` | 1.07848 | 0.405748 |
| `insight::get_sigma()` | 1.07848 | 0.405748 |
| `variables(fit)` | `sigma` | `sigma` (same name, now brms's quantity) |
| pearson residuals, sd | 1.00995 | 1.00996 (marginal sd, as before) |
| `simulate()` residual sd (400 draws) | 1.09073 | 1.09073 (marginal, as before) |
| logLik | -32.6326 | -32.6326 |
| `sigma ~ x2`: fixef sigma rows | -0.869391, 0.214814 | -0.871871, 0.214813 (intercept shifts by log c, slope unchanged) |
| `sigma ~ x2`: `conditional_effects(dpar = "sigma")` range | 0.241189 to 0.794320 | 0.240592 to 0.792350 |
| `predict(dpar =)` | refused (brms has none) | refused |
| `hypothesis(fit, "sigma > 0")` | refused, class b (brms's default) | same |
| `frmtmb_autocor_pars(fit)` | 0 rows (lists registered terms only) | same |
| `frmtmb_covstruct_info()` | not a reader of sigma | same |
| `conditional_effects(method = "posterior_predict")` band, half-width / 1.96 (reviewer's AR(1) 0.84 fit, `dev/review-sigmainn-readers.R`) | 0.854937 (marginal sd) | 0.469741 (innovation sd, brms's band; added in fix round 1) |
| frmtmb.sample `conditional_effects(method = "posterior_predict")` on draws | at the draws' sigma (marginal) | at the draws' sigma (innovation, brms's) |
| `cooks.distance()` sum (same fit) | 0.641959 | 0.935047 (the coefficient vector holds log innovation sd; added in fix round 1) |
| `dfbeta()`, `dfbetas()`, `influence()` | read log marginal sd | read log innovation sd |

On the `sigma ~ x2` data of the readers probe (iid errors, so `ar` is
near 0) the two `conditional_effects()` curves differ by one factor,
0.99753 at both ends (0.240592 / 0.241189, 0.792350 / 0.794320): the
factor `c` of the fitted `ar`. The slope moves by 1e-6 (optimizer).

Readers with no change of code that now report the innovation sd by
construction: `summary()`, `print()`, `VarCorr()`, `sigma()`,
`fixef()`/`coef()` and the parameter table, `confint()`,
`confint_varcorr()`, `vcov()` and ML SEs, `fitted()`/`predict()`/
`frm_linpred()` with `dpar = "sigma"`, `conditional_effects()`,
emmeans and marginaleffects on dpar sigma, insight's `get_sigma()`,
`hypothesis()` names, `prior_summary()`, `default_prior()`/
`validate_prior()` tables (the sigma rows are unchanged text; they now
mean brms's parameter), frmtmb.sample's draws, `variables()`,
`posterior_summary()`, `as_draws*()`, `summary()` of draws.

Readers that read the marginal sd (changed code): the objective,
`simulate()`/`posterior_predict()` (core's `sim_autocor_rows()`, which
frmtmb.sample's draws reach through `draws_fit_at()`), `log_lik()`
(`autocor_row_loglik()`), pearson and deviance residuals, the SE
check's grid scale.

## 5. Tests that pin the change, seen to fail on rellib-r14

One file per process, `FRMTMB_BRMS_FIT_TESTS=true`, the lane's test
files against rellib-r14 (tier `sigmainn-before1`; every `lib:` line
names rellib-r14) and against the lane library (`sigmainn-f2`,
`sigmainn-f3`):

    RESULT test-autocor.R pass=184 fail=8 err=0 skip=0 warn=0            (r14)
    RESULT test-autocor-sigma.R pass=26 fail=23 err=0 skip=0 warn=1      (r14)
    RESULT test-autocor-loglik.R pass=53 fail=13 err=0 skip=0 warn=0     (r14)
    RESULT test-priors-autocor-classes.R pass=163 fail=6 err=1 skip=0 warn=0  (r14)
    RESULT test-brms-likelihood.R pass=474 fail=8 err=0 skip=0 warn=0    (r14)
    RESULT test-default-priors-autocor.R pass=91 fail=7 err=0 skip=0 warn=0   (r14, frmtmb.sample)
    RESULT test-loglik-autocor.R pass=103 fail=36 err=0 skip=0 warn=0    (r14, frmtmb.sample)
    RESULT test-loo.R pass=143 fail=3 err=0 skip=0 warn=0                (r14, frmtmb.sample)
    RESULT test-default-priors-fit.R pass=26 fail=2 err=0 skip=0 warn=0  (r14, frmtmb.sample)

The failures are behavioral: brms's log density at frmtmb's estimates
misses frmtmb's logLik by 9.34, 8.69, 14.28 and 13.64 nats (ar, ma,
arma, ar2) with gradients 46 to 54 where the lane's is below 1e-3 of
that; brms's `log_prob()` at frmtmb's translated parameters misses by
5.28 nats (gradient 58) for ar, 1.91 (37) ma, 5.35 (58) arma, 5.28
(58) ar + `sigma ~ x`; the sigma prior density moves with theta by
0.209, 0.089, 0.153; a bound on sigma is refused; the frm_sample sigma
mean is 65.4 Monte Carlo standard errors from brms's; the dense
conditional of each group differs at every draw; the saved-fit refusal
is absent. On rellib-r14 test-autocor-sigma.R also lets one warning
escape (a singular-convergence message of one fit), which the lane
build does not raise.

## 6. frm_sample() against brms's own draws

`dev/sigmainn-agree.R` (`agree-summary.txt`; 4 chains of 4000, seed
3, each side; brms's raw stanfit from `stancode()`/`standata()`
through the Stan cache; nothing weighted or converted). User priors:
`normal(0, 0.5)` on sigma with `ub = 1` and `normal(0.5, 0.3)` on ar
(ar_user); `normal(-1, 0.5)` on the sigma intercept and `normal(0,
0.5)` on the sigma slope under `sigma ~ x2` (sigx).

    #### after (lane build)
    ar     ar[1] z 0.14 (sd ratio 1.010); sigma 0.95 (1.007); b_Intercept -0.82; b_x 0.37
    ma     ma[1] 1.42 (1.037); sigma -0.44 (0.977); b_Intercept -0.54; b_x -0.80
    arma   ar[1] 1.36; ma[1] -1.53; sigma 0.51 (1.006); b_Intercept -0.69; b_x 1.41
    ar_user ar[1] 0.47; sigma 0.32 (0.979); b_Intercept -1.80; b_x 0.64
    sigx   ar[1] -0.55; b_sigma_Intercept 1.54 (1.001); b_sigma_x2 0.15; b_Intercept -0.58; b_x 0.88

    #### before (rellib-r14, draws as reported)
    ar     sigma z 122.36 (sd ratio 5.421): mean 1.1540 against brms 0.4324
    ma     sigma 76.54; arma sigma 185.39
    ar_user refused: "class = "sigma" takes no ub here ..."
    sigx   ar[1] -7.95; b_sigma_Intercept 199.58 (-0.5843 against -0.8481)

Every |z| after is at most 1.80 and every sd ratio within 0.969 to
1.037; 0 divergent transitions on either side in every run. Bulk ESS
of sigma before/after: ar 2967/5982, ma 3757/3513, arma 3052/3797.
The sigx "before" row is the 0.76.0 sigma-intercept prior defect: the
prior read at the log marginal sd moved `ar` by 7.95 Monte Carlo
standard errors.

## 7. log_lik() and the ported brms rows

- `log_lik()` against a dense reference: core `test-autocor-loglik.R`
  (ar, ar(2), ma, arma, cosy, ar + RE, ar + `sigma ~ x`; gaussian and
  student) with the reference covariance built from
  `stats::ARMAacf()` times the psi-weight process variance; frmtmb.sample
  `test-loglik-autocor.R` (the dense conditional at hand-set draws, and
  brms's own `log_lik()` at the same parameter vectors with brms's
  sigma now passed unchanged; under `sigma ~ x` brms's
  `Intercept_sigma` with no shift).
- The four rows the coordinator gave this lane, ported through the
  bin-2 route (`dev/brmsport-bin2-adapters.R`: `brms_adapt_time_model()`,
  `brms_adapt_acframe()`, `brms_adapt_ac_any()`; spec:
  `gaussian_time`/`student_time` in `brms_adapt_ac_spec_funs`, the two
  blocks in `bin2_adapted`, `brms_adapt_time_subst`, the four ids in
  `bin2_adapted_input` with their own reason).
  `dev/sigmainn-bin2-brmsprep.R` (`bin2-brmsprep.txt`) shows why an
  adaptation is needed: brms's preps carry no acframe, so
  `has_ac_class(NULL, "arma")` is FALSE and brms's own
  `get_cov_matrix_ac()` builds `sigma^2 I + diag(se^2)` (draw 1:
  off-diagonal 0, diagonal 134.5 = 4.26 + 130.24); `ar` and `ma` never
  enter. The port drops `se` (frmtmb refuses `se()` beside a cov =
  TRUE term; the adapter refuses a prep that still carries it) and
  adds an acframe naming arma, so brms's value check reads the ARMA
  covariance. An `|ma| > 1` (posterior_predict's `rnorm(ns, 0.2, 1)`)
  is mapped to `1 / ma` with sigma times `|ma|`, the same covariance.
- `Rscript dev/brmsport-bin2-gen.R`: "wrote 6 files; 242 assertion
  ids; 4 without a verdict; 157 adapted; 73 value checks". Recorded on
  the lane build (`dev/sigmainn-bin2-record.sh`, records in
  `dev/sigmainn-bin2-log/`):

        RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=101 fail=0 err=0 skip=0 blocks=16
        RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=63 fail=0 err=0 skip=0 blocks=15
        assert frmtmb log_lik:130 unclassified TRUE
        value  frmtmb log_lik:130 HOLDS: 80 values (0 NA, 0 infinite), max abs diff 1.35e-13, max rel 2.39e-15
        assert frmtmb log_lik:133 unclassified TRUE
        value  frmtmb log_lik:133 HOLDS: 80 values (0 NA, 0 infinite), max abs diff 1.42e-14, max rel 3.4e-15
        assert frmtmb posterior_predict:93 unclassified TRUE
        assert frmtmb posterior_predict:97 unclassified TRUE

  On rellib-r14 the same files: the four assertions hold (shapes), and
  the value checks differ (log_lik:130 max rel 2.23, :133 max rel
  1.31): 0.76.0 read brms's sigma as the marginal sd.
- Dry-run ledger (`FRMTMB_BIN2_REC_DIR=dev/sigmainn-bin2-log
  FRMTMB_BIN2_LEDGER_OUT=dev/sigmainn-ledger-out`, `ledger-dry.txt`):
  against the tracked bin-2 ledger exactly four rows move:

        log_lik:130            cannot transfer / absent      -> pass / adapted-input
        log_lik:133            cannot transfer / absent      -> pass / adapted-input
        posterior_predict:93   cannot transfer / not-ported  -> pass / adapted-input
        posterior_predict:97   cannot transfer / not-ported  -> pass / adapted-input

  Bins 1 and 2: 619 -> 623 of 823 (75.7%); without the decisions 623 of
  706 (88.2%); gaps 4 -> 2, not ported (cost) 8 -> 6.
- Tier run with the dry-run verdicts and value-check verdicts swapped
  in temporarily (then restored; `git diff` of
  `dev/brmsport-bin2-verdicts.tsv` and `-valuechecks.tsv` is empty, and
  the generated files were regenerated from the tracked verdicts, so
  they carry the four ids as "unclassified" for the consolidation):

        RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=101 fail=0 err=0 skip=0 blocks=16
        RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=63 fail=0 err=0 skip=0 blocks=15

- Bin 1 (`dev/brmsport-verdicts.tsv`): no row reads `sigma` under a
  cov = TRUE term (grep of the generated `test-brms-suite-*.R`: the one
  `cov = TRUE` is an error-message check in `brm:124`), so no bin-1 row
  flips. Not edited, not regenerated.

## 8. Models without these terms

`dev/sigmainn-unchanged.R` (`unchanged-before.txt` against
`unchanged-after.txt`): gaussian, random intercept, `sigma ~ z`,
student + RE, poisson + RE, cosy, cosy + `sigma ~ z`, unstr, `ar()` and
`arma()` cov = FALSE, fcor, multivariate rescor; logLik, every internal
estimate, every SE (vcov full) and sigma to 12 digits; MAP fits with
sigma, b and dpar-sigma intercept priors (objective and estimates); and
`frm_sample()` draws (1 chain of 300, seed 5) of `y ~ x` and
`y ~ x + cosy(t, g)`: 69 lines, identical past the library line.

<!-- the suite runs, sessions, OpenBLAS and checks are appended below -->

## 9. Every suite, one file per process

Plain tier (`NOT_CRAN=true`; `dev/release/run-par.sh sigmainn-plain`,
8 processes; summed by `dev/sigmainn-sum.R`, `plain-sum.txt`):

    frmtmb           files=282 pass=20128 fail=0 err=0 skip=205 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1838 fail=0 err=0 skip=5 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=430 fail=0 err=0 skip=13 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=70 pass=3280 fail=0 err=0 skip=14 warn=0
    frmtmb.spatial   files=8 pass=478 fail=0 err=0 skip=6 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs without a RESULT line: 0
    logs whose frmtmb is not the lane lib: 0

Gated tier, every gate on (`dev/ciharden-scan.sh sigmainn-gated ref
<lane lib> 8`; `gated-sum.txt`):

    frmtmb           files=282 pass=22705 fail=0 err=0 skip=0 warn=0
    frmtmb.coupling  files=11 pass=557 fail=0 err=0 skip=0 warn=0
    frmtmb.eam       files=32 pass=1862 fail=0 err=0 skip=0 warn=0
    frmtmb.latent    files=12 pass=687 fail=0 err=0 skip=0 warn=0
    frmtmb.learn     files=15 pass=515 fail=0 err=0 skip=0 warn=0
    frmtmb.ode       files=11 pass=551 fail=0 err=0 skip=0 warn=0
    frmtmb.sample    files=70 pass=3619 fail=0 err=0 skip=0 warn=0
    frmtmb.spatial   files=8 pass=561 fail=0 err=0 skip=0 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs without a RESULT line: 0
    logs whose frmtmb is not the lane lib: 0

Against the 0.76.0 release's gated counts (`dev/round-20261027.md`,
r076-ref2: frmtmb 22637, sample 3610, spatial 561, the other six
equal): core +68, sample +9, every other suite equal. Every changed
file accounted for, HEAD's copy against the lane's, both gated
(`dev/sigmainn-head-tests/`, tier `sigmainn-head` on rellib-r14, and
the gated tier above):

| file | 0.76.0 | lane | why |
|---|---|---|---|
| core `test-autocor.R` | 192 | 192 | five comparisons with gls()/lme() and simulate() now carry the marginal factor |
| core `test-autocor-sigma.R` | (new) | 49 | dense brms reference at the estimates, the readers, the saved-fit refusal |
| core `test-autocor-loglik.R` | 66 | 66 | the dense reference reads the marginal sds |
| core `test-priors-autocor-classes.R` | 166 | 179 | the two sigma-prior tests rewritten: six formulas, the bound taken (+13) |
| core `test-brms-likelihood.R` | 476 | 482 | row 18 hands brms sigma unchanged; row 18e (ma, arma, ar + sigma ~ x) +6 |
| sample `test-default-priors-autocor.R` | 89 | 98 | the default on sigma itself, five formulas, the sigma ~ x intercept (+9) |
| sample `test-loglik-autocor.R` | 139 | 139 | dense and brms references read brms's sigma unchanged |
| sample `test-loo.R` | 146 | 146 | the dense reference's covariance over 1 - rho^2 |
| sample `test-default-priors-fit.R` | 28 | 28 | sigma's draws compared with brms's as they are |

No test was dropped. The helper changes (`helper-brms.R` in core and
frmtmb.spatial) remove the translator's sigma scale rule; frmtmb.spatial's
gated count is unchanged (561).

Plain tier against 0.76.0's (`dev/round-20261027.md`, r076-plain:
frmtmb 20066 pass, 204 skip): +62 passes (test-autocor-sigma.R 49,
test-priors-autocor-classes.R 13) and +1 skip (the new gated block
"row 18e" of test-brms-likelihood.R, skipped without
`FRMTMB_BRMS_FIT_TESTS`). frmtmb.sample 3280 against 3271 (0.76.0's
3263 plus the release review's 8): +9, test-default-priors-autocor.R.
frmtmb.spatial 478: 0.76.0's 466 plus the release review's 12.

## 10. One session per suite, three orders

Plain tier, all nine suites (`dev/release/run-sessions.sh sigmainn-s1 8`,
`dev/release/sigmainn-s1.log`); every session loads frmtmb and
frmtmb.sample from the lane library and the other seven from
rellib-r14 (27 `lib:` lines each):

    SESSION frmtmb forward files=272/282 pass=20128 fail=0 err=0 skip=205 warn=0
    SESSION frmtmb reverse files=272/282 pass=20278 fail=0 err=0 skip=205 warn=0
    SESSION frmtmb preload files=272/282 pass=20371 fail=0 err=0 skip=207 warn=0
    SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample forward files=65/70 pass=3280 fail=0 err=0 skip=14 warn=0
    SESSION frmtmb.sample reverse files=65/70 pass=3280 fail=0 err=0 skip=14 warn=0
    SESSION frmtmb.sample preload files=65/70 pass=3280 fail=0 err=0 skip=14 warn=0
    SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
    SIGMAINN-S1 ran 27 of 27
    LEAK lines: 46, untagged: 0

The core reverse and preload sessions pass more (+150, +243) as at
0.76.0, where some core tests assert more with extensions loaded; the
preload skips (+2) are the arms that need an extension absent.

## 11. OpenBLAS with LAPACK

`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack` (La_version
3.12.0), then `dev/ciharden-scan.sh sigmainn-ob26 ob0.3.26-lapack-t4`
and `sigmainn-ob32 ob0.3.32-lapack-t4` on the 12 changed or numeric
files (`dev/sigmainn-ob.jobs`), every gate on, 4 threads; each log's
`blas:` line names the OpenBLAS copy (12 of 12 each):

    test-autocor.R 192, test-autocor-sigma.R 49, test-autocor-loglik.R 66,
    test-priors-autocor-classes.R 179, test-autocor-seam.R 80,
    test-autocor-cond.R 57, test-se-check.R 134, test-brms-likelihood.R 482,
    sample test-default-priors-autocor.R 98, test-loglik-autocor.R 139,
    test-loo.R 146, test-default-priors-fit.R 28

fail 0, err 0, skip 0, warn 0 under both versions, the same RESULT
lines as the reference BLAS gated tier. No new assertion sits on a
rounding-decided branch: the tolerances are ratios of order 1e-10 to
1e-3 against differences that are 0 or of order 1 before the change.
The OpenBLAS copies of R (`dev/ciharden-out/`) are deleted.

Gated, every gate on, frmtmb, frmtmb.sample and frmtmb.spatial (the
packages this lane changes; `sigmainn-s2g`, 6 at a time):

    SESSION frmtmb forward files=282/282 pass=22705 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb reverse files=282/282 pass=22855 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb preload files=282/282 pass=22948 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.sample forward files=70/70 pass=3619 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.sample reverse files=70/70 pass=3619 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.sample preload files=70/70 pass=3619 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=561 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=561 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=561 fail=0 err=0 skip=0 warn=0
    SIGMAINN-S2G ran 9 of 9
    LEAK lines: 157, untagged: 0

## 12. Floors, checks, and what is left

Floors (`dev/release/run-par.sh`, copies of installed builds, gated):
the lane's frmtmb.sample over frmtmb 0.76.0 (rellib-r14): its four
changed files fail 47 assertions (default-priors-autocor 6,
loglik-autocor 36, loo 3, default-priors-fit 2). frmtmb.sample 0.23.0
over the lane's core: 28 errors `could not find function
"autocor_innov_sd"`. So frmtmb.sample needs the new frmtmb, and the old
frmtmb.sample does not work on it.

R CMD check --as-cran (`dev/sigmainn-check.sh`):
- frmtmb.spatial: Status: 1 NOTE (V8 math rendering), tests FAIL 0 WARN
  0 SKIP 7 PASS 410.
- frmtmb: NOT COMPLETED. The build failed ("copying to build directory
  failed ... Permission denied") because it ran beside the
  frmtmb.spatial check, whose open files under dev/sigmainn-check/ the
  core build tried to copy (R CMD build copies the tree before
  .Rbuildignore applies). Rerun alone: `bash dev/sigmainn-check.sh .`.
- frmtmb.sample: started; its status was not read before the report.
  Log: dev/sigmainn-check/frmtmb.sample/check.log.

Readers of sigma that stay on the marginal scale, by design: the
objective, simulate()/posterior_predict(), log_lik(), pearson and
deviance residuals, the SE check's grid scale. No reported value or
draw column of sigma is the marginal sd.

### Checks, completed after the first report

Run one at a time (`dev/sigmainn-check.sh`):
- frmtmb: Status: 1 NOTE (V8 math rendering); tests FAIL 0 WARN 0
  SKIP 233 PASS 19881 (0.76.0: SKIP 232, PASS 19819; +62 passes from
  test-autocor-sigma.R and test-priors-autocor-classes.R, +1 skip from
  the gated row 18e).
- frmtmb.sample: Status: 1 NOTE (V8); FAIL 0 WARN 0 SKIP 22 PASS 3266
  (0.76.0: 3257; +9, test-default-priors-autocor.R).
- frmtmb.spatial: Status: 1 NOTE (V8); FAIL 0 WARN 0 SKIP 7 PASS 410.
The .Rcheck trees and tarballs are deleted; the build and check logs
stay in dev/sigmainn-check/. dev/sigmainn-head-tests/ and
dev/sigmainn-ledger-out/ are deleted.

## Fix round 1 (review `dev/reviews/2026-10-27-sigmainn.md`: FIX, no blocker)

1. **Major: `conditional_effects(method = "posterior_predict")`.**
   Decision: keep brms's behavior. brms 2.23.0 drops every
   autocorrelation term in `conditional_effects()` (`incl_autocor =
   FALSE`) and draws iid at its sigma, the innovation sd; frmtmb's grid
   path draws at the dpar sigma, which is now that sd. Changes: the
   comment above the `autocor_cond_strip()` call in
   `conditional_effects()` and that function's own doc say so;
   `?frmtmb-autocor` ("The scale of `sigma`") says so; NEWS bullets in
   both packages (Breaking changes); reader-table rows (section 4).
   Tests, seen to fail on rellib-r14:
   - core `test-autocor-sigma.R`, "conditional_effects(method =
     'posterior_predict') draws at brms's sigma": the band's half-width
     over 1.96 against the innovation sd read from the residuals and
     the fitted ar (the reviewer's data, seed 20261028, 4000 draws).
     rellib-r14: ratio 1.57 (bound: within 0.06 of 1).
   - frmtmb.sample `test-loglik-autocor.R`, "conditional_effects(method
     = 'posterior_predict') on draws draws at brms's sigma" (gated by
     `skip_sampler()`): against the innovation sd of the simulated data,
     0.5. rellib-r14: ratio 1.49 (bound: within 0.17 of 1).
2. **Minor: influence measures.** NEWS bullet (`cooks.distance()`,
   `dfbeta()`, `dfbetas()`, `influence()`; Cook's sum 0.642 to 0.935 on
   the reviewer's fit); reader-table rows.
3. **Minor: saved 0.76.0 fits.** The refusal now follows one rule:
   refuse where the answer reads `sigma`'s value or its coefficients,
   answer where it reads only the mean. `autocor_check_fit()` left
   `require_fit_packages()` (which mean-only methods also call) and is
   called by: `summary()`/`print()`, `sigma()`, `VarCorr()`,
   `confint()`, `vcov()` (always with `full = TRUE`, else only under
   `sigma ~ x`), `fixef()` and `coef()` (only under `sigma ~ x`),
   `predict()` (its draws read sigma; it used to warn about the
   covariance first), `influence()` and so `cooks.distance()`,
   `dfbeta()`, `dfbetas()`, `hypothesis()` naming sigma,
   `insight::get_parameters()` and `get_varcov()`; through
   `autocor_check_dpar()` by `fitted()`, `frm_linpred()`, `emmeans()`
   and `conditional_effects()` with `dpar = "sigma"`; and by
   `conditional_effects(method = "posterior_predict")`. Pearson and
   deviance residuals, `simulate()`, `log_lik()` and the objective
   refuse at the point of use; `residuals(type = "response")` no longer
   reads `autocor_marg_dpars()` and answers. New
   `.get_sigma.frmtmb_fit()` (registered on insight's internal
   `.get_sigma()`, which `insight::get_sigma()` dispatches on) returns
   `sigma()` as insight's default does when `sigma()` answers, and lets
   its refusal through; before, insight's fallback gave 1.23088 on the
   reviewer's saved fit. It depends on an internal generic of insight
   1.6.0; if insight renames it, the registration is lost and
   `get_sigma()` falls back as before. Both NEWS files say what refuses
   and what answers, and that saved 0.76.0 draws still report the old
   sigma in `summary()`. The saved-fit test now checks every refusing
   reader, the mean-only readers answering with identical values,
   `sigma ~ x` (fixef, vcov, coef refuse), and the absent cases
   `cosy()` and `ar(cov = FALSE)` (struct "ar" too).
4. **Minor: bin-2 reason.** `bin2_adapted_input_time_reason`, the spec
   header and the adapter comment now say: brms's rows exercise the se
   covariance (`sigma^2 I + diag(se^2)`, no acframe); the port swaps it
   for the ARMA covariance; brms 2.23.0 fits `se()` beside `cov = TRUE`
   (`normal_time_hom_se_lpdf()`) and frmtmb refuses it, a parity gap.
   The class stays adapted-input. The regenerated bin-2 files are
   byte-identical (the ledger alone reads the reason). Filed in
   `dev/test-backlog.md`, new section "Filed in the round opened
   2026-10-27, after 0.76.0". The compat row (`R/compat.R`, `se()` of
   `ar()`, `ma()`, `arma()`) and `?frmtmb-autocor` ("What cannot be
   combined with it", "Divergence from brms") say brms fits it.
5. **Nit.** `autocor_is_innov()` cites `dev/sigmainn-brms-code.R`.

Runs (one file per process, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`,
lane library; `sigmainn-fr1`, then `sigmainn-fr2` after the `predict()`
check): every touched file and the reviewer's list, plus the files
that exercise insight, influence, emmeans, interop, extension-home,
`mo()` and sugar:

    test-autocor-loglik.R 66, test-autocor-sigma.R 71, test-autocor.R 192,
    test-brms-likelihood.R 482, test-brms-suite-emmeans.R 11,
    test-compat.R 767, test-emmeans.R 110, test-extension-home-mixed.R 359,
    test-extension-home.R 31, test-influence-plot.R 21, test-interop.R 45,
    test-message-uniqueness.R 6, test-mo-coding.R 122,
    test-priors-autocor-classes.R 179, test-sugar.R 45;
    sample test-conditional-effects-draws.R 76,
    test-default-priors-autocor.R 98, test-default-priors-fit.R 28,
    test-evidence-ratio.R 36, test-laplace-draws.R 93,
    test-loglik-autocor.R 140, test-loo.R 146; spatial test-sar.R 126.
    fail 0, err 0, skip 0, warn 0 in every file, except one escaped
    predict() covariance warning in test-autocor-sigma.R in fr1, which
    the predict() check removed (fr2: 71, warn 0).

On rellib-r14 (`sigmainn-fr3b`): test-autocor-sigma.R pass 37 fail 34,
frmtmb.sample test-loglik-autocor.R pass 103 fail 37 (the new band
test among them).

## For consolidation (from the review, "Consolidation: what must change and be tested")

reslot (`frmtmb-wt-reslot`):
- `extensions/frmtmb.sample/R/sample.R` `flat_spec()` sets
  `s$innov_sd <- NULL`; the field no longer exists after this lane.
  Drop the line. This lane's removals of the `innov_sd` blocks in
  `default_priors_for()` and `announce_default_priors()` apply cleanly.
- reslot's `extensions/frmtmb.sample/tests/testthat/test-flat-prior.R:115`
  asserts `expect_false(is.null(sg[[1L]]$innov))`, which FAILS after
  the merge. Rewrite it: (flat) on class `sigma` under `ar(cov = TRUE)`
  is flat on the innovation sd itself. The entry has no `innov` field,
  and the negative log prior at a fixed `betad` is the same at two
  `thetaac` values (the log-Jacobian of `exp` only, as brms's
  `real<lower=0> sigma` with no prior).
- `R/priors.R`: reslot's hunk at base line 3839 (the ordinal Intercept
  bound message) ends at `for (tg in target_coefs(s)) {`, where this
  lane's removal of the `innov` block begins. Expect a textual
  conflict, and resolve it to this lane's loop body. Base line 3264
  (`target_coefs()` loses `resp`/`dpar`) is near reslot's hunk at
  3245; reslot's new code reads neither field (grep).
- Both lanes edit `tests/testthat/test-priors-autocor-classes.R`,
  frmtmb.sample's `test-default-priors-autocor.R`, `man/set_prior.Rd`,
  `frm_sample.Rd` and both NEWS files: merge, re-roxygenise, rerun both
  test files gated on the merged build.

priordraws (`frmtmb-wt-priordraws`):
- No overlapping hunks in `R/objective.R`, `R/sampling-api.R`,
  frmtmb.sample `R/sample.R`, `test-loo.R` or `brms-migration.Rmd`;
  NEWS conflicts at the top.
- priordraws review finding 2 (sigma under `sample_prior = "only"` z
  -7.95; `prior_sigma` against a marginal `sigma` column in
  `Evid.Ratio`) is what this lane removes. On the merged build rerun
  the priordraws reviewer's O3 (`review-priordraws-only.R brms 2` and
  `frm 2`), add the test the priordraws review asks for (`"only"`
  sigma under `ar(cov = TRUE)` with `exponential(1)` matches its
  prior), and check that `hypothesis(x, "sigma = 1", class = "")`
  divides densities of one quantity.
- The frmtmb.sample floor must be the merged core: an old
  frmtmb.sample calls `autocor_innov_sd()` (28 errors measured,
  section 12).

No reader of sigma's reported value stays on the marginal scale. The
marginal sd is read only where the model needs it: the objective, the
simulator (`simulate()`, `posterior_predict()`, `predict()`),
`log_lik()`, pearson and deviance residuals, and the SE check's grid.

### Fix round 1: sessions and checks

`test-arg-refusal.R` ("no registered S3 method has a `...` it never
touches") failed in the first gated core sessions (`sigmainn-s3g`, 1
failure each order) on the new `.get_sigma.frmtmb_fit`; it joins the
file's exemption list with its reason (insight forwards its caller's
dots), as `get_parameters.frmtmb_fit` does. Rerun one file per process:
122 pass, 0 fail.

One session, three orders, every gate on:

    sigmainn-s4g (after the fix):
    SESSION frmtmb forward files=282/282 pass=22727 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb reverse files=282/282 pass=22877 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb preload files=282/282 pass=22970 fail=0 err=0 skip=2 warn=0
    sigmainn-s3g:
    SESSION frmtmb.sample forward files=70/70 pass=3620 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.sample reverse files=70/70 pass=3620 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.sample preload files=70/70 pass=3620 fail=0 err=0 skip=0 warn=0
    LEAK lines untagged: 0 (113 and 148 tagged)

R CMD check --as-cran, one at a time, after fix round 1:
- frmtmb: Status: 1 NOTE (V8); FAIL 0 WARN 0 SKIP 233 PASS 19903.
- frmtmb.sample: Status: 1 NOTE (V8); FAIL 0 WARN 0 SKIP 22 PASS 3267.
- frmtmb.spatial: no file changed in fix round 1; its check before it
  stands (1 NOTE, PASS 410).
The .Rcheck trees and tarballs are deleted.
