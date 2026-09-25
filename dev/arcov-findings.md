# Lane arcov: brms's `cov = FALSE` form of `ar()`, `ma()`, `arma()`

Round of 2026-09-25. Gap closed: `ma()` (and `ar()`, `arma()`) without
`cov = TRUE`. Before this lane every such term was refused with "only
the residual-covariance formulation is implemented, so the call needs
cov = TRUE"; that message was also what the brms-suite row brm:110
matched for `ma(x)` on poisson, where brms says something else.

## 1. What brms fits (read from brms 2.23.0, not recalled)

`brms::make_stancode()` for `y ~ x + arma(time, g)`, `ar(time, g)`,
`ma(time, g)` and `ar(time, g, p = 2)` (gaussian and student), plus the
functions `stan_ac()`, `data_ac()`, `order_data()`, `frame_ac.btl()`,
`.predictor_arma()` and `predictor_ac()` printed from the installed
namespace (`/tmp/lanes/arcov/arcov-stancode.R`, `arcov-brmsfuns.R`,
`arcov-combos.R`). The facts the implementation rests on:

- The model block, per row `n` of the SORTED data:

      mu[n] += Err[n, 1:Kma] * ma;
      err[n] = Y[n] - mu[n];
      for (i in 1:J_lag[n]) Err[n + 1, i] = err[n + 1 - i];
      mu[n] += Err[n, 1:Kar] * ar;

  then the family's ordinary per-row density at the shifted `mu`. The
  MA part enters `err`, the AR part does not. With
  `e_t = r_t - sum ma_i e_{t-i}` this is `(1 - ar B)(y - m) =
  (1 + ma B) eps`, a standard ARMA, conditioned on `e_s = 0` before a
  group's first row.
- Sort order: `order_data()` sorts by `order(gr, time)` (raw values;
  `time = NA` means data order), refuses duplicate `(gr, time)`.
- `J_lag` (`data_ac()`): `J_lag[n] = sum(tgroup[n:max(1, n + 1 -
  max_lag)] == tgroup[n + 1])`, 0 on the last row. So the lag is
  counted in ROWS within the group; a missing time point is not a gap.
  frmtmb's `cov = TRUE` form counts lags in the GLOBAL set of time
  levels (nlme's reading). The two readings agree only when no group
  skips a level. Documented in `?frmtmb-autocor` and the vignette.
- Parameters: `vector[Kar] ar; vector[Kma] ma;` with NO bounds and a
  flat default prior (`get_prior()` shows `(flat) ar`). A user bound
  `lb/ub` becomes `vector<lower=, upper=>`. frmtmb therefore uses the
  identity: `thetaac` IS `c(ar, ma)`.
- Link: with `gaussian(link = "log")` the Stan code applies `mu =
  exp(mu)` BEFORE the ARMA loop, so the term is added on the response
  scale. brms's R-side `.predictor_arma()` instead adds it to `eta`
  (link scale), so brms's own predictions disagree with its likelihood
  for a non-identity link. frmtmb follows the likelihood everywhere.
- Families: `has_natural_residuals()` is TRUE only for gaussian and
  student (the "residuals" special, checked over every `.family_*`).
  For any other family `stan_ac()` refuses an MA part ("Please set cov
  = TRUE when modeling MA structures for this family") and gives an AR
  part LATENT residuals (`err = sderr * zerr`, an N-vector of iid
  effects added to `eta` with ar-weighted lags of itself).
- Combinations brms generates code for under `cov = FALSE`:
  `weights()`, `cens()`, `trunc()`, `mi()` (the residual is against
  `Yl`), `sigma ~ x`, predicted `nu`, random effects, `rescor = TRUE`
  (multi_normal at the shifted mu). Refused: `se()` ("Please set cov =
  TRUE in ARMA structures when including known standard errors"),
  mixtures and any dpar but mu ("Explicit covariance terms can only be
  specified on 'mu'"), nonlinear formulas (the term is not parsed as
  autocorrelation there).
- Post-fit: `predictor_ac()` runs the same recursion over the OBSERVED
  `Y` for `posterior_epred()`, `posterior_linpred()` and
  `posterior_predict()` (one-step conditional mean); on newdata it uses
  newdata's response, and draws recursively where it is NA.
  `conditional_effects()` and `emmeans()` pass `incl_autocor = FALSE`.

## 2. What changed

- `R/autocor.R`: `parse_autocor_call()` keeps `cov` instead of refusing;
  `autocor_is_cond()`; `check_autocor_cond()` (families, `se()`,
  matrix response, mixture) with brms's own sentences quoted;
  `autocor_cond_positions()` / `autocor_cond_block()` (per-position row
  lists, longest group first, and per-lag gather indices);
  `autocor_cond_shift()` (AD-safe); natural-scale names and "raw"
  interval types; `autocor_matrix()` refusal; post-fit helpers
  (`autocor_cond_mu()`, `autocor_cond_dpars()`, `linpred_arma_cond()`,
  `autocor_cond_fd_se()`, `autocor_cond_newdata_y()`,
  `autocor_cond_strip()`); roxygen: a new section "The
  residual-regression form (cov = FALSE)", the `cov` argument, an
  example.
- `R/objective.R`: the shift is applied to `mu` before the rescor and
  rowwise densities; the joint `cov = TRUE` branch is taken only for
  covariance blocks.
- `R/predict.R`: `frm_linpred()` (so `fitted()`, `predict(type=)`)
  routes mu of a `cov = FALSE` response to `linpred_arma_cond()`; the
  truncated-mean route and `residual_values()` use the one-step mean;
  `residuals(type = "osa")` refusal of its own.
- `R/predict-brms.R`: `predict()` draws rowwise around the one-step
  mean (in sample and on newdata that carries the response).
- `R/families.R`: `sim_autocor_cond()`, the sequential simulator,
  reached by `simulate()`, `frm_simulate()`, `dharma_residuals()`.
- `R/simulate-newdata.R`: `simulate(newdata =)` rebuilds the positions
  on newdata; every group starts from an empty past.
- `R/conditional-effects.R`: drops the term (brms's `incl_autocor =
  FALSE`).
- `R/priors.R`: identity map for `class = "ar"/"ma"` on a `cov = FALSE`
  block; `lb`/`ub` allowed at any order there.
- `R/compat.R`: rows naming `ar()`, `ma()`, `arma()` against aterms,
  `se()`, `student`, `rescor`, `fitted`, `predict`, `simulate`,
  `residuals`; group notes updated so none says `cov = FALSE` is
  refused.
- `R/fit.R`, `R/importance.R`, `R/sandwich.R`: the three refusals
  (quadrature, importance, `vcov_cluster()`) still refuse, with wording
  true of both forms.
- `extensions/frmtmb.sample/R/methods-draws.R`: `posterior_predict()`
  draws around the one-step mean for a `cov = FALSE` term. Without
  this, the extension drew fresh sequential replicates (what
  `simulate()` does), which is not brms's `posterior_predict()` and not
  core `predict()`.
- Docs: `man/frmtmb-autocor.Rd` (roxygen 8.1.0; the diff of `man/`
  touches that topic only), `vignettes/brms-migration.Rmd` (table row
  and the "Within-group residual correlation" paragraphs),
  `vignettes/compatibility.Rmd` (one sentence), `README.md`, `SPEC.md`
  (the stale "remaining deferrals: ar()/ma()" sentence).
- `tests/testthat/test-brms-likelihood.R`: row 18d was an exemption
  asserting the refusal; it is now a log_prob identity against brms's
  compiled program, and the header's exemption count says so.
- Tests: new `tests/testthat/test-autocor-cond.R`; the refusal block of
  `test-autocor.R` flipped to assert the two forms differ; a new block
  in `extensions/frmtmb.sample/tests/testthat/test-simulators.R`.

## 3. Validation

Script: `dev/arcov-validate.R`, output `dev/arcov-log/validate.txt`
(lane build, this machine). The reference side is an R transliteration
of the model block above, fed brms's OWN Stan data
(`brms::make_standata()`: Y in brms's order and brms's `J_lag`); it
never calls frmtmb's recursion or ordering. The data are ragged (8 rows
removed from 12 x 8, several interior gaps) and shuffled. Each model is
compared at the ML optimum and at a seeded perturbation of the full
parameter vector (random effects and imputed values included), using
the joint density `obj$env$f`, so no Laplace step is involved.

```
                          case     point  loglik_frmtmb loglik_brms_recursion          rel_diff
                gaussian ar(1)   optimum -119.886196324        -119.886196324 2.37072409517e-16
                gaussian ar(1) perturbed -120.363911909        -120.363911909 1.18065743210e-16
                gaussian ma(1)   optimum -120.099812664        -120.099812664 1.18325369540e-16
                gaussian ma(1) perturbed -120.653823099        -120.653823099 0.00000000000e+00
            gaussian arma(1,1)   optimum -118.745380073        -118.745380073 0.00000000000e+00
            gaussian arma(1,1) perturbed -120.537410601        -120.537410601 4.71583208709e-16
            gaussian arma(2,2)   optimum -118.581642704        -118.581642704 3.59520775505e-16
            gaussian arma(2,2) perturbed -118.887024193        -118.887024193 1.19532428469e-16
                gaussian ar(3)   optimum -118.632944340        -118.632944340 0.00000000000e+00
                gaussian ar(3) perturbed -119.121091884        -119.121091884 1.19297552519e-16
       gaussian ma(2), no time   optimum -119.812080933        -119.812080933 2.37219061792e-16
       gaussian ma(2), no time perturbed -120.322984595        -120.322984595 1.18105902733e-16
 gaussian arma(1,1), sigma ~ x   optimum -118.739564352        -118.739564352 1.19680872948e-16
 gaussian arma(1,1), sigma ~ x perturbed -119.572042874        -119.572042874 1.18847636736e-16
           gaussian(log) ar(1)   optimum -119.470776827        -119.470776827 2.37896749191e-16
           gaussian(log) ar(1) perturbed -120.747078752        -120.747078752 3.53073255157e-16
             student arma(1,1)   optimum -118.745380083        -118.745380083 9.57400091205e-16
             student arma(1,1) perturbed -119.009233911        -119.009233911 8.35867770405e-16
     student arma(1,1), nu ~ x   optimum -118.745380080        -118.745380080 5.98375057016e-16
     student arma(1,1), nu ~ x perturbed -119.050703785        -119.050703785 1.19368086566e-16
                 weights ma(1)   optimum -155.871611685        -155.871611685 0.00000000000e+00
                 weights ma(1) perturbed -156.633172817        -156.633172817 5.44361879145e-16
                cens arma(1,1)   optimum -106.740174347        -106.740174347 0.00000000000e+00
                cens arma(1,1) perturbed -109.390096610        -109.390096610 2.59819767157e-16
                   trunc ma(1)   optimum -120.096932453        -120.096932453 1.18328207265e-16
                   trunc ma(1) perturbed -120.650342923        -120.650342923 2.35570896376e-16
        (1 | subj) + arma(1,1)   optimum -104.338308792        -104.338308792 5.44799120464e-16
        (1 | subj) + arma(1,1) perturbed -105.820399404        -105.820399404 1.34292204483e-16
      mi() response, arma(1,1)   optimum -117.424138675        -117.424138675 2.42043158681e-16
      mi() response, arma(1,1) perturbed -118.894024207        -118.894024207 3.58576172603e-16
         rescor, ar(1) on both   optimum -233.857071680        -233.857071680 3.64603591752e-16
         rescor, ar(1) on both perturbed -234.763728201        -234.763728201 1.21065164743e-16
max rel_diff: 9.57e-16 
```

```
 model    par         frmtmb      arima_css          rel_diff
 AR(1)   mean 2.013701089591 2.013700634030 2.26230551843e-07
 AR(1)    ar1 0.620546209033 0.620546494256 4.59631509525e-07
 AR(1) sigma2 0.996013394657 0.996013420185 2.56302338636e-08
 MA(1)   mean 1.012341635144 1.012341512162 1.21482713097e-07
 MA(1)    ma1 0.480731353238 0.480730802991 1.14460517776e-06
 MA(1) sigma2 0.984477686672 0.984477692159 5.57420156315e-09
              frmtmb    arima_resid          rel_diff
AR(1) -2832.46557774 -2832.46557774 1.60548235594e-15
MA(1) -2822.23302619 -2822.23302619 6.44521337063e-16
```

Against brms's COMPILED Stan program (Stan compiles on this machine
since the round update): `dev/arcov-brms-stan.R`, output
`dev/arcov-log/brms-stan.txt`. It uses test-brms-likelihood.R's own
translator (`stan_pars_from_fit()`, flat priors) and prints the two
quantities `brms_lp_check()` asserts: `rstan::log_prob()` at frmtmb's
estimates against `logLik()` (the joint density for the random-effect
model), and the largest gradient entry there. Ragged groups with gaps,
shuffled rows, N = 194.

```
                      model  stan_log_prob         frmtmb          rel_diff      max_abs_grad
                      ar(1) -211.376605269 -211.376605269 0.00000000000e+00 1.99401469999e-04
                      ma(1) -219.202592179 -219.202592179 5.18638199446e-16 3.87822477875e-04
                  arma(2,1) -207.574956842 -207.574956842 2.73845264022e-16 6.31588541179e-04
          student arma(1,1) -206.973994375 -206.973994375 2.19712313259e-15 1.80936014802e-05
              weights ma(1) -270.932142064 -270.932142064 6.29420545245e-16 1.53283642229e-04
             cens arma(1,1) -184.939793239 -184.939793239 3.07361751981e-16 2.77863396977e-04
                trunc ma(1) -219.201736657 -219.201736657 6.48300279548e-16 8.05489622063e-04
 (1 | g) + arma(1,1), joint -217.554623953 -217.554623953 2.61283432308e-16 8.88178419700e-16
```

The same identity is now row 18d of `tests/testthat/test-brms-likelihood.R`
(ar(1), ma(1), arma(2,1) on ragged, shuffled data), where it used to be
an exemption asserting the refusal.

The `stats::arima(method = "CSS")` rows compare one series of 2000.
Both optimizers are run tight (nlminb `rel.tol = 1e-14`, which makes
nlminb print "singular convergence (7)" and "false convergence (8)" in
the log: it cannot make that tolerance, and the objective check below
is what the comparison rests on).
The objectives agree to 1e-15 at arima's own estimates (frmtmb's
log-likelihood there against the gaussian log-likelihood of arima's
CSS residuals). The estimates agree to 1e-7 to 1e-6 relative, which is
optimizer tolerance on a flat optimum, with both optimizers run tight.
For AR(1), CSS drops the first row while brms keeps it with no lagged
term; weight 0 on row 1 removes its density and keeps its residual as
the lag of row 2, which is exactly CSS. A pure MA conditions on nothing
in either, so no weight is needed. ARMA(1,1) is NOT arima's CSS: arima
conditions on the first p rows with `e_p = 0`, brms sets `e_1 = r_1`,
so only the brms-recursion check covers it.

Speed (same script, N = 2000 in 100 groups of 20): per `fn + gr`
evaluation of the taped objective, blocks >= 1.2 s, arms interleaved,
minimum of 5 rounds, with the same object timed twice as a control.

```
         model us_per_fn_plus_gr ratio_to_control
       control             66.64            1.000
 control_again             68.52            1.028
            ar             85.08            1.277
            ma             89.53            1.343
          arma            102.89            1.544
 arma_cov_true            373.13            5.599
 ma_one_series             91.33            1.370
```

The MA recursion runs over the within-group POSITION (max group length
iterations of vectorized ops, 20 here); pure AR is p gathers of one
vector. The worst shape for the recursion, one series of 2000 rows
(2000 positions of one row each), costs about the same per evaluation,
because the tape holds O(N (p + q)) scalar operations either way; the
position loop matters for taping time and R overhead, not for the
evaluation. `cov = TRUE` ARMA on the same data is about 4 to 6 times
the plain model.

## 4. Tests run (lane build unless marked)

One file per process, `NOT_CRAN=true`, lane library first. Final
runs (`dev/arcov-log/tests-final.txt`):

```
RESULT frmtmb test-autocor-cond.R pass=46 fail=0 err=0 skip=0
RESULT frmtmb test-autocor.R pass=192 fail=0 err=0 skip=0
RESULT frmtmb test-brms-likelihood.R pass=17 fail=0 err=0 skip=34
RESULT frmtmb test-brms-suite-brm.R pass=23 fail=0 err=0 skip=0 (FRMTMB_BRMS_FIT_TESTS=true)
RESULT frmtmb test-brms-suite-standata.R pass=87 fail=0 err=0 skip=0 (FRMTMB_BRMS_FIT_TESTS=true)
RESULT frmtmb test-brms-suite-priors.R pass=36 fail=0 err=0 skip=0 (FRMTMB_BRMS_FIT_TESTS=true)
RESULT frmtmb.sample test-simulators.R pass=54 fail=0 err=0 skip=0
```

Earlier in the lane, on the same R code (`dev/arcov-log/tests.txt`):
test-compat.R pass=576, test-priors-autocor-classes.R pass=62,
test-simulate-newdata.R pass=44, test-adefects.R pass=85,
test-portability.R pass=95, test-prior-mv-resp.R pass=36,
test-prior-location-dpar.R pass=55, test-importance.R pass=226,
test-conditions.R pass=150, test-conditions-census.R pass=11,
test-method-residue.R pass=64, all fail=0 err=0 skip=0. The one failure
in that batch was test-brms-likelihood.R's row 18d exemption block,
which asserted the old refusal; it is now the identity test below.

Gated, with Stan (`/opt/rlib/stan` first, `FRMTMB_BRMS_FIT_TESTS=true`),
only the blocks this lane changes, through `test_file(desc =)`
(`/tmp/lanes/arcov/arcov-row18.R`): "row 18d: brms's default cov = FALSE
ARMA is an identity" pass=12 fail=0; "row 18: ar(p = 1), cosy and unstr
residual correlation" pass=30 fail=0. The rest of the gated
test-brms-likelihood.R (about 30 Stan compiles) was not run: no other
row reaches an autocorrelation block, and the objective's only change
outside `cov = FALSE` blocks is the branch condition.

Pins seen failing:

- `test-autocor-cond.R` against the base build: `pass=2 fail=0 err=10`,
  every fitting block erroring on the old refusal
  (`dev/arcov-log/base-pin.txt`). That is a refusal, so two behavioral
  mutants were also run against the lane build
  (`dev/arcov-log/mutants.txt`): the AR part put inside `err`
  (`err[[t]] <- res[rows] - sma - sar`) fails 9 assertions in 4 blocks;
  sorting a group by data order instead of time order fails 15 in 3
  blocks.
- The frmtmb.sample block against base frmtmb.sample with the lane's
  frmtmb: `fail=2` in "a cov = FALSE term: fresh recursion, or brms's
  one-step draw" (`dev/arcov-log/sample-pin.txt`).

## 5. brms-suite ports

brm:110 (`brm(y ~ ma(x), dat, poisson())`, brms "Please set cov =
TRUE") now holds as brms wrote it: frmtmb's refusal quotes brms's
sentence. The harness reported "STALE OWN-WORDS: brms's assertion holds
as written"; the row moves from own-words to plain pass and out of
`dev/brmsport-verdicts-own.tsv`. The defect its note filed ("frmtmb
gives this message for ma(x) on every family, including gaussian") is
fixed. brm:106 and brm:108 keep their own-words passes; their notes said
the grammar check "runs before the cov = TRUE one", which no longer
exists. standata:310, :315, :316 stay "cannot transfer" (brms's
`old_order` attribute); their reasons said the setup was refused, and
the setups now run.

The generator: the verdict files are edited by `dev/arcov-verdicts.R`,
which applies, idempotently, exactly the rows `dev/brmsport-ledger.R`
would write (the ledger builder itself still needs a full recording of
every topic). The generated files were then produced by
`/tmp/lanes/shared/brmsport-gen.R` (the round's copy pointed at the
brms 2.23.0 mirror); its output was identical to a hand edit made
before the generator was available, and `git diff` of the generated
files touches only brm:106, :108, :110 and standata:310, :315, :316.
The generator also copies `tests/testthat/helper-brms-suite.R` into
frmtmb.sample, which carries the fixture note below. The recordings are
`dev/arcov-log/rec-brm.tsv` and `rec-standata.tsv`, from the lane build;
the tier's own `dev/brmsport-log/rec-frmtmb-{brm,standata}.tsv` take
those rows (the three assertions and the three standata setup lines,
which now run).

Not changed: brms fixture 1 (`helper-brms-suite.R`) still adds
`cov = TRUE` to brms's `arma(visit, patient)`. It can now be brms's
formula exactly, but fixture 1 feeds 94 assertions in the methods file
and 8 in emmeans, all verdicted on the covariance form; switching it is
a re-verdict of that tier, not a row flip. Its `changed` note now says
so instead of "frmtmb has only the residual-covariance ARMA".

## 6. Decisions and what is NOT done

- Non-gaussian families under `cov = FALSE` are refused. For an MA part
  this is brms's own refusal. For an AR part brms fits N latent
  residuals; that model is not implemented and the refusal says what
  brms fits and names the random-effect spelling of a latent AR
  process. It is feasible under Laplace (an N-dimensional random
  effect), but it is a different model from the one the gap names.
- `newdata` for `fitted()` / `frm_linpred()` / `predict()` requires the
  response and refuses without it (brms would draw the missing response
  recursively, which a point estimate cannot do); `simulate(newdata =)`
  is the forecasting route. frmtmb.sample's `posterior_predict(newdata
  = )` and `re_formula =` still refuse as "structured", as they do for
  `cov = TRUE`.
- `se.fit` / `Est.Error` of the one-step mean is a central-difference
  delta method over the outer parameters and the group effects
  (`fit_fd_se()`), one level at a time: the shift of a row reads other
  rows' residuals, possibly of other levels, so the batched attribution
  `fitted()` uses elsewhere would be wrong.
- `residuals(type = "osa")`, `autocor_matrix()`, quadrature,
  importance, `vcov_cluster()`, frmtmb.sample `log_lik()` / `loo()`:
  refused as before, with messages that are true of `cov = FALSE`.
  `log_lik()` is the most useful of these to add later: the likelihood
  is a product of rows given the observed past, which is what brms's
  `log_lik()` returns.
- The `max 300 time levels` cap and the non-consecutive-level warning
  do not apply to `cov = FALSE` (no matrix is built; lags are rows).
  A new refusal: a lag order no group is long enough to reach.

## 7. Found, not fixed

- On the base build, a fit whose bound is active warns "Large maximum
  absolute gradient" (16.7 for `class = "b"` with `ub = 0.1`, 32.9 for
  `cov = TRUE` `class = "ar"` with `ub = 0.1`;
  `/tmp/lanes/arcov/arcov-bound.R`): the convergence check does not
  exclude parameters sitting on a bound. Pre-existing, not specific to
  this lane.
- At N = 2000 the default nlminb tolerance stops with gradients of
  1.3e-3 to 2.2e-3 and warns, for `cov = TRUE` as well as `cov = FALSE`
  (`/tmp/lanes/arcov/arcov-grad.R`). Also pre-existing.

## NEWS entry

- `ar()`, `ma()` and `arma()` now fit brms's default `cov = FALSE`
  form, which was refused before: `mu` gains a regression on the
  group's earlier residuals and each row keeps the family's density,
  conditional on each group's first rows. It matches brms 2.23.0's
  likelihood to 1e-15 relative (gaussian and student, with `weights()`,
  `cens()`, `trunc()`, `mi()`, random effects and `rescor = TRUE`) and
  `stats::arima(method = "CSS")` on one series. The coefficients are
  unconstrained, as in brms, and priors and bounds of class `ar`/`ma`
  act on them directly. `fitted()`, `predict()` and `residuals()` use
  brms's one-step mean; `simulate()` runs the recursion over its own
  draws. Other families are refused in brms's words. `brm:110` of the
  ported brms suite now passes as brms wrote it.
