# frmtmb.sample (development version)

Needs the development version of frmtmb, which exports
`predict_rescor_draw()` and gives `rescor_row_loglik()` its `resp`.

## Breaking changes

* `posterior_predict()` of a `set_rescor(TRUE)` model draws the
  responses of a row jointly, as brms does: per draw, a multivariate
  normal with covariance `diag(sigma) R diag(sigma)`, or for `student()`
  the multivariate t with the shared `nu`. Before, no call returned that
  draw: without `resp` the method returned the first response only, and
  draws of two responses taken one call each correlated at -0.05 to
  0.08 on a fit whose residual correlation was 0.82. On a fit of three
  responses (200 rows, 2000 draws), the correlations of the predictive
  residuals are now 0.750, 0.554 and 0.315, within 1.3 Monte Carlo
  standard errors of what the draws imply and within 0.7 of brms's own
  draw from the same parameter draws. Under `sigma ~ x` and under
  `student()`, the same holds within 1.4 and 0.8. A `sigma` that varies
  by row is drawn at that row's value.

* Without `resp`, `posterior_predict()`, `posterior_epred()`,
  `posterior_linpred()` and `predictive_error()` of a multivariate
  model return every response, as brms does: a draws x observations x
  responses array with brms's response names. `fitted()`, `predict()`
  and `residuals()` summarize it into brms's observations x 4 x
  responses array. Before, they returned the first response only.
  Without `rescor`, each layer is the draw that `resp =` gives for its
  response.

* `resp` takes several responses, in the order you name them, in
  `posterior_predict()`, `posterior_epred()`, `posterior_linpred()`,
  `predictive_error()`, `fitted()`, `predict()`, `residuals()`,
  `log_lik()`, `loo()`, `waic()` and `loo_predict()`, `loo_epred()`,
  `loo_linpred()` and `loo_predictive_interval()`. Each takes brms's
  spelling of a name (`ya` for `y_a`) as well as frmtmb's, refuses a
  name given twice, as brms does, and reads an empty `resp` as every
  response, as brms does. Under `set_rescor(TRUE)`, a `resp` that is
  not in the model's order is answered, with the correlations in the
  order named; brms 2.23.0 stops there, because it looks up the
  correlations by the names in the order given.

* Where a response predicts more than one value per row (the category
  probabilities of an ordinal or categorical response), the stacked
  layout has no place for it. `posterior_epred()` and `fitted()`
  without `resp` refuse such a model; brms binds the category columns
  beside the other responses there. `predict()`
  without `resp` refuses to summarize it, as brms stops; with
  `summary = FALSE` it returns the draws. Name the response with
  `resp =`. `predictive_interval()` without `resp` is refused, as brms
  refuses it.

* `log_lik(resp =)` of a `set_rescor(TRUE)` model is brms's: one
  response gives its own density, and several give their joint
  density. Before, every `resp` gave the joint density of all the
  responses, so `loo(resp = "y1")` on three responses reported -724.7
  where the density of y1 gives -280.8. `loo()`, `waic()`,
  `loo_predict()`, `loo_R2()` and the other `loo_*()` functions read
  the corrected matrix through their `resp`. `log_lik(combine = FALSE)`
  returns the joint matrix, as in brms, where it was refused.

* One response of a `student()` model with `set_rescor(TRUE)` is drawn,
  and its density taken, as a t with the shared `nu`, which is the
  margin of the multivariate t. brms 2.23.0 reads `nu = 1` there (its
  one-response preparation finds no `nu` of its own), which is not the
  margin of its own model. This is a deliberate divergence from brms.

# frmtmb.sample 0.24.0

Needs frmtmb 0.77.0. That release reads brms's `(flat)` as a prior
density and gives `build_objective()` its `prior_only` field. Its
`sigma` is the innovation sd under the covariance form of `ar()`,
`ma()` and `arma()`.

## Breaking changes

* Under the covariance form of `ar()`, `ma()` and `arma()`, the draws of
  `sigma` are of the innovation sd, as brms's are. Before, they were of
  the marginal sd. On an AR(1) with `ar` near 0.92, the posterior mean
  of `sigma` was 1.154 against brms's 0.432. It is now 0.433.
  `summary()`, `posterior_summary()`, `variables()` and the draws read
  the same parameter as brms's.

* The default `student_t(3, 0, s)` on `sigma` and a prior of your own
  are densities on the innovation sd itself. A `lb` or `ub` on `sigma`
  is taken, as in brms. Under `sigma ~ x`, the default and a prior of
  your own on the intercept of `sigma` are read at the log innovation
  sd, as brms's `Intercept_sigma`. With `normal(-1, 0.5)` there, the
  posterior mean of `ar` was 7.95 Monte Carlo standard errors from
  brms's and is now 0.55.

* We compared the draws with brms's own (4 chains of 4000 on each side):
  `ar()`, `ma()` and `arma()` under the defaults, `ar()` under user
  priors, and `ar()` under `sigma ~ x`. Every posterior mean is within
  1.8 Monte Carlo standard errors of brms's. Every sd ratio is within
  0.969 to 1.037. The bulk effective sample size of `sigma` is 2.0, 0.94
  and 1.2 times what it was under `ar()`, `ma()` and `arma()`.

* `conditional_effects(method = "posterior_predict")` on draws draws its
  band at the draws' `sigma`. Under these terms that is now the
  innovation sd. brms does the same, because its `conditional_effects()`
  drops every autocorrelation term. The band is narrower than before by
  `sqrt(1 - ar^2)` for an AR(1).

* Draws saved with frmtmb.sample 0.23.0 for a model with one of these
  terms still report the old marginal sd as `sigma`. `summary()`,
  `posterior_summary()`, `as.matrix()` and `variables()` show it. Their
  `log_lik()`, `loo()` and `posterior_predict()` are refused. Sample
  such a model again.

* The message of the default priors no longer marks `sigma` with "on
  brms's innovation sd": it is brms's sigma.

* `prior = "flat"` is brms's `(flat)` in every slot: no density on the
  parameter brms declares, inside the bounds brms declares. That is the
  model brms fits with `""` in every row of its prior table; brms has no
  opt-out of its priors. Until now the opt-out added no density at all,
  which is flat on frmtmb's internal parameters: on a log standard
  deviation, improper at `sd = 0`, and on a first-order correlation,
  `(1 - rho^2)^(-3/2)`. A standard deviation, `sigma` and the other
  positive distributional parameters are now flat on their own scale. A
  correlation matrix (`cor`, `rescor`, `cortime`, `lncor`) is flat on
  its Cholesky factor, and the ordered thresholds and `delta` are flat
  on themselves. The defaults that are already flat on brms's parameter
  stay: `ar`, `ma`, `cosy`, `lagsar`, and `dirichlet(1)` on a simplex, a
  mixture's weights included. Against brms 2.23.0 with `""` in every
  slot, 4 chains of 3000 draws each. On `y ~ 1 + (1 | g)` with 6 groups
  the posterior mean of `sd_g__Intercept` is 1.458 in brms, 1.455 now
  and was 1.197. On `y ~ x + (1 + x | g)` with 8 groups the
  correlation's mean is -0.283 in brms, -0.265 now, and was -1.000 (sd
  below 0.0001, 529 divergent transitions). On `ar(cov = TRUE)` with 6
  series of 5 the mean of `ar[1]` is 0.481 in brms, 0.480 now and was
  0.541. A random-effect block now non-centers under `prior = "flat"`,
  since flat on the standard deviation closes the tail at `sd = 0`, and
  `prior_summary()` lists the `(flat)` rows.

* A bounds-only specification that replaces a default,
  `set_prior("", class = "sd", lb = 0.1)`, is brms's `(flat)` on that
  slot: flat on the standard deviation above 0.1, where it was flat on
  its log. So is `set_prior("", class = ...)` with no bound, and a table
  row with an empty prior: on `cor` it is flat on the Cholesky factor,
  on the intercept no density. The `prior_summary()` of a
  `prior = "flat"` run lists the intercepts as `(flat)` too. Passed back
  as `prior =`, as the object or as its table through `as.brmsprior()`,
  it gives the same density: log density and gradient to 1e-12 relative
  at 40 points.

* The named list of prior objects
  (`prior = list(theta_1 = prior_normal(0, 1))`) is refused, as in
  frmtmb. It took over the internal parameters it named after the class
  priors. A name inside a correlation block removed the default `lkj(1)`
  from the whole block: `list(theta_4 = prior_normal(0, 1))` on
  `(1 + x + z | g)` left two correlations with no density. The draws'
  `prior_summary()` no longer carries an `overrides` attribute.

* frmtmb.spatial's `car()` dependence parameter is the one parameter
  brms has a class for that is not yet brms's `(flat)`: brms declares it
  on (0, 1) and leaves it flat. Here it carries no density, which is
  flat on its logit and not integrable at 0. Under `prior = "flat"` such
  a posterior is improper; give it a prior with class `"theta"` on the
  coef `default_prior()` lists.

* `prior =` refuses a `default_prior()` or `validate_prior()` table of
  `route = "fit"`, the default route, edited, cut down to some rows or
  columns, or not. It refuses a priorlist made from one with
  `as.brmsprior()`, `+` or `c()` too. Its `(flat)` rows are `frm()`'s
  defaults. Read as brms's `(flat)`, they replaced brms's defaults on
  the intercept, sigma, the standard deviations and the correlations
  without a word. The message points to `route = "sample"`, whose table
  is brms's defaults, and to `prior = "flat"`. `as.data.frame()` drops a
  table's route. So `as.brmsprior()` marks a data frame that holds
  frmtmb's printed `(flat)`. The priorlist is refused where its `(flat)`
  rows would replace one of brms's proper defaults. A `route = "sample"`
  table never does that. A table built by hand in brms's spelling, with
  `""` for `(flat)`, is read as written, and so is brms's own.

* `delta`, the distance between equidistant ordinal thresholds, is flat
  on `delta` by default, brms's `(flat)` with `lb = 0`. frmtmb holds
  `log(delta)` for an ordered family, and the default was flat on that.

* `prior_summary()` and the default-prior message print brms's `(flat)`
  with `lb` and `ub` for `cosy`, for `lagsar` and for a first-order
  `cov = TRUE` `ar` or `ma`. They printed `uniform(-1, 1)`. The density
  is unchanged.

* `set_prior(class = "theta")` on a standard deviation, a correlation or
  `rescor` is refused, as in frmtmb (`?set_prior`). One such prior used
  to retire the default `lkj(1)` on every correlation of its block.

* **The evidence ratio of a point hypothesis is brms's.** `hypothesis()`
  divides the kernel density of the posterior draws at the point by the
  kernel density of the prior draws there, as brms's `evidence_ratio()`
  does. It returns the prior draws of each hypothesis as
  `prior_samples`. Without prior draws for every parameter a hypothesis
  names, its `Evid.Ratio` is `NA` with no message, as in brms. Before,
  the denominator was the prior density read off the prior
  specification. So draws made without prior draws gave a number where
  brms gives `NA`. And `sigma`, a standard deviation, a correlation, a
  nonlinear function and the intercept were refused with a warning. brms
  gives a ratio for all but the intercept. We compared with brms 2.23.0
  on `y ~ x + (1 | g)` (4 chains of 6000 a side, seed 11). The ratios of
  `x = 0`, `x = 0.2`, `sigma = 1` and `sd_g__Intercept = 0.5` are within
  1.7 Monte Carlo standard errors of brms's.

* `as_draws*(inc_warmup = TRUE)` and `as.mcmc(inc_warmup = TRUE)` return
  the warmup instead of refusing, and `update()` of draws samples again
  instead of refusing.

## New features

* `frm_sample()` samples a `mo()` term on two values, as brms 2.23.0
  does (`Jmo = 1`). The draws report the weight `simo_<term>1[1]` as a
  constant 1, as brms's do. On 200 rows, `y ~ mo(x) + z` with `x` in 1
  and 2 agrees with brms's own Stan program (2 chains of 2000 each):
  every posterior mean within 1.6 Monte Carlo standard errors, every
  posterior sd within 3.3 percent.

* `frm_sample()` refuses a `mo()` term on one value. That term is 0 on
  every row, so its coefficient has an improper posterior under the
  default flat prior. brms 2.23.0 cannot sample it either: its Stan
  program rejects `Jmo = 0`.

* `frm_sample(sample_prior = )` takes brms's `"no"` (the default),
  `"yes"` and `"only"`, and `TRUE` and `FALSE`.
  * `"yes"` adds draws from the priors to the draws, under brms's names,
    as `brm(sample_prior = "yes")` does. That is one column per prior
    statement of brms's Stan program for the model: `prior_b` for a
    class-wide prior, `prior_b_x` for each coefficient of a vector with
    a coefficient-specific prior, `prior_Intercept`, `prior_sigma`,
    `prior_sd_g`, `prior_cor_g`, `prior_simo_mox1[k]`, and so on. A flat
    prior has no column. A bounded parameter is drawn truncated to its
    bounds. We checked 19 models. They cover population-level and
    monotonic coefficients, intercepts of the location and of a
    distributional parameter, nonlinear parameters, `sigma`, `nu`, `phi`
    and `hu`. They also cover smooth standard deviations in two
    predictors, group-level standard deviations, correlations, `rescor`
    and `mo()` simplexes. And they cover mixture weights, ordinal
    thresholds, `ar` of two responses, and multivariate and bounded
    priors. On each, the names are brms 2.23.0's, in brms's order. Each
    column is drawn directly from its prior, independent of the
    likelihood, as brms draws it in generated quantities. On 39 columns
    of 9 models the one-sample Kolmogorov-Smirnov p-value against the
    exact (truncated) prior is between 0.125 and 0.986. Against brms's
    own draws of two models (14 columns, 4000 draws a side), it is
    between 0.033 and 0.962. An ordinal model with a threshold prior of
    its own gets brms's `prior_Intercept__1`, `prior_Intercept__2`, ...,
    one per threshold. A `gp()` term's `sdgp` and `lscale` get no prior
    draws, and neither do brms's `meanme` and `sdme`, which frmtmb has
    no prior slot for.
  * `"only"` samples the prior of the whole model with the likelihood
    left out, as brms does. So `posterior_predict()`, `pp_check()` and
    `conditional_effects()` give the prior predictive distribution.
    Every parameter needs a proper prior: a flat one is refused with
    brms's message ("Sampling from priors is not possible as some
    parameters have no proper priors. Error occurred for parameter
    'b'."). Unlike brms, a vector with some coefficients flat is refused
    too. So is a model with missing values of an `mi()` term, which have
    no density without the likelihood. A `mi(sdy)` response keeps its
    measurement density, as in brms. A model with `me()` terms is
    refused, naming brms's classes `meanme` and `sdme`: brms needs
    proper priors there, and frmtmb has no prior slot for them yet. We
    checked `y ~ x + (1 | g)` with a `normal(0, 1)` prior on `b` (4
    chains of 4000). We compared each mean and each 5, 50 and 95 percent
    quantile of `b_Intercept`, `b_x`, `sigma`, `sd_g__Intercept`,
    `r_g[1,Intercept]` and of a prior predictive draw. Each is within
    3.0 Monte Carlo standard errors of brms's. At two more seeds every
    quantile of `b_x`, `sigma` and `sd_g__Intercept` is within 1.3 of
    the exact prior's (brms's: 2.4).

* New `prior_draws()`, brms's, with brms's `variable` and its deprecated
  `pars`; brms's deprecated `prior_samples()` reaches it while brms is
  loaded. A name gets the prior column whose name is the longest match
  at its start, in random order (`b_x` reads `prior_b`). The overall
  intercept `b_Intercept` is left out with brms's warning.

* The warmup iterations rstan saves (`save_warmup = TRUE`, the default
  of rstan and of brms) stay in the draws' stanfit, as in a brmsfit.
  `inc_warmup = TRUE` on `as_draws()`, `as_draws_array()`,
  `as_draws_df()`, `as_draws_list()`, `as_draws_matrix()`,
  `as_draws_rvars()`, `as.mcmc()`, `as.matrix()`, `as.array()`,
  `as.data.frame()`, `rhat()` and `neff_ratio()` converts them to the
  draws' variables. It puts them before each chain's draws, as brms
  does. `as.mcmc()` refuses it with `combine_chains = TRUE` and on
  thinned draws, with brms's messages, and every accessor refuses a
  warmup that no longer fits the draws. Every other method reads the
  post-warmup draws alone, as before.

* `update()` of draws samples the same model again with changed
  arguments of `frm_sample()`, as brms's `update()` does for
  `sample_prior`, `prior` and the sampler's options:
  `update(ds, sample_prior = "yes")` adds prior draws. A new formula or
  new data is refused, as before.

# frmtmb.sample 0.23.0

Needs frmtmb 0.76.0, which exports `autocor_flat_priors()`,
`autocor_bound_chart()`, `autocor_innov_sd()` and
`autocor_row_loglik()` and codes `mo()` from its minimum.

## Breaking changes

* With frmtmb 0.76.0, which codes an integer `mo()` variable from its
  minimum as brms 2.23.0 does, `frm_sample()` samples brms's model on
  such data. On 300 rows with values 1 to 4, the draws have the 3
  weights `simo_mok1[1]` to `simo_mok1[3]`, where they had 4. Every
  posterior mean agrees with brms's within 2 Monte Carlo standard
  errors, where the Intercept was 12 of them away. The code of this
  package did not change; its tests need frmtmb 0.76.0.

* **The default priors of a residual-correlation term are brms's.**
  Before, an `ar()`, `ma()` or `cosy()` coefficient was flat on frmtmb's
  internal parameter. For a first-order `cov = TRUE` `ar` or `ma` that
  is `(1 - rho^2)^(-3/2)` on `rho`, with its mass at -1 and 1. For
  `cosy` it is a density flat on the logit of its whole window. And the
  default on `sigma` was on frmtmb's marginal sd, where brms's `sigma`
  under a first-order `cov = TRUE` `ar()`, `ma()` or `arma()` is the
  innovation sd. Now each coefficient is flat on itself, as brms's
  `(flat)` is: on (-1, 1) for `ar` and `ma` under `cov = TRUE`, on
  brms's [0, 1] for `cosy`. The default `student_t()` on `sigma` is read
  on brms's innovation sd, `sigma sqrt(1 - ar^2)`,
  `sigma / sqrt(1 + ma^2)` or
  `sigma sqrt((1 - ar^2) / (1 + ma^2 + 2 ar ma))`, with the Jacobian of
  that factor. The reference is brms 2.23.0's own draws, 4 chains of
  4000 on each side (seed 3). Against them, the means of `ar[1]` were
  42.3 Monte Carlo standard errors apart, on an AR(1) of 0.9 in 8 groups
  of 6 time points. For `sigma` they were 19.6 apart. An MA(1) of 0.6
  had every `ma[1]` draw at 1. Now every parameter of the ar, ma and
  arma cases is within 1.6 standard errors of brms's, the sds within 4
  percent. For `cosy` (0.03 in 30 groups of 4), 48 percent of the draws
  were below 0, where brms has none. The mean was 0.011 against brms's
  0.076. Now there are none, and the means are 0.075 to 0.078 against
  brms's 0.075 to 0.078 on three seeds. A `cov = TRUE` term of order two
  or more, which brms does not have, is flat on its stationary or
  invertible region. Its `sigma` default stays on the marginal sd. The
  `cov = FALSE` coefficients, brms's default form, are flat as before,
  as in brms. `frm_sample()`'s message and
  `default_prior(route = "sample")` list the new rows.

* **`sigma` under a first-order `cov = TRUE` `ar()`, `ma()` or `arma()`
  is the marginal sd. Every prior on class `sigma` is read at the
  innovation sd, as brms's is.** The draws' `sigma` column, `summary()`
  and `variables()` report frmtmb's marginal sd: on the AR(1) above, a
  posterior mean of 1.16. brms's `sigma` is the innovation sd (0.43);
  the model is the same. A prior of your own on class `"sigma"` is read
  where the default is, at the innovation sd. So a
  `default_prior(route = "sample")` table passed back as `prior =`, or
  the default's text typed, gives the default's density exactly. Before,
  a user prior there was read on the marginal sd. That density was up to
  0.72 nats from brms's with the same text (three points on the AR(1)
  data). `frm()` with the same prior is the mode of the density
  `frm_sample()` samples. A `lb` above 0 or a `ub` on `sigma` there is
  refused.

* **A prior of your own keeps the bound brms declares for its slot.**
  brms keeps a declared bound under a user prior: `real<lower=1> nu`,
  `real<lower=0> sigma` under an identity link,
  `real<lower=0,upper=min_Y> ndt`, `real<lower=0,upper=1> cosy`. Before,
  a user prior of the same slot dropped the default's bound: with
  `set_prior("normal(2, 2)", class = "nu")` 3.55 percent of the draws
  were below 1. Now the bound stays on each side the prior leaves unset,
  and a `lb` or `ub` of your own replaces that side, as in brms. With
  `set_prior("normal(0, 0.1)", class = "cosy")` the posterior mean is
  0.054 against brms's 0.053 (z 1.5); it was 0.009, with draws below 0.
  A bounds-only prior on a residual-correlation class keeps the flat
  density as well: `set_prior("", class = "ar", lb = 0.5)` is flat on
  [0.5, 1] in brms. Here it was flat on frmtmb's internal parameter
  (mean 0.946 against brms's 0.916, z 43.8); now the mean is 0.916 (z
  0.5).

* A prior on `class = "theta", coef = "thetaac_1"` is refused where a
  class reaches that parameter (`ar`, `ma`, `cosy`, `cortime`); see
  frmtmb's NEWS. The message for the retired `lower`/`upper` arguments
  names `coef = "theta_1"` as its example.

## New features

* `log_lik()` of draws no longer refuses a residual covariance term
  (`ar(cov = TRUE)`, `ma()` and `arma()` with `cov = TRUE`, `cosy()`,
  `unstr()`, `fcor()`, frmtmb.spatial's `sar()` and other registered
  terms). A column is brms 2.23.0's leave-one-out conditional density of
  the row given the other rows of its group, normal or Student-t. It
  equals brms's `log_lik()` at the same draws to 2.0e-15 relative on 18
  cases. `loo()`, `waic()`, `psis()` and the `loo_*()` functions take
  these models. The values do not sum to the log-likelihood, as in brms.
  On a group with a missing time point, the conditional is the fitted
  model's. Its `ar()`, `ma()` and `arma()` lags count time levels, where
  brms's count positions.

## Bug fixes

* A bounded residual-correlation parameter is sampled through one
  logistic onto its interval, as Stan's own bound is in brms. Stan's log
  bound on frmtmb's internal parameter, stacked on the term's own map,
  gave 5 to 9 divergent transitions per 4 chains of 2000. That was on a
  `cosy` near its bound at 0 (seeds 3 to 5, mode and random starts),
  where brms has none. Now there are none.

* When the ML estimate of `cosy` is below brms's bound at 0, the warning
  says that the bound is the default, brms's declaration. It names
  `set_prior("", class = "cosy", lb = -1/(d - 1))` to sample the whole
  window `frm()` fits. It said "the requested bound(s) on thetaac_1".

# frmtmb.sample 0.22.0

Needs frmtmb 0.75.0, which exports `frmtmb_covstruct_info()` and
`frmtmb_autocor_pars()` and gives coefficients brms's names.

## Breaking changes

* The draws of a model with an `exp()` term carry brms's names:
  `b_expx`, `sd_g__expx` and `cor_g__Intercept__expx`, as brms 2.23.0
  names them, where they were `b_.frm_expx`, `sd_g__.frm_expx` and
  `cor_g__Intercept__.frm_expx`. `fixef()` of the draws names the row
  `expx`, and `hypothesis(draws, "expx > 0")` works. The change is in
  frmtmb, which builds the names, and the new test of it needs frmtmb
  0.75.0. The same holds for a nonlinear parameter of a multivariate
  model, `b_yn_a_Intercept` where the draws held `b_a_yn_Intercept`.

* `fixef()` of the draws lists the rows in brms's order, which is the
  order of `fixef()` of the fit: every centered intercept first. For
  `bf(y ~ x, sigma ~ exp(z))` brms 2.23.0 gives `Intercept`,
  `sigma_Intercept`, `x`, `sigma_expz`. The draws' `fixef()` listed
  `Intercept`, `x`, `sigma_Intercept`, `sigma_expz`, the order of the
  draws' columns, which is unchanged.

## New features

* The draws of a model with frmtmb.spatial's `sar()` carry brms's
  `lagsar` or `errorsar` column, derived from the sampled internal
  parameter, as `ar[1]` and `cosy` are. The names come from core's
  registry (`frmtmb::frmtmb_autocor_pars()`), so a term another package
  registers gets its own column the same way.

* The default prior on `lagsar` and `errorsar` is brms's: flat on `rho`
  over its interval, which `frm_sample()` announces as
  `uniform(lower, upper)`, the ends shown to 6 digits. The same holds
  for a bounded parameter of a term another package registers.
  `set_prior(class = "lagsar")` and `class = "errorsar"` replace it.

## Bug fixes

* The derived `sd_` and `cor_` columns asked core's
  `frmtmb_covstruct_info()` which structures are fields, where they kept
  their own list of five names. A field another package registers,
  beside a block on the same grouping factor, was mapped onto that
  block's names. It failed the map's own check and sent every draw
  through the slow path: on one model, 0.70 ms a draw against 0.053 ms
  on the fast path now. The values were right either way. Needs frmtmb
  0.75.0, which exports `frmtmb_covstruct_info()`.

# frmtmb.sample 0.21.1

Needs frmtmb 0.74.0. The package itself loads on frmtmb 0.73.0, but
four of its test files need what frmtmb 0.74.0 and frmtmb.latent 0.7.0
add: the covariance structure registry, `frmtmb_family(simplex_pars =)`
and `mixture_mvn()` in frmtmb.latent.

## Bug fixes

* The note of default priors from `frm_sample()` names every simplex
  that a family declares with `frmtmb_family(simplex_pars =)`, by its
  brms name. It named only the `sbhaz` simplex of `cox()`, so the
  `dirichlet(1)` default on the simplex of another family was not
  announced. For `cox()` the note now reads "sbhaz: dirichlet(1), brms's
  default, on the simplex the cox() family declares".

## Other changes

* The test of `mixture_mvn()` draws reads the family from frmtmb.latent,
  where it moved, and skips without it. The help of
  `posterior_predict()` and `pp_mixture()` links it there.

* frmtmb.survival is a suggested package. The tests of `cox()` draws
  read `cox()` from it, and skip without it.

# frmtmb.sample 0.21.0

Needs frmtmb 0.73.0. It registers its draws methods on the `loo_*()`,
`kfold()`, `reloo()`, `loo_moment_match()` and `loo_subsample()`
generics that frmtmb now defines, and re-exports them. It calls the
sampling-API export `check_weights_unbounded()`.

## Breaking changes

* **On a subset of the draws (`ndraws`, `draw_ids`), `loo()`, `psis()`
  and the weights of every `loo_*()` function take the relative
  efficiency on one chain.** brms 2.23.0 does the same, and so did the
  `loo_*` types of `pp_check()`. They used `r_eff = 1`. The same draws
  are now weighted the same way everywhere. On 60-row gaussian and
  random-intercept models with 4000 draws, `loo(ndraws =)` at 500, 1000
  and 3000 draws moved by at most 0.0065 in `elpd_loo`. Its standard
  error was 5.0. On every draw nothing changes, and `frm_sample(thin =)`
  draws keep their chains.

* **`posterior_average(weights = "stacking")` and `"pseudobma"` take
  their weights from `loo_model_weights()`**: one `log_lik()` per model
  and its relative efficiency, as brms's `model_weights()` does. Before,
  they read a `loo` criterion that `add_criterion()` stored, which can
  come from other draws. Without a stored criterion the weights are
  identical to the old ones.

* The `kfold()`, `reloo()`, `loo_moment_match()` and `loo_subsample()`
  generics, and those of the six `loo_*()` functions under New features,
  are frmtmb's, re-exported here. On a maximum-likelihood fit each
  refuses, with frmtmb alone attached too: the six name `frm_sample()`,
  the other four what to use instead. `loo_moment_match()` and
  `loo_subsample()` had no method for a fit, even with this package
  loaded.

* **`frm_sample()` from a formula refuses a posterior that does not
  integrate.** Suppose that every row with a positive weight has a
  random effect of its own in `mu`. The likelihood grows like
  `sigma^-(sum(w) - n)` as `sigma` shrinks to 0. A prior that is finite
  at 0, as brms's default `student_t(3, 0, s)` on `sigma` is, then
  leaves the posterior improper once `sum(w) - n` reaches 1. Such a run
  used to return chains that had not mixed (largest Rhat 2.5 on 2 chains
  of 1000, `sigma` near 1e-4). The check is core's, against the priors
  the chains sample under. A prior on `sigma` that vanishes at 0 fast
  enough makes the posterior proper but not informative:
  `lognormal(0, 1)` gave a posterior median `sigma` of 7.2e-05 and a
  largest Rhat of 2.83. A usable model has weights that sum to the
  number of rows or fewer, or no observation-level term.

## New features

* **brms's leave-one-out functions on draws**: `loo_epred()`,
  `loo_linpred()`, `loo_predict()`, `loo_predictive_interval()`,
  `loo_R2()` and `loo_model_weights()`, with brms 2.23.0's arguments,
  defaults and output shapes. A multivariate model without `resp` gives
  one layer per response, named as brms names them, and a category
  margin stays third. Each computes one `log_lik()` and one
  `loo::psis()` per call, on the same draws as the prediction (per
  response for `loo_R2()`, per model for `loo_model_weights()`). The
  weights are those of `psis()` and of `loo(save_psis = TRUE)`.
  `psis_object =` reuses weights and computes no `log_lik()`. On the
  same posterior draws as a brms 2.23.0 fit, every `loo_epred()` and
  `loo_linpred()` value is identical to brms's under R's reference BLAS.
  So is every `loo_model_weights()` weight (a stacking pair of 0.155 and
  0.845). The `loo_R2()` bootstrap draws at one seed agree to rounding
  (largest relative difference 1.1e-14). Under OpenBLAS, where brms's
  and frmtmb's predictors round apart in the last bits, they agree to
  that rounding.

* `add_criterion(x, "loo_R2")` stores the `loo_R2()` bootstrap draws, as
  brms stores them, and `loo_R2()` reads them. It was refused.
  `add_criterion()` gives each criterion the arguments it takes, as brms
  does: `add_criterion(x, c("loo", "loo_R2"), seed = 1)` seeds the
  bootstrap and no longer stops in `loo()`.

* `psis()` honors `draw_ids`, which it ignored, and refuses an argument
  it does not read, which `loo::psis()` swallowed.

* `loo()` and the `loo_*()` functions accept brms's default values by
  name (`allow_new_levels = FALSE`, `sample_new_levels = "uncertainty"`,
  `newdata = NULL`), which `loo()` refused, and `loo()` and
  `loo_model_weights()` take `draw_ids`.

* `resp` takes brms's spelling (`ya` for a column `y_a`) in all six
  functions, and the layers of a multivariate result carry it.
  `loo_epred()`, `loo_linpred()`, `loo_predict()`,
  `loo_predictive_interval()`, `loo_model_weights()` and `psis()` also
  take frmtmb's `y_a`, which `log_lik()` and `posterior_epred()` take.
  `loo_R2()` takes brms's spelling only, as `bayes_R2()` and brms do.

* **`brmsfamily("negbinomial2")` and `dirichlet_multinomial()` fits
  sample under brms 2.23.0's default priors**: `sigma` of negbinomial2
  gets `gamma(0.4, 0.3)` from brms's family hook (where negbinomial's
  `shape` gets `inv_gamma(0.4, 0.3)`), and `phi` gets
  `gamma(0.01, 0.01)`. On the same data the posterior means agree with
  brms's to 0.05 posterior sd, with and without a group term.

* **A `dirichlet_multinomial()` fit's draws take the `multinomial()`
  paths**: `posterior_linpred()` with no `dpar` stacks the `mu<cat>`
  predictors into brms's `draws x observations x (K - 1)` array.
  `posterior_epred()` and `posterior_predict()` return
  `draws x observations x K` arrays that read `trials()` from `newdata`.
  `predictive_error()` and `pp_check(type = "error_binned")` refuse the
  family as polytomous.

## Divergences from brms

* `newdata`, `re_formula` and the new-level arguments, given a value
  other than brms's default, are refused by the `loo_*()` functions, as
  `loo()` refuses them: a weight comes from the observation's own column
  of `log_lik()`, so it exists only for a row the model was fitted on.
  brms computes the weights on the new rows' own response.

* `ndraws` and `draw_ids` choose the draws once, for the weights and the
  prediction. brms draws a separate random subsample for each, so its
  weights belong to other draws: on a 60-row gaussian model with
  `ndraws = 100` the two subsamples shared 2 draws.

* An argument the `loo_*()` functions have no use for is an error that
  names it. `loo_linpred(scale = "response")`, which brms 2.23.0 reads
  nowhere (its result is the default's, on a logit link too), is one.

* `loo_R2()` refuses an ordinal or categorical family, as `bayes_R2()`
  does. brms refuses the categorical families and computes the ordinal
  ones with a warning.

# frmtmb.sample 0.20.0

Needs frmtmb 0.72.0. It calls the sampling-API export `ord_lp_blocks()`
and writes `uniform()` and `lognormal()` priors. It reads the
per-response `lncor_<resp>` names and the per-coefficient `sd` priors of
that release.

## Breaking changes

* **`posterior_linpred()` with no `dpar` on a family with several
  location predictors returns brms's `draws x observations x P` array**,
  every margin unnamed: every `mu<cat>` of `categorical()`,
  `multinomial()` and the simplex families (`K - 1` layers with a
  reference category, `K` without), and every component's `mu<k>` of a
  mixture. It returned the first predictor's matrix. On these families
  `transform = TRUE` without a `dpar` is refused in brms's words, a
  two-category `categorical()` with its one predictor included.

* **A default steps aside only for a prior on its own slot**, as in
  brms. A narrower prior keeps the defaults on the parameters it does
  not reach. `set_prior(class = "sd", group = "g")` replaces the default
  on the standard deviations of `g`. It keeps `student_t(3, 0, 2.5)` on
  the other groups. Adding `coef = "x"` narrows the prior to `sd_g__x`.
  Before, the defaults of the whole class stepped aside, and the other
  standard deviations were flat.

* **A class-wide prior does not replace a default written per
  coefficient.** For example, `set_prior("inv_gamma(2, 1)", class =
  "lscale")` on a `gp(x)` model keeps brms's tuned default on
  `lscale_gpx`, as brms keeps it. Before, the class-wide prior replaced
  it. The default-prior message now names a specification that applies
  to no parameter. That happens for this reason, or when a more specific
  specification of your own covers its parameters. The message names the
  coef to write.

* **`frm_sample()` puts brms 2.23.0's default prior on every
  distributional parameter that brms gives one.** They were flat on the
  link scale: `shape` (`gamma(0.01, 0.01)`, and `inv_gamma(0.4, 0.3)`
  for the negative binomials), `phi` (`gamma(0.01, 0.01)`), `nu`
  (`gamma(2, 0.1)`) and `kappa`. So were the exgaussian `beta`, `alpha`,
  `zi`, `hu`, `zoi`, `coi` and `quantile` (`beta(1, 1)`), and
  frmtmb.eam's `bs`, `bias` and `ndt`. A parameter with a formula has
  brms's default on its intercept, on the link scale. On small data the
  posteriors move. A beta `phi` near 300 on 10 rows moved by 30 Monte
  Carlo standard errors. A negative binomial `shape` and a Student-t
  `nu` had no proper posterior before.

* **The ordinal thresholds take brms's `student_t(3, 0, 2.5)`**, class
  `Intercept`, and a mixture's weights without a formula take brms's
  `dirichlet(1)` on the simplex. Both were flat.

* **`sigma ~ 1` takes brms's default on sigma's intercept.** It got a
  natural-scale default that addressed a slot the model does not have.
  `sigma` without a formula now reads the response on the identity
  scale, as brms does. Under a log `mu` link it read it on the log
  scale.

* An `unstr()` residual correlation takes brms's `lkj(1)`, class
  `cortime`.

* **`frm_sample(seed = )` fixes every chain.** The init jitter of chains
  2 and up came from the session's random stream, so two calls with
  `seed = 3` differed unless `set.seed()` came first. The jitter is now
  drawn from `seed` and no longer from the session's stream. rstan's
  sampler still advances the session's stream, as under `brm()`.

* A natural-scale default carries brms's bound on its class where
  frmtmb's link does not hold the parameter inside it. A log-linked
  `ndt` is bounded above by the smallest response. An identity- or
  log-linked `nu` is bounded below by 1. brms declares them so, and Stan
  then transforms the draws inside rather than meeting the density's
  edge. Where the link holds `ndt` inside a bound of its own,
  frmtmb.eam's, the uniform is over the link's image and no bound is
  written: `uniform(0, min_Y)` for an ungrouped `wiener()`, `lba()` or
  `rdm()`, and `uniform(0, 1)` for the fraction `ndt_group()` makes of
  `ndt`.

* **`posterior_predict(newdata =)` reads every addition term from
  newdata.** It read only `trunc()`, `rate()` and `thres(gr = )`, and
  gave the simulator no other term, from the method's first version on.
  A binomial with `y | trials(N)` and newdata `N = 1000` drew 0 or 1
  (column mean 0.595 where `posterior_epred()` gave 576.3). A gaussian
  with `se(s)` drew at sd 0.033 where newdata had `s = 5`. `predict()`,
  `predictive_error()` and `residuals()` of draws, which call it, change
  with it. Newdata without the `trials()` or `se()` variable is now
  refused, naming the term, as brms refuses it.

## New features

* Draws of a `logistic_normal()` response inside `mvbf()` report its
  correlations as `lncor_<resp>__b__c`, with `<resp>` brms's Stan name
  of the response. The sampler's default puts brms's `lkj(1)` on each
  such response.

* **The draws carry brms's `sd_<group>__<coef>` and
  `cor_<group>__<coef1>__<coef2>`**, and a smooth's `sds_`, `rescor__`,
  `ar[1]`, `cosy` and `cortime__`, as columns of the draws matrix. So
  `as_draws_df(draws)$sd_g__Intercept` works as on a brms fit. They are
  derived from each draw's `theta_` columns, which stay. The names are
  core's, the ones `variables()` of the fit lists.

## Bug fixes

* `posterior_predict(newdata =)` and `posterior_predict(re_formula =)`
  hand the simulator each distributional parameter on its natural scale.
  They handed it the reporting scale, which for a mixture's `theta` is
  already the softmax probability. The simulator applied the softmax
  again, so the mixing weights were wrong. The defect is in every
  release from 0.2.0 (frmtmb 0.52.0, whose afd00a43 gave a mixture's
  weights their softmax reporting scale). With frmtmb's fix of the same
  release, a mixture also draws each class at its own rows' `trials()`.

# frmtmb.sample 0.19.0

Needs frmtmb 0.71.0: it reads the sampling-API exports `ce_method()`,
`ce_predict_family_check()`, `ce_sim_one()` and `pp_check_resp()`,
`ce_plan_eval(one = )`, and the `gen_extreme_value()` chart
`gev_xi_chart_frame()`.

## Breaking changes

* **`loo()` and `waic()` of draws take brms's arguments in brms's
  order**, `(x, ..., compare, resp, ndraws, model_names)`. Draws given
  without a name are further models, which are compared as brms
  compares them. So `ndraws` must be named, as it is in brms.

* **`loo()` and `waic()` of draws refuse by name the arguments they
  dropped.** These are `newdata`, `re_formula` and the other arguments
  `log_lik()` refuses, which brms reads. `loo(newdata = )` returned the
  in-sample elpd in silence (lane remeasure, R2). `moment_match = TRUE`
  and `reloo = TRUE` refit, so they are refused too. `pointwise` and
  `k_threshold` change nothing and are accepted.

## New features

* **The draws of core's `dirichlet()` and `logistic_normal()` fits take
  brms's shapes.** `posterior_epred()` on a `dirichlet()` fit is the
  `draws x observations x categories` array of category means.
  `posterior_predict()` draws rows on the simplex for both families,
  with layers named by category, as brms's are.

* The correlations of a `logistic_normal()` fit are the draws column
  `lncor__b__c`, brms's name. They take brms's `lkj(1)` default, and
  each `sigma<k>` takes the default that brms gives `sigma`.

* As in brms, `posterior_epred()` refuses `logistic_normal()`.
  `pp_check()` and `predictive_error()` refuse both families.

* `frm_sample()` samples a `gen_extreme_value()` `xi` without a formula
  as brms samples it, through `tmp_xi`. It maps `tmp_xi` into the
  support that each draw's `mu` and `sigma` leave, and gives it brms's
  default `normal(0, 2.5)`. On its `log1p` scale the support's edge was a wall in the
  sampled space. The review fitted n = 40 data sets at `xi` -0.45,
  0.20 and -0.85 (4 chains x 2000). There 1187, 263 and 3856 of 4000
  transitions diverged, with Rhat 2.30 on the last, against 0, 1 and 0
  in brms.
  Through `tmp_xi` none diverge, and the posterior mean of `xi` is
  within 0.005 of brms's at the first two. One difference
  from brms remains: frmtmb holds `xi` above -1 on every route, where
  brms's range goes on (`dev/famuni-out/tmpxi.txt`, `tmpxi-sample.txt`).
  The hold is material where the data put mass there. On the xi -0.85
  data brms puts 29 percent of its posterior below -1. There
  `frm_sample()`'s posterior is brms's conditioned on `xi > -1` (median
  of `xi` -0.876 against -0.877; the review's round 2). A prior on
  class `xi` is then a prior on `tmp_xi`, as in brms. A `cens()` or
  `trunc()` response, a second GEV response, or GEV components of a
  `mixture()` keep the `log1p` scale, and `frm_sample()` warns that they
  do.

* `frm_sample()` puts brms's default `normal(0, 4)` on the intercept of
  a predicted `xi`, `xi ~ 1` included (`normal(0, 2.5)` under
  `link_xi = "identity"`).

* **`loo()` and `waic()` of several draws objects return brms's
  comparison.** It holds the criterion of each model (`loos`),
  `loo::loo_compare()` of them (`diffs`) and brms's pairwise
  `ic_diffs__`, printed as brms prints its `loolist`. `compare =` and
  `model_names =` are brms's. Models with different numbers of
  observations are refused in brms's words. Models of different
  responses get loo's warning, through brms's `yhash` attribute.

* **`logLik()` of draws is `log_lik()`**, the draws by observations
  matrix, as brms's `logLik()` of a fit is. It had no method.

* **`bayes_R2()` of draws takes `newdata` and `re_formula`**, which
  brms hands to `posterior_epred()`. With `newdata` it is the
  R-squared of newdata's rows, against newdata's own response column.
  A newdata without that column is refused in brms's words.

* **`model.frame()` of draws is the model frame of their fit.** stats'
  default rebuilt a frame from the formula's environment, which holds
  none of the data. It stopped at "object 'count' not found" on brms's
  ported row `brmsfit-methods:120`.

* **`conditional_effects(method = "predict")` on draws** (brms's
  `"posterior_predict"`) draws one response per grid row and posterior
  draw, at the parameters of that draw. That is the posterior
  predictive band brms draws. Under the default `robust = TRUE`,
  `estimate__` is its median and `se__` its MAD. A category-valued
  response is refused, as the fit method refuses it.

* **`transform =` applies to the predicted draws**, in
  `conditional_effects(method = "predict")` and `posterior_predict()`,
  with brms's warning that the argument is deprecated.

* **`frm_sample()` of a `frm_multiple()` result samples the fit of every
  imputation and combines the draws**, as brms's `brm_multiple()`
  combines its chains. `nchains()` is the number of imputations times
  `chains`. A `seed` gives imputation `j` the seed `seed + j - 1`.
  Methods that read the data in sample read the first imputed data set,
  as brms's do.

* **The diagnostics of combined imputations are brms's.** `summary()`
  leaves out `Rhat`, `Bulk_ESS` and `Tail_ESS`, as brms's does for a
  `brm_multiple()` fit. `loo()` takes a relative efficiency of 1, as
  brms does. Chains of different imputations sample different
  posteriors, so their R-hat measures the imputations.

## Bug fixes

* **`pp_check()` of multivariate draws refuses a missing, plural or
  unknown `resp` in brms 2.23.0's words**, as the fit method now does.
  A missing one was "pp_check() is not supported yet for multivariate
  fits".

* **`summary(draws, waic = TRUE)` is refused by name, naming
  `waic(x)`.** Old brms printed WAIC there, and brms 2.23.0 ignores the
  argument. It got the generic "has no argument `waic`". `waic = FALSE`
  is accepted.

# frmtmb.sample 0.18.0

Needs frmtmb 0.70.0: it reads `simplex_frame_terms()` and the
Dirichlet prior entries frmtmb now resolves, and a `cox()` fit's
baseline through `mo_chart_frame()`.

## Breaking changes

* **A `cox()` baseline is sampled as brms samples it**: in softmax
  coordinates under brms's default `dirichlet(1)` (class `sbhaz`), and
  the draws carry the weights as `sbhaz[k]`, brms's names, in place of
  the `sbhaz_raw_k` coordinates. With no density on the coordinates two
  of them ran to -3e11 and -9e12 on a right-censored model, with 239
  divergences and R-hat 3.07; now the weights' means and sds agree
  with brms 2.23.0's within 2.1 Monte Carlo standard errors on two data
  sets, with 0 divergences (`dev/moparity-cox-sample.R`).

* **An exact `gp()` field is sampled non-centered**, as brms's
  `sdgp * L zgp`, when its sd and length scales have priors (the
  defaults do). The centered field is a funnel: on lane surface's
  `y ~ gp(x)` the formula route ran 1286 of 2000 transitions into the
  maximum tree depth with R-hat 1.52 and a bulk ESS of 7; non-centered,
  none, R-hat at most 1.016, bulk ESS 312 or more, in a third of the
  time. The factor is a Cholesky factorization of the kernel on the
  tape, `n^3 / 3` per gradient for `n` positions.
  `reparameterize = FALSE` keeps the centered field.

## New features

* **A `simo` or `sbhaz` prior row replaces `dirichlet(1)` on its
  simplex**, with the softmax's Jacobian. With `dirichlet(c(5, 1, 1))`
  the draws of `simo_` agree with brms 2.23.0's on a weakly and a
  strongly identified simplex: the means within 1.71 Monte Carlo
  standard errors and the sds within 2.51 (`dev/moparity-dir-sample.R`).

* **`mo()` crossed with a factor and group-level `mo()` terms sample**,
  with brms's names (`simo_mox:fb1[k]`, `bsp_mox:fb`, `r_g[, mox]`). On
  `mo(x) * f2`, with a two-level `f2`, and on `mo(x) + (mo(x) | g)` the
  means and sds agree with brms 2.23.0's within 2.12 Monte Carlo
  standard errors (`dev/moparity-dir-sample.R`).

* **`as_tmbstan()` names only the simplexes with no density on the
  fit's tape**: a MAP fit with a `simo` row has its Dirichlet there.
  The function still samples the fit's own tape, as documented; it
  adds no prior.

## Bug fixes

* **`parnames()` on draws warns once whichever generic dispatched it**,
  as on a fit.

# frmtmb.sample 0.17.0

Needs frmtmb 0.69.0: it imports and re-exports frmtmb's `stancode()`,
`standata()`, `pp_mixture()` and `add_criterion()` generics, and it
holds a `mo()` simplex through frmtmb's `mo_simplex()`, `mo_coords()`,
`mo_chart_frame()`, `mo_frame_terms()` and `summary_mo_frame()`.

## Breaking changes

* **`stancode()`, `standata()` and `pp_mixture()` are frmtmb's
  generics**, re-exported here; this package keeps its methods for
  draws. frmtmb defines them because it now has methods for a fit.

* **`pp_mixture()`'s summary names its components `P(K = k | Y)` and
  numbers its rows**, as brms does and as frmtmb's method for a fit
  does. They were `class1`, `class2`. On draws of a model that is not a
  mixture it stops in brms's words, "Method 'pp_mixture' can only be
  applied to mixture models.", as on a fit.

* **`frm_sample()` puts brms's default `dirichlet(1)` on every `mo()`
  simplex and reports brms's `simo_` weights.** The sampler put no
  density on a simplex, so the fit's coordinates were the prior. In
  frmtmb's old softmax coordinates that prior is `1 / prod(w)`, and a
  weakly identified simplex sat at a vertex (44 divergences); in the
  coordinates frmtmb now fits in it is not integrable at the
  barycenter, and the draws sat there with sd 0. The simplex is now
  sampled in softmax coordinates under `dirichlet(1)`, which is the
  flat density on the simplex (so `prior = "flat"` keeps it), and the
  draws carry the weights as `simo_<term>1[k]`, brms's names, in place
  of the `zeta<j>_k` coordinates. Against brms 2.23.0 on a weakly
  identified simplex and on brms_monotonic's `fit1` (4 chains of 2000
  each side), every weight's mean and sd agree within 1.7 Monte Carlo
  standard errors, with 0 divergences and R-hat at most 1.0033.
  Models without `mo()` terms sample exactly as before. A `simo` row in
  `prior =` is still refused: `set_prior()` has no Dirichlet density.

* With more than one `mo()` term, each term's weights are adjacent
  columns, term by term, as brms orders them.

* `check_laplace()` on a `mo()` fit compares each simplex on the
  weight scale: the ML weights and their delta-method errors (as
  `summary()` gives them) against the posterior of the `simo_`
  weights. The coordinates are not comparable, since the fit's sheet
  and the draws' need not be the same. The simplex is sampled under
  `dirichlet(1)` there too, the flat density on the simplex, under
  which the weights' posterior is the fit's likelihood; a MAP fit
  keeps its own prior. A weight at 0 or 1 is reported in a message,
  with no Wald error.

* `as_tmbstan()` says in a message, on a fit with `mo()` terms, that
  the simplex has no density on that route (it is sampled flat in the
  fit's own coordinates, which a weakly identified simplex collapses
  in) and that `frm_sample()` samples it under `dirichlet(1)`.

## New features

* **`plot()` of draws is brms's display**: a histogram and a trace per
  parameter through `bayesplot::mcmc_combo()`, `nvariables` to a page,
  with brms's `pars`, `variable`, `regex`, `fixed`, `combo`, `bins`,
  `theme`, `plot`, `ask` and `newpage`, and `N` as brms's deprecated
  alias with brms's warning. It refused every call before.

* **`add_criterion()` on draws** stores `loo`, `waic` and `bayes_R2` in
  `x$criteria` as brms does, and `loo()` and `waic()` with no further
  argument return the stored object, brms's `use_stored`.
  `overwrite = TRUE` clears the stored criterion before it computes, as
  brms does, so it recomputes. `kfold`,
  `loo_subsample`, `loo_R2` and `marglik` are refused with their own
  reasons.

## Bug fixes

* `posterior_samples()` on draws warns once that it is deprecated.
  While brms was loaded it warned twice, once from brms's generic and
  once from the method.

* The compatibility pre-flight of `frm_sample()` no longer loads R's
  base-priority packages to ask whether a name is one of their
  functions; it reads them without loading. On a headless Linux
  machine, loading tcltk warned "no DISPLAY variable so Tk is not
  available" into an unrelated call.

* `check_laplace()` leaves a parameter whose `z_shift` or `sd_ratio` is
  not finite out of its "Laplace/Wald approximation questionable"
  message. The test passed `NA` to `if (any(...))`, so such a parameter
  stopped the call ("missing value where TRUE/FALSE needed") when no
  other parameter was flagged, and was listed as `NA` when one was. The
  returned table is unchanged.

# frmtmb.sample 0.16.0

Needs frmtmb 0.68.0: for `ord_thres_linpred()`; for
`gp_brms_natural()`, the prior classes `"sdgp"` and `"lscale"` and the
kriging draw; and for ordinal mixtures, `hurdle_cumulative()` with
`thres(gr = )` and `cs()`, and the names of the draws of a model with
no location column.

## Breaking changes

* **Through frmtmb, a smooth's null-space draws (`bs_sx_1`) are on
  brms's scale**, and a prior with `coef = "sx_1"` is about brms's parameter:
  frmtmb now builds `s()` with brms's `diagonal.penalty = TRUE` basis.

* **A `gp()` term's draws carry brms's names on brms's scales:**
  `sdgp_<term>` and `lscale_<term>` where they were `theta_1` and
  `theta_2` on the internal log scale, as brms names them (`sdgp_gpxfa`
  for level `a` of `gp(x, by = f)`). The lengthscale is on the inputs
  brms divides by their largest distance.

* **The default priors include brms's for a GP:** `student_t(3, 0, s)`
  on class `"sdgp"` per prefix, and on class `"lscale"` an inverse
  gamma per sub-GP, tuned to that sub-GP's own distances as brms's
  `def_lscale_prior()` tunes it. `get_prior(route = "sample")` lists
  them; on nine `gp()` designs every density row equals brms
  `default_prior()`'s (63 rows, `dev/gpby-priors.R`). A Hilbert-space
  `gp(k = )` block, whose lengthscales the defaults left flat, is now
  non-centered.

## New features

* `posterior_linpred(incl_thres = TRUE)` returns brms's linear
  predictor with the thresholds included on every ordinal family: a
  draws by observations by thresholds array with brms's dimnames,
  layer `k` being `disc * (thres_k - mu)` for `cumulative()` and
  `sratio()` and `disc * (mu - thres_k)` for `cratio()` and `acat()`,
  under every threshold structure, with `cs()` terms, and `NA` past the
  thresholds of a `thres(gr = )` level that has fewer. As in brms it is
  ignored with `dpar`, `nlpar` or `transform = TRUE` and on a family
  that is not ordinal. `hurdle_cumulative()` is refused, a divergence
  the user decided on 2026-10-06: brms returns the hurdle probability
  beside the threshold predictors times `1 - hu` there, which is the
  predictor of nothing (`dev/upstream-bugs.md`, brms-21). An ordinal
  mixture is refused in brms's words, "'incl_thres' is not supported
  for mixture models."

* **Draws of an ordinal mixture** carry brms's names, each component's
  thresholds `b_mu1_Intercept[k]` and `b_mu2_Intercept[k]`. Under
  `order = "mu"` both come from the one shared vector, which is read
  back through the first component that reports it uncentered.
  `posterior_epred()` and `log_lik()` agree with the theta-weighted sum
  of brms's R-side category probabilities at the stored columns to the
  double epsilon.

* **brms's default `normal(0, 1)` on the intercept of a modeled
  `disc<k>`** of an ordinal mixture, as on `disc`.

* Through frmtmb 0.68.0, draws of a `cs()` model on `cumulative()`, as
  brms samples it: the per-threshold coefficients are
  `bcs_x[k]`, and a draw's `posterior_epred()` equals brms's own
  R-side density at that draw to 1.9e-16 (`dev/rel068-cs-probe.R`).
  A row whose offsets cross two thresholds has `NaN` probabilities and
  an `NA` draw, as in core.

## Bug fixes

* **A prediction at a position an exact `gp()` did not see is a draw of
  the field there**, from its conditional law given the draw's own
  values at the fitted positions, as brms's `posterior_epred()` draws
  it. It was that law's mean, so `posterior_epred(newdata = )` past
  the data had no spread from the field itself. At the same 1500
  parameter draws of a brms fit, the spread of the predicted field at
  six new positions is 0.978 to 1.001 of brms's, where the mean alone
  reached 0.861 to 1.000 (`dev/gpby-brms-epred2.R`, final build in
  `dev/gpby-p1-nugget-brms-1e-6.txt`). The law carries frmtmb's
  nugget, `1e-6` of the variance, as white noise per distinct
  position, as the fit's latent values do; two rows at one position
  take one value exactly.

* The draw's factor costs under 1 ms per posterior draw at 300 new
  positions: a pivoted Cholesky factor that builds only the kernel
  columns it pivots on (rank 7 to 15 there) and stops where what is
  left is below `1e-12` of the prior variance, which joins the white
  part (`dev/gpby-p1-krigfactor.R`). On `posterior_epred()` at 300
  new positions over 300 draws, CPU time is 0.78 s, against 5.76 s for
  a dense Cholesky of the same law and 0.63 s for the previous
  release's conditional mean, on one machine in one session
  (`dev/gpby-p1-krigdraw.R`).
# frmtmb.sample 0.15.0

Needs frmtmb 0.67.0: for `ord_delta_info()`, the ordinal families'
`disc` and threshold structures and `hurdle_cumulative()`'s inverse
threshold map; for the `parnames()` generic, which it imports, and the
conditional-effects engine that draws crossed and multi-membership new
levels; and for `arma_cond_fill_dpars()`, `arma_cond_fill_epred()` and
`response_codes_newdata()`.

## New features

* **Draws of an ordinal fit with equidistant or sum-to-zero
  thresholds** carry brms's names: every threshold `b_Intercept[k]`
  and, for equidistant thresholds, `delta`. The sampler holds fewer
  parameters than that (two per vector, or one fewer than the
  thresholds); the columns past them are added after the sampled ones,
  before `lp__`, as a mixture's last weight is, and dropped again
  whenever a draw is handed back to the model.

* **brms's default `normal(0, 1)` on the intercept of an ordinal
  family's `disc`**, which is what places that intercept, since the
  likelihood cannot tell it from the scale of the thresholds. Under
  `link_disc = "identity"` brms uses `lognormal(0, 1)`, which
  `set_prior()` does not carry, and the default-prior message names the
  gap.

* `nsamples(incl_warmup = TRUE)` counts the iterations the sampler
  saved, warmup included, from the stanfit that `frm_sample()` keeps,
  as brms counts them. It was refused. Draws with no stanfit behind
  them refuse it by name. `nsamples(subset = )` takes brms's check
  that the subset fits the count.

* `posterior_samples(pars = )` returns the coefficients in brms's
  order, as brms's `variables()` lists them: the intercept of every
  distributional parameter first (`b_Intercept`, `b_sigma_Intercept`,
  `b_x`, `b_sigma_x`), a nonlinear parameter's intercept staying with
  its own coefficients (`b_a_Intercept`, `b_a_z`, `b_b_Intercept`). Without
  `pars` the columns are those of `variables()`, as before. A `pars`
  that matches nothing returns `NULL`, as in brms.

* `parnames()` is frmtmb's generic, re-exported here. With brms loaded
  it warns once, as brms's does; brms's generic and the method each
  warned before.

* `conditional_effects()` on draws answers crossed
  `(1 | g) + (1 | h) + (1 | g:h)` at an unseen combination of observed
  levels, and a new member of an `mm()` term with a `by` variable.

## Bug fixes

* A `hurdle_cumulative()` fit's draws are named `b_Intercept[k]`, as
  the other ordinal families' are. They were `tau_raw_k`.

* **`posterior_predict(newdata = )` of a model with `ar()`, `ma()` or
  `arma()` and `cov = FALSE`** was refused as a structured draw, the
  same refusal as a `cov = TRUE` block, though core's `predict()`
  answers it. It answers now, and with brms's treatment of a missing
  response: a row of newdata whose response is `NA`, or every row when
  the column is absent, is filled with a draw at that draw's
  parameters before the rows after it read its residual.
  `re_formula = ` is no longer refused on such a model either.

* **`posterior_epred(newdata = )` under `cov = FALSE` fills a missing
  response with a draw**, as brms's `posterior_epred()` does, so the
  expected-response draws on the rows after it carry the spread of the
  unobserved past. Core's `fitted()` on a maximum likelihood fit fills
  with the expected value instead, because it has no draws to carry
  the fill through.

* **`predictive_error(newdata = )` of a `bernoulli()` model** codes
  newdata's response as the fit coded its own, so a response held as
  two values other than 0 and 1 (brms's level-order coding, new in
  frmtmb) is compared with the draws on the same scale.

# frmtmb.sample 0.14.0

Needs frmtmb 0.66.0: for `subset_resp_check()` and `subset_newdata()`,
for `cs_build()`, `ce_level_plan()` and the other engine functions of
`conditional_smooths()` and `conditional_effects()`, for the ordinal
families' inverse threshold maps that name the draws, and for equated
parameters.

## Breaking changes

* **An ordinal fit's draws carry brms's names.** The thresholds are
  `b_Intercept[1]`, `b_Intercept[2]` (per group under `thres(gr = )`,
  per response in a multivariate model) and hold the thresholds, not
  `cumulative()`'s first threshold and log increments; a `cs()` column
  is `bcs_<column>[k]`. They were `tau_raw_k` and `bcs<j>_k`.
  `variables(ds)` now lists what `variables(fit)` lists. A draws object
  saved with the old names still works. Code that selected
  `tau_raw_1` from the draws must select `b_Intercept[1]` now. The
  names come from each ordinal family's inverse threshold map, which
  frmtmb 0.66.0 declares.

## Bug fixes

* **`posterior_predict()`, `posterior_epred()` and `posterior_linpred()`
  on `frm_sample(laplace = TRUE)` draws returned silently wrong
  numbers.** Such draws hold no random effects, and every draws method
  read them in the layout of full draws: `sigma` and `theta` landed in
  the random effects and `NA` past the last column. On
  `y ~ x + ar(t, g) + (1 | g)` all 4500 cells of `posterior_epred()` were
  NaN; at `newdata` the values were finite and wrong, each group effect
  read off `theta_1` or `lp__`. `fitted()`, `predict()`, `residuals()`
  and `bayes_R2()` returned the same kind of numbers; `pp_check()`,
  `predictive_interval()`, `pp_mixture()` and `hypothesis()` on `sd_`
  quantities stopped in an internal error from `quantile()` or
  bayesplot. Laplace draws are now read in their own layout. A quantity
  that reads what was integrated out (the group-level coefficients, a
  smooth's coefficients, `mi()` values) is refused, naming the function
  you called, as `log_lik()` already refused; one that does not, such as
  `re_formula = NA` on a model with group-level terms alone, is computed
  exactly and equals the full draws' answer bit for bit. The new
  "Laplace draws" section of `?frm_sample` lists both, and says why the
  missing values are not filled in (user decision, 2026-09-29).
* **`conditional_effects()` on draws under `re_formula = NULL` reads
  each grid row's groups by brms's rule** (found by lane sampfix). With
  `conditions = data.frame(g = <a level>)`, every draw's effects for a
  NEW group were drawn into the first level and predicted there: moving
  `r_g[1,]` or `r_g[2,]` by 10 left the curve where it was, and levels
  1 and 2 drew the same curve. The rule is now core's, per row and per
  group-level term: a term reads its draws where the row sets every one
  of its grouping variables to an observed level, in `conditions` or as
  an effect (`"x:g"`), and draws a new level otherwise, from each
  draw's own covariance: an unset variable, an unseen level (`"99"`),
  a row of a mixed `conditions` column, or a nested `(1 | g / h)` with
  `g` set and `h` unset. At levels 1, 2 and 3 every column and the
  spaghetti frame equal brms's at the same five draws, with a
  difference of 0 (`dev/postfit2-celevel-brms.R` in the frmtmb
  repository), and `"x:g"` equals brms's to 1.6e-8.

  On a `gr(g, by = f)` term a row reads only the term of its own `f`
  level; every `re_formula = NULL` display of such a model was refused
  before, even at an observed `g`. An `mm(g1, g2)` term's unset members
  are new levels instead of the first observed one. At the same draws
  and seed, the `gr()` band at an unseen `g` and the `mm()` band with
  members unset equal brms's `sample_new_levels = "gaussian"` bands to
  8.9e-16 (`dev/postfit2-p2-brms.R` in the frmtmb repository). Two
  different unseen members are two new levels here; brms 2.23.0 draws
  them once, an artifact of how it numbers unseen values. `mm(g1, g2)`
  beside `(1 | g1)` with `g1` unset or unseen now draws both new
  levels, where it was refused. One new level is drawn once per draw
  and shared by every panel that names it. The core NEWS has the
  details.
* `predictive_error()`, `pp_check()` and `bayes_R2()` of draws no
  longer read a missing `mi()` response as an observation. The frame
  holds the placeholder 0 there, and the predictive error was
  `0 - yrep`; it is `NA` now, as brms's is, `bayes_R2()` is `NA`, as
  brms's is, and `pp_check()` leaves the row out with brms's warning.
  A `pp_check()` loo type on such a model is refused by name, since
  `log_lik()` is not defined there; it died inside bayesplot.
* `frm_sample(laplace = TRUE)` on a model with nothing to integrate out
  died in tmbstan with "invalid argument to unary operator". It now says
  so and samples the model itself; the draws equal those of the call
  without `laplace`. On a REML fit sampled as it stands
  (`prior = "flat"`) it is refused: that objective integrates the
  coefficients out too, and the draws carried `b_Intercept` over other
  parameters' values.
* `frm_sample()` on a model with one outer parameter (`y ~ 0 + x`) died
  in rstan with "no more scalars to read": a length-one init was read as
  a scalar. Each init is now passed as a one-dimensional array.
* A draws object with `stanfit = NULL` works with every method that
  reads the chain count (`nchains()`, `summary()`, `VarCorr()`,
  `as_draws_*()`, `posterior_summary()`, `mcmc_plot()` and more), as one
  chain. `nuts_params()` and `log_posterior()` refuse it by name.

## New features

* **`pp_check()`'s `loo_pit_overlay`, `loo_pit_qq`, `loo_intervals` and
  `loo_ribbon` types work.** They failed on every model, because nothing
  built the PSIS weights bayesplot needs. They are built as brms builds
  them, `loo::psis()` on `log_lik()` of the same draws as the
  predictions, and match brms 2.23.0 bitwise on brms's own draws. Where
  `log_lik()` refuses, these types refuse and say so. A type that
  `bayesplot::available_ppc()` does not list, such as the deprecated
  `loo_pit`, is refused with brms's message.
* **`posterior_average()`**, brms's model-averaged draws: draws from
  several `frm_sample()` objects in proportion to their weights, given
  as numbers or computed by `"stacking"` (the default), `"pseudobma"`,
  `"loo"` or `"waic"`. `"kfold"` and `"bma"` are refused, as `kfold()`
  and `post_prob()` are. On the same draws, weights and seed its output
  equals brms's to 1.4e-17, with the same rows, names and attributes
  (`dev/postfit2-brms-compare.R` in the frmtmb repository).
* **`conditional_smooths()` on draws**: brms's median, MAD and
  quantiles of each smooth term's curves, and brms's `spaghetti`. On
  ten draws that a `fixed_param` brms fit shares, every column and the
  spaghetti frame equal brms's to 5.3e-15.
* **`conditional_effects()` on draws takes `spaghetti`, `surface`,
  `too_far` and `select_points`**, which it refused before. The
  spaghetti frame equals brms's to 5.4e-15 on the same ten draws, for a
  one- and a two-variable effect. On laplace draws it draws its default
  curves instead of refusing every model with a group-level term.
* `log_lik()`, `posterior_epred()`, `posterior_linpred()`,
  `posterior_predict()` and `predictive_error()` on a multivariate
  model whose responses use `subset()` need one `resp =`, as brms's
  do, and on `newdata` keep the rows where the response's subset is
  `TRUE`. `log_lik(resp = )` returns one column per row of that
  response.
* `posterior_predict(newdata = )` reads the exposure of a `rate()`
  response from `newdata`.
* `nobs()` passes brms's `resp =` to the fit's method, which counts
  that response's rows.
* The draws of a model with an equated parameter,
  `bf(y ~ x, sigma1 = "sigma2")`, carry a `sigma1` column that is a copy
  of `sigma2`, as brms's draws do, and the default priors skip the
  equated parameter, which owns no coefficient.
* `as.array()`, `as.matrix()` and `as.data.frame()` of draws take
  brms's `chain` and `iteration`, which brms hands to
  `posterior::subset_draws()`: `as.array(ds, chain = 1)` is the first
  chain. They were refused by name. Any other name is still refused.
* `conditional_effects(categorical = FALSE)` on a `hurdle_cumulative()`
  fit scores the expected category by the category codes, so the hurdle
  scores 0, as core's does. Every other ordinal family is unchanged.
* The default-prior disclosure names `kappa`, which brms gives a gamma
  prior, beside `shape`, `phi` and `nu`: `xbeta()` and `von_mises()`
  carry it.

# frmtmb.sample 0.13.0

Needs frmtmb 0.65.0, for `arma_cond_resp()` and `arma_cond_dpars()`.
With an older frmtmb the package installs and loads, and then
`log_lik()`, `loo()` and `posterior_predict()` stop at their first call
with "could not find function".

* **`log_lik()` and `loo()` now work on brms's default `cov = FALSE`
  form of `ar()`, `ma()` and `arma()`**, which was refused before. A
  row's log-density is taken at brms's one-step conditional mean, the
  regression on the group's earlier OBSERVED residuals, so a column is
  an observation and the columns multiply to the whole likelihood. It
  matches `brms::log_lik()` cell for cell at the same parameter values
  (gaussian `ar(1)` bitwise; student `arma(1,1)` and
  `arma(1,1) + (1 | g)` to 4.4e-15 absolute), and `loo()`, `waic()` and
  `psis()` follow. `loo()`'s help says what that makes the estimate:
  PSIS-LOO over these columns leaves out one CONDITIONAL density, the
  one-step-ahead density of a row given its observed past, and not the
  series. The `cov = TRUE` covariance form, `cosy()` and `unstr()` stay
  refused, because there the smallest independent unit is a whole
  group; the refusal now names the form that does factor.

* A draws object with no `stanfit` behind it no longer errors in
  `log_lik()`: `x$stanfit@sim` is an error and not a `NULL` on such an
  object, so the chain identifier could not be read at all.

# frmtmb.sample 0.12.0

Needs frmtmb 0.64.0, for `rescor_row_loglik()` and `fit_extras(resp =)`.

* **`log_lik()` on a Student-t `rescor` model used the gaussian joint
  density.** It now uses the multivariate t through frmtmb's
  `rescor_row_loglik()`; on the test fixture the old value was 38.6
  log units off. On a multivariate model with an ordinal response,
  each response reads its own thresholds.

* Draws of the models frmtmb 0.64.0 adds: `thres(gr = )` thresholds
  per group, `gr(g, by = f)` blocks merged in `ranef()`, `coef()` and
  `VarCorr()` as brms merges them, and `posterior_predict()` for
  brms's `cov = FALSE` ARMA drawn around each row's one-step mean.

* An intercept written `0 + Intercept`, or under `bf(center = FALSE)`,
  gets no default prior, as brms leaves class `"b"` flat.

* The flat-prior advice no longer reads `class = "Intercept",  = ""`
  for a univariate intercept.

# frmtmb.sample 0.11.0

* The `pp_check(type = "error_binned")` refusal on draws of a
  multinomial fit calls the response a set of counts over categories,
  not a category, as the fit method now does.
* Two test files no longer leave `Rplots.pdf` in the tests directory:
  `test-draws-methods.R`, whose device guard covered one block and not
  the `pairs()` call, and the generated `test-brms-suite-methods.R`,
  whose generator now opens a null device for the file.

* **`posterior_predict()` on draws now draws a `cs()` term.** On an
  `sratio()`, `cratio()` or `acat()` model with a category-specific
  effect it drew every row as if the term were absent, in sample and at
  `newdata`, with no warning. `predict()` and `pp_check()` on draws go
  through it and change with it. `posterior_epred()` and `log_lik()`
  were right and are unchanged. Needs the frmtmb that exports
  `cs_offsets_add()`.

* On draws, a `newdata` with a grouping level the fit did not see is
  refused by the draws methods themselves, with or without
  `allow_new_levels`, and the refusal says that
  `allow_new_levels = TRUE` is refused too. Without the flag, core's
  message reached the caller and recommended "Use
  allow_new_levels = TRUE", which the draws methods then refused. The
  same holds for `sample_new_levels` alone and for a `newdata` without
  the grouping column.

* `re_formula = NA` with `allow_new_levels = TRUE` and an unseen level
  now answers, as brms does, instead of being refused: the check asked
  about the level without the `re_formula` that drops it.

* The default priors of a categorical, multinomial, mixture or
  latent-class model carry the `dpar` of each location (`mub`, `mu1`), as
  brms writes them. They used to be written without it, which frmtmb now
  refuses, so this is what keeps those models sampling. The resolved
  defaults are the same parameters with the same densities. With no
  prior of your own, a categorical or mixture fit samples exactly as
  before. With a `dpar` prior of your own on a mixture, the default for
  that slot now steps aside instead of being written and then
  overridden, which changes the order the priors are applied in: the
  draws at a fixed seed differ (by up to 8.74 on one measured fit), and
  `prior_summary()` no longer lists the overridden defaults. The
  resp-, dpar- and nlpar-qualified defaults of the other families with
  several location dpars (`lba()`, `rdm()`, `hmm()`, `lca()`,
  `mixture_mvn()`) follow the same rule.

* The default-prior announcement names each slot with its `resp`,
  `dpar` and `nlpar`, so every line is a spelling `set_prior()` takes.
  On a multivariate model it used to print `Intercept` with no response.

* `hypothesis()`'s advice for a coefficient with no proper prior gives a
  full `set_prior()` call, with `resp`, `dpar`, `nlpar` and `coef` as
  the model needs them. It used to say `set_prior(class = "b")`, which a
  multivariate model refuses.

# frmtmb.sample 0.10.0

* **The default priors are brms's on a multivariate model**: each
  response gets its own `Intercept`, `sigma` and `sd` defaults, read off
  that response, and `set_rescor(TRUE)` gets `lkj(1)` on the residual
  correlation. Until now a multivariate model was sampled with no
  default priors at all. This moves the draws of every multivariate
  model sampled without `prior = "flat"`.

* **The default intercept location takes off the mean `offset()`**, as
  brms does: `y ~ 1 + offset(off)` with `y` near 2 and `off = 10` gets
  `student_t(3, -8, 2.5)`, where it got `student_t(3, 2, 2.5)`. A
  location is now written to 15 significant digits, as brms writes it.
  This moves the draws of a model with an offset in its location
  formula.

* **The `sd` default is written once per response, distributional
  parameter and nonlinear parameter**, to follow frmtmb's brms scoping
  of class `"sd"`. With no `prior =`, every standard deviation gets the
  density it got before. With `prior = set_prior(..., class = "sd")`,
  that prior now replaces the location parameter's default only, and
  a `phi` or `sigma` block keeps its own default, as in brms. A SMOOTH's
  smoothing standard deviation is keyed the same way, so `sigma ~ s(z)`
  gets its own default and its own `sd` row. This moves the draws of a
  distributional model sampled with a class-wide `sd` prior.

* **BREAKING: `predictive_error()` and `residuals()` on draws refuse an
  ordinal, categorical or multinomial model**, as brms refuses them.

* `pp_check()` on draws refuses `type = "error_binned"` for a category
  response in the same words the fit method uses, and ignores `resp` on
  a model with one response, which is what brms does with it. A
  multivariate model still selects with `resp` and still refuses a name
  none of its responses has.

* **A one-sided `re_formula` on draws keeps the terms it names**, in
  `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`,
  `fitted()` and `predict()`. Each draw is evaluated through
  `frmtmb::frm_linpred()`, which kept every term for any formula; the
  fix is in frmtmb's development version, so it needs that frmtmb. On
  a `(1 + x | g) + (1 | h)` fit, `re_formula = ~ (1 | g)` used to be
  off the draw's own intercept-only construction by up to 2.98; it now
  equals it. A term the fit does not have is an error that names it.
  The draws already carry each draw's sampled group effects, so the
  uncertainty at a known level was in them before and nothing is added.

# frmtmb.sample 0.9.0

* **`log_lik()` is frmtmb's generic now, re-exported here**, rather than
  a second generic of the same name defined in this package. frmtmb owns
  the name because it has a method for it: `log_lik(fit)` on a
  maximum-likelihood fit refuses and names the route to draws, as
  `loo()` and `waic()` do, where it used to be "could not find
  function". `log_lik(draws)` is unchanged, and the `frmtmb_draws`
  method is still registered on rstantools' generic as well. The method
  page is `?log_lik.frmtmb_draws`; `?frmtmb::log_lik` documents the
  generic.

* `posterior_epred()` and `posterior_predict()` on draws refuse every
  name in `...` they do not read, where they used to ignore it, as the
  other methods here already do. The new-level arguments pass through
  to the refusal described below.

* Requires frmtmb 0.61.0, for `log_lik()`, `frm_linpred()`,
  `fixef_by_dpar()`, `brms_fixef_rows()` and the brms shape helpers.

* **`fitted()`, `predict()` and `residuals()` on draws exist** (item
  2.6f). They are brms's three summarizing methods, each a summary of
  the draws method beside it: `posterior_epred()`,
  `posterior_predict()` and `predictive_error()`. Before this there
  was no method at all, so `fitted(ds)` reached
  `stats::fitted.default()`, read `ds$fitted.values` and returned
  `NULL`, a wrong answer with nothing said. `summary = FALSE` gives
  the draws, `robust` and `probs` are brms's, and `predict()` on an
  ordinal or categorical response returns brms's `P(Y = k)`
  proportions.

* **BREAKING: `nsamples()`, `posterior_samples()` and `parnames()`
  answer instead of refusing.** All three were refused as "the
  deprecated brms spelling". brms keeps all three live with a
  deprecation warning, and the refusal stopped a ported script dead
  where brms's warning does not. Each now does what brms's does and
  warns as brms warns:

  - `posterior_samples(x)` is `as.data.frame(x)`, with brms's `pars`
    (a regular expression unless `fixed = TRUE`), `add_chain`,
    `subset`, `as.matrix` and `as.array`;
  - `nsamples(x)` is `ndraws(x)`, with brms's `subset`;
  - `parnames(x)` is `variables(x)`.

  `nsamples(incl_warmup = TRUE)` is still refused: `frm_sample()`
  discards the warmup rather than storing it, so there is nothing to
  count.

* **`fixef()` on an ordinal draws object reports the thresholds**,
  `Intercept[1]`, `Intercept[2]`, then the coefficients and any `cs()`
  terms, in brms's order and on the model's own scale, which is what
  the fit's `fixef()` reports. It reported the slopes alone, because
  the sampler stores the thresholds as `tau_raw`.

* `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`,
  `fitted()`, `predict()` and `residuals()` on draws refuse
  `allow_new_levels = TRUE` when `newdata` holds a grouping level the
  fit did not see, with a message that names the function called and
  says why. That case used to fail in the design builder with "Use
  allow_new_levels = TRUE", which the caller had just done. Predicting
  an unseen level from draws is not implemented;
  `predict(fit, allow_new_levels = TRUE)` on the maximum-likelihood fit
  does it. Every other use answers as before: `allow_new_levels =
  FALSE`, brms's default, and `TRUE` with no `newdata` or with levels
  the fit saw.

* **BREAKING: `hypothesis()` on draws no longer returns an object of
  class `brmshypothesis`.** The class is `"frmtmb_hypothesis"` alone;
  the shape and every element are unchanged, and what stops working is
  `is(x, "brmshypothesis")` and `inherits(x, "brmshypothesis")` in a
  ported script. frmtmb owns `print()` and `plot()` for its own class,
  so neither changes.

* **`point_estimate` and `ndraws_point_estimate` are honored** by
  `posterior_epred()`, `posterior_predict()`, `posterior_linpred()`
  and `log_lik()`. They were accepted and ignored, so
  `posterior_epred(ds, point_estimate = "median",
  ndraws_point_estimate = 2)` returned all 25 draws where brms returns
  2. As in brms it is a PARAMETER-space operation: the draws are
  collapsed to their mean or median first and the method runs once at
  that one vector, repeated `ndraws_point_estimate` times.

* `posterior_epred()`, `posterior_linpred()` and
  `posterior_predict()` build each draw's prediction with
  `frmtmb::frm_linpred()` rather than `predict()`, which is brms's
  predictive summary in frmtmb 0.61.0.

# frmtmb.sample 0.8.0

* Requires frmtmb 0.60.0, for `frm_stop()` and the other condition helpers.
* `tmbstan (>= 1.2.1)` in Suggests, and `check_tmbstan_build()` now names
  that version as the remedy. The CI step that pinned a dated snapshot of
  the Stan trio is removed.

* Requires the frmtmb release that exports `frm_stop()`; the
  `frmtmb (>= 0.59.0)` floor must move to it.
* **BREAKING:** every error, warning and message that frmtmb.sample
  raises is classed. An error has the class
  `c("frmtmb_sample_error", "frmtmb_error", "error", "condition")`,
  and warnings and messages follow the same pattern, so
  `tryCatch(frmtmb_error = )` catches any refusal. The class vector no
  longer contains `simpleError`, `simpleWarning` or `simpleMessage`.
  See `?frmtmb::frmtmb-conditions`.
* **BREAKING:** a value that matches none of the choices of
  `loo_compare(criterion =)`, `hypothesis(scope =)`,
  `pp_check(prefix =)` and `predictive_error(method =)` is refused
  with a `frmtmb_sample_error` that names the argument, the
  value and the choices. The old text was `'arg' should be one of ...`.
  A partial value still matches.
* `posterior_interval()`, `as.matrix()`, `rhat()`, `neff_ratio()`
  and the other accessors refuse an exact variable name that the
  draws do not carry as a `frmtmb_sample_error`. The posterior
  package used to raise it as a plain error. The text keeps
  posterior's words, `missing in the draws object`.

# frmtmb.sample 0.7.0

* Requires frmtmb 0.59.0, for `brms_par_labels()`, `brms_coef_table()`,
  `brms_stan_name()`, `expand_b()`, `varcorr_layout()`, the
  `hypothesis()` helpers (`hyp_eval_in()`, `hyp_expr_vars()`) and
  `check_prior_slots()` on `?frmtmb::frmtmb-sampling-api`.

* **BREAKING: the draws carry brms's names.** `b_Intercept`, `b_x`,
  `b_sigma_Intercept` for the coefficients and `r_g[1,Intercept]` for a
  group-level coefficient (`r_gs[lev.1,Intercept]` for a level with a
  space, `r_g:h[1_p,Intercept]` for an interaction group), in
  `ds$draws`, `variables()`, every
  `as_draws_*()`, `as.matrix()`, `as.array()`, `as.data.frame()`,
  `rhat()`, `neff_ratio()`, `summary()` and `posterior_summary()`. They
  were `Intercept`, `x`, `sigma_Intercept` and `b[1]`. A covariance
  parameter brms does not sample keeps the name `confint()` gives it on
  the fit, `theta_1`; a block with no brms counterpart (reduced rank,
  smooth, GP, CAR, SPDE) keeps `b[i]`. What stops working:
  `ds$draws[, "x"]`, and anything that selected `b[` columns.

* **BREAKING: a distributional parameter nobody wrote a formula for is
  stored on its natural scale, under brms's name.** The column is
  `sigma` (or `shape`, `nu`, and `sigma_ya` for response `y_a` of a
  multivariate model) and holds sigma, as brms's does. It used to be
  `b_sigma_Intercept` holding log sigma, which brms reports only when
  `sigma ~ 1` is written out, so `fixef(ds)` had a row, `coef(ds)$g` a
  slice, and `posterior_summary(ds)` and `as_draws_df(ds)` a column that
  brms does not have. They are gone. What stops working:
  `ds$draws[, "b_sigma_Intercept"]` on such a model; use `sigma`, or
  `log(ds$draws[, "sigma"])` for the old values.

  A mixture's weights are stored as brms's simplex: `theta1 ...
  theta<K-1>` hold the mixing probabilities, and brms's `thetaK` is a
  derived column before `lp__`. They used to be the log ratio under the
  name `theta1_Intercept`. Every reader that hands a draw to the model
  maps all `K` back to the log ratios together.

* **Draws labels a level repeats are suffixed as brms suffixes them.**
  Levels `lvl 1` and `lvl.1` are both `r_gd[lvl.1,Intercept]` in brms,
  which names the later one `r_gd[lvl.1,Intercept]__1`. Without the
  suffix every accessor failed on "Duplicate variable names".
  `hypothesis()` on draws refuses a name that matches two columns
  instead of reading the first.

* **BREAKING: the accessors return brms's objects.** `as.matrix()`,
  `as.array()` and `as.data.frame()` return brms's unclassed draws
  objects and take brms's `pars`, `variable`, `draw` and `subset`, with
  brms's deprecation warnings; `as_draws_*()` take `variable`, `regex`
  and `inc_warmup`; `as_draws()` is brms's `as_draws_list()`.
  `posterior_summary(x, pars, variable, probs, robust)` and
  `posterior_interval()` cover every variable by default, as brms's do.
  Compared against brms's own installed methods on the same draws: 88 of
  90 calls `identical()`, the other 2 `rvars` objects that differ only by
  a cache environment (`dev/brmsnames-findings.md`).

* **BREAKING: `fixef()`, `ranef()`, `coef()` and `VarCorr()` are brms's.**
  `summary = FALSE` returns the raw draws, as in brms, where it used to
  return a summary; `robust`, `probs`, `pars` and `groups` are brms's.
  `coef()` broadcasts every population-level coefficient, as brms's
  does, and `VarCorr()` has brms's `sd`/`cor`/`cov`/`residual__`
  structure, with one `residual__` row per response on a multivariate
  model. `ranef()` on a model with no group-level effects is refused
  with brms's message. `ranef()` and `coef()` on a reduced-rank block,
  whose draws are factor scores rather than coefficients, compute the
  coefficients per draw; on `frm_sample(laplace = TRUE)` draws, which
  carry no group-level draws, they are refused by name rather than
  returned as `NA`.

* **BREAKING: `hypothesis()` returns brms's `brmshypothesis` object**,
  takes `class = "b"` by default and brms's argument order, supports
  `robust`, `seed` and `scope = "ranef"` or `"coef"`, and reports
  `Evid.Ratio` and `Post.Prob` in brms's columns. `evid_ratio`,
  `post_prob`, `z`, `p` and `attr(h, "draws")` are gone; the draws are
  `h$samples`. The string is read with brms's renaming, so `x:fe` is the
  interaction coefficient and not R's `:` operator, and a stored column
  such as `r_g[1,Intercept]` can be named with `class = NULL`, as in
  brms.

* **BREAKING: `bayes_R2()` on a multivariate model returns one row per
  response, `R2ya` and `R2y2`, as brms does.** It returned the first
  response alone, named `R2y_a`. `resp` takes brms's spelling of a
  response (`ya`).

* **BREAKING: brms's positional slots on six more methods.**
  `bayes_R2(ds, NULL, TRUE, TRUE)` now asks for brms's `robust` summary;
  the fourth slot used to be `probs` and returned a `Q100` column.
  `posterior_summary()`'s second slot is `pars`, `hypothesis()`'s third
  is `class`, `pairs()`'s second is `pars`, `conditional_effects()`'s
  third is `conditions` and `pp_check()`'s fourth is `prefix`. The
  formals audit of `dev/brmsmatch-beyond.R` now finds 0 of 68 methods
  diverging, from 6.

* **BREAKING: `summary()` has brms's columns and slots**: `Estimate`,
  `Est.Error`, `l-95% CI`, `u-95% CI`, `Rhat`, `Bulk_ESS`, `Tail_ESS`,
  with `prob`, `robust` and `mc_se`.

* **BREAKING: `conditional_effects()` defaults to brms's
  `robust = TRUE`**: `estimate__` is the median of the drawn curves and
  `se__` their MAD. `robust = FALSE` gives the mean and SD it used to.

* `pp_check()` takes brms's `prefix`, `group`, `x`, `newdata`, `resp`
  and `draw_ids`, and brms's default draw count with brms's message.
  `log_lik()` takes brms's `pointwise`, `combine`, `add_point_estimate`
  and `cores`, and refuses an argument it does not have, which it used
  to accept and ignore.

# frmtmb.sample 0.6.0

* **BREAKING, and a published number moves. `rhat()` and
  `neff_ratio()` now answer brms's question.** They called
  `bayesplot::rhat()` and `bayesplot::neff_ratio()` on the `stanfit`,
  which is `rstan::summary()`'s classic split-R-hat and rstan's
  `n_eff`. brms computes neither. `rhat.brmsfit` is
  `posterior::summarise_draws(rhat = posterior::rhat)`, the
  rank-normalized split-R-hat, and `neff_ratio.brmsfit` is
  `min(ess_bulk, ess_tail) / ndraws`. Both now do the same, computed on
  the draws array rather than on the `stanfit`.

  Measured on 4 chains of 500 (`dev/brmsmatch-findings.md`): the old
  `rhat()` differed from brms's by 0.00563108 at most, against a whole
  signal of `max|rhat - 1| = 0.0100318`, so by 1.15 times the entire
  excess over 1 that the diagnostic reports. `neff_ratio()` differed by
  0.43165, which is 77 times larger and is the number a user reads to
  decide whether to sample longer.

  These two no longer agree with `ds$stanfit`, which carries a
  different definition of the same idea.
  `rstan::summary(ds$stanfit)$summary` still has the sampler's own
  numbers for anyone who wants them.

* **BREAKING. `summary()`'s convergence columns move with them.** It
  reported rstan's `n_eff` and classic split-R-hat; it now reports
  `Rhat`, `Bulk_ESS` and `Tail_ESS` from posterior, which is what
  `brms:::summary.brmsfit()` reports, under brms's own column names.

  **The `n_eff` column is gone.** `summary(ds)[, "Rhat"]` is now
  exactly `rhat(ds)` on those rows and `neff_ratio(ds)` is exactly
  `pmin(Bulk_ESS, Tail_ESS) / ndraws(ds)`, so the package no longer
  disagrees with itself. It did: on the fit measured in
  `dev/reviews/20260915-brmsmatch.md`, `sigma_Intercept` read
  0.99990396 in the `Rhat` column against 1.0037623 from `rhat(ds)`,
  BELOW 1 in one place and above it in the other, and `x` read 1600.37
  effective draws against 1028.60, a 55.6 percent overstatement of the
  number a user reads to decide whether to sample longer.

* **`rhat()` and `neff_ratio()` now report frmtmb's parameter names.**
  They reported Stan's, because they read the `stanfit`: 4 of 11
  entries on the measured fit (`beta[1]`, `beta[2]`, `betad`, `theta`)
  were not addressable by the names `variables(ds)` lists, so
  `rhat(ds)["x"]` was `NA`. It is now a number. This falls out of the
  change above rather than being separate work.

* Both take brms's `pars`, and it is brms's `pars` FOR THESE TWO, which
  is not the `pars` that `as.mcmc()`, `mcmc_plot()` and
  `posterior_interval()` take. `brms:::rhat.brmsfit()` passes
  `variable = pars` to `as_draws_array()` rather than through
  `extract_pars()`, so `NULL` is the default and means every variable,
  a string is an EXACT variable name, `regex = TRUE` makes it a regular
  expression, and a name that is not there is an error. The other three
  keep the `extract_pars` rule, where `NA` means every variable and a
  string is a regular expression. `?draws-diagnostics` documents both
  under *Two different `pars` rules, both brms's*.

* **BREAKING. Ten methods take brms's arguments in brms's POSITIONS.**
  A positional call ported from brms used to answer a different
  question, twice without saying so. Every method below now matches
  brms's method name for name, as far as this package's own arguments
  go, and `tests/testthat/test-draws-spellings.R` asserts that against
  the installed brms.

  - `as.mcmc(x, pars, fixed, combine_chains, inc_warmup, ...)`, was
    `(x, combine_chains, ...)`. `as.mcmc(ds, TRUE)` used to return the
    pooled draws; brms reads that slot as `pars` and refuses a `pars`
    that is neither `NA` nor character, and so does this now.
  - `posterior_interval(object, pars, variable, prob, regex, fixed,
    ...)`, was `(object, prob, variable, ...)`.
    `posterior_interval(ds, 0.9)` used to return a 90% interval; it is
    now the same refusal brms gives.
  - `log_lik(object, newdata, re_formula, resp, ndraws, draw_ids,
    ...)`, was `(object, ndraws, resp, ...)`.
  - `mcmc_plot(object, pars, type, variable, regex, fixed, ...)`, was
    `(object, type, variable, ...)`.
  - `posterior_epred(object, newdata, re_formula, re.form, resp, dpar,
    nlpar, ndraws, draw_ids, ...)`, was `(object, newdata, resp,
    re_formula, re.form, ndraws, ...)`.
  - `posterior_linpred(object, transform, newdata, re_formula, re.form,
    resp, dpar, nlpar, incl_thres, ndraws, draw_ids, ...)`.
  - `posterior_predict(object, newdata, re_formula, re.form, transform,
    resp, negative_rt, ndraws, draw_ids, ...)`.
  - `pp_mixture(x, newdata, re_formula, resp, ndraws, draw_ids, log,
    summary, robust, probs, ...)`, was `(x, summary, ndraws, ...)`.
  - `predictive_error(object, newdata, re_formula, re.form, method,
    resp, ndraws, draw_ids, ...)`, was `(object, resp, re_formula,
    re.form, ndraws, ...)`.
  - `psis(log_ratios, newdata, resp, model_name, ndraws, ...)`, was
    `(log_ratios, ndraws, resp, ...)`.

  **A call that named its arguments is unaffected. A call that passed
  them positionally must be re-read.**

* New along the way, because brms puts them in those positions:
  `draw_ids` on the predictive methods, on `log_lik()` and on
  `pp_mixture()` names the draws to use by row index; `dpar` on
  `posterior_epred()`; `log`, `robust` and `probs` on `pp_mixture()`;
  `method` on `predictive_error()`, which also takes `newdata` now and
  re-evaluates the response term on it.

* Eight brms arguments are carried in brms's position and refused, each
  with the reason and the replacement: `log_lik(newdata =)`,
  `log_lik(re_formula =)`, `psis(newdata =)`, `pp_mixture(newdata =)`,
  `pp_mixture(re_formula =)`, `as.mcmc(inc_warmup = TRUE)`,
  `posterior_linpred(incl_thres = TRUE)` and
  `posterior_predict(negative_rt = TRUE)`. `re_formula = NULL` is
  brms's own default and is refused nowhere: it is the no-op these
  methods already do.

* Two spellings of the random-effect switch survive on the five
  methods where brms itself ACCEPTS both, which is not the same as the
  four that declare both. `posterior_epred()`, `posterior_linpred()`,
  `posterior_predict()` and `predictive_error()` declare `re_formula`
  and lme4's `re.form` on `brmsfit`. `predictive_interval()` declares
  neither and honors both anyway, because its whole brms body is
  `posterior_predict(object, ...)`, so the alias reaches a formal one
  frame down. `tests/testthat/test-draws-spellings.R` derives "brms
  accepts it" from brms's own code rather than restating a list.

* **BREAKING, one method.** `pp_check()` no longer takes `re.form`.
  Note that this is the one place the package is knowingly narrower
  than brms: `pp_check.brmsfit()` forwards its dots to
  `posterior_predict()` through `do_call()`, so brms accepts the alias
  there too. `dev/argspell-brms-accepts.R` in the monorepo has the
  measurement, and the decision is open rather than settled.

* **BREAKING.** An argument that lands in a method's `...` is an error
  naming it, rather than being swallowed. Every method of this package
  that took a `...` it never read now refuses one; `dev/argspell-report.R`
  in the monorepo reports 0 still swallowing. See the frmtmb NEWS entry
  for why the refusal is the part that prevents the class.

# frmtmb.sample 0.5.0

* **frmtmb.sample no longer breaks brms, rstantools, loo,
  bridgesampling, bayesplot, posterior, coda or gratia.** This package
  defined its own generic for 28 names those packages own, such as
  `log_lik()`, `posterior_epred()`, `psis()`, `rhat()` and
  `as.mcmc()`. `UseMethod()` reads the method table of the namespace
  where the generic it reached was defined, so after
  `library(brms); library(frmtmb.sample)`, `log_lik(fit)` on a
  `brmsfit` found no method. Measured on the previous release: all 28
  lost after brms, all 28 with brms only loaded, all 28 inside a
  package that imports frmtmb.sample, and every gratia method for
  `posterior_samples()` after gratia. It is now 0 in each of those
  orders. The package now calls `frmtmb::frm_install_generics()` from
  its `.onLoad()`, so each name resolves to the owner's generic while
  the owner is loaded and to this package's own while it is not.

* Two methods reached a foreign `.default` in silence and now reach
  their own. With posterior loaded after this package, `rhat()` on
  draws went to posterior's `rhat.default`; with gratia loaded after
  it, `posterior_samples()` went to gratia's. Both are now also
  registered on those owners' generics.

* **BREAKING.** Three generics take their owner's first argument, so a
  call that NAMED the old one fails:

  - `bridge_sampler(samples, ...)`, was `(x, ...)`, as bridgesampling.
  - `bayes_factor(x1, x2, log = FALSE, ...)`, was `(x, ...)`, as
    bridgesampling.
  - `rhat(x, ...)`, was `(object, ...)`. posterior and bayesplot both
    define `rhat` with different first arguments; this follows
    posterior, which is the one brms uses.

  Two only GAIN arguments and break nothing: `post_prob()` gains
  `prior_prob` and `model_names` from bridgesampling, and
  `posterior_samples()` gains `pars` from brms. `posterior_linpred()`'s
  generic gains `transform`, which its method already had.

  **If you write an S3 method on any of these names, match the
  owner's formals.**

* While brms is loaded, `parnames()` and `posterior_samples()` on draws
  give brms's own deprecation warning before this package's refusal,
  because the generic that runs is brms's.

* `gratia (>= 0.9.0)` joins Suggests, the first gratia with
  `posterior_samples()`, and `posterior` gains the floor `(>= 1.0.0)`,
  whose `rhat()` this package now registers a method on. As in frmtmb,
  a partial install of either on the library path that lacks the
  generic, loaded before this package, stops it loading.

* The `frmtmb` floor moves to the release that exports
  `frm_install_generics()`.

# frmtmb.sample 0.4.2

* **Requires frmtmb 0.56.0, and the requirement is a hard one.**
  `hypothesis()` on a `frmtmb_draws` object now arms the note about a
  covariate that shadows a distributional parameter itself, through
  `frmtmb::hyp_shadow_arm()` and `hyp_shadow_disarm()`, which no earlier
  frmtmb exports. Before this release the draws method relied on core's
  generic to arm that note; once frmtmb 0.56.0 stopped owning the
  `hypothesis` generic, the note fired on a fitted model and never on
  draws, in every session and whether or not brms was loaded.

* The generics this package re-exports from frmtmb no longer break
  brms's methods when brms is loaded. With brms loaded but not
  attached, 26 of 26 re-exports lost their brms method in 0.4.1 and 0
  of 26 do now. That repair comes from frmtmb 0.56.0 and reaches this
  package for free.

* Three method signatures now carry their generics' full formals,
  because an S3 method must accept every argument its generic declares
  and the generics are now the owners' rather than frmtmb's.

* This package still DEFINES 28 generics of its own that break brms's
  methods in the same way frmtmb's did. That is not fixed here.

# frmtmb.sample 0.4.1

* **A broken `tmbstan` now skips the sampler tests instead of erroring
  them**, and `skip_sampler()` asks the package's own detector rather
  than only whether `tmbstan` is installed. 111 test blocks reach a
  sampler; 47 of them had been using a file-local copy of
  `skip_sampler()`, so searching for the helper's name looked like
  coverage. All 111 go through it now.

* **The `\donttest` examples carried the same weak guard.** Fifteen
  topics called a sampler behind `requireNamespace()`, which a broken
  build satisfies. An example cannot skip, so it must not run.

* **A check that does not depend on the marker string.** On an affected
  build `rstan::grad_log_prob()` disagrees with `-fit$obj$gr()`, since
  the patched log-density overload is `require_not_st_var` and the
  missed one is `require_st_var`. With infinite bounds the two are the
  same number bitwise, 0 ulp over four model shapes and five points,
  so this catches an upstream that renames the placeholder, which the
  marker cannot.

* The check workflow pins the Stan trio to a dated snapshot and fails
  if `model.hpp` is still unpatched, reporting the state it measured
  rather than a boolean. Public RSPM rebuilt `tmbstan` the day after
  StanHeaders 2.39.1 landed, and four consecutive green runs of this
  package's checks sampled a standard normal instead of the model.

# frmtmb.sample 0.4.0

BREAKING. `frm_sample(control =)` is the SAMPLER's control list now,
brms's spelling and brms's meaning, passed to tmbstan unchanged. The
fit-time options move to `fit_control =`. `as_tmbstan()` already meant
`control` that way, so the two entry points agree rather than offering
a third convention, and brms code that tightens `adapt_delta` ports
across unchanged.

* That closes a silent discard. On the fitted-object route
  `control = list(adapt_delta = 0.99)` was IGNORED:
  `stan_args[[1]]$control` came back NULL where 0.97 was asked for,
  and the mean `accept_stat__` was 0.9266 against 0.9836 once the
  setting arrives. `fit_control`, `start`, `data2`, `na.action` and
  `REML` were silently ignored on that route too, and are refused by
  name. An abbreviation is refused as well: `contro =` fell into
  `...`, reached `rstan::sampling()` and partial matched onto ITS
  `control`, bypassing every check.

* `as_tmbstan()` no longer returns an empty `stanfit` without saying
  so. That was worse than a missing refusal: it bricked the session,
  and every later model failed with `empty_nested() must be true`.
  Both doors share one check, and a run that produced nothing is
  worded apart from a model that cannot be sampled at all, because
  those want different fixes.

* `frm_sample()` reads the compatibility registry before taping and
  refuses a fit the sampler cannot run, in about half a second with
  rstan never loaded, against the roughly three seconds the failure
  used to take. Where the formula position settles it the match is
  made on position; where it cannot, the match is made on names, and
  `?frm_sample` says that is the largest set justifiable that way
  rather than a sound one.

* **`frm_sample(control =)` is the SAMPLER's control list now, which is
  what brms means by the name.**
  `frm_sample(fit, control = list(adapt_delta = 0.99, max_treedepth =
  12))` ports straight from a `brm()` call and reaches rstan through
  tmbstan. The fit-time options `frm()` takes are `fit_control =`, and
  they still apply only on the formula route. This is a breaking
  rename, and the old shape is refused rather than reinterpreted: a
  `frmtmb_control()` list arriving in `control` is caught by the field
  names it carries and the refusal names `fit_control`. The two
  vocabularies share no name, which is what makes the two lists
  decidable, and a test asserts they stay disjoint. An option rstan
  does not have is also refused by name, because rstan answers an
  unknown one by declining to sample and returning an empty fit, which
  used to be reported as a solver failure. So is an ABBREVIATION of the
  name: `control` follows `...` in the formals, so only the exact
  spelling binds to it, and a shorter one would land in `...` and then
  partial-match rstan's own `control`, reaching the sampler unchecked.
  `vignette("brms-posterior")` carried a row saying `adapt_delta` had
  no route through this function; it now says the opposite.

* **`frm_sample()` refuses an assembly argument on a fitted object
  instead of accepting and discarding it.** `fit_control`, `start`,
  `data2`, `na.action` and `REML` build a model, and a fitted object is
  already built, so on that route they were read by nothing and said
  nothing: `frm_sample(fit, REML = TRUE)` sampled the model it had.
  `data` and `family` were already refused there; these five now are
  too, and the message points at `control` for the caller who meant
  the sampler. `fit_control` matters most, because it is the name the
  `control` refusal sends a confused user to.

* **Breaking, and a return-contract change: `as_tmbstan()` never hands
  back an empty `stanfit` now.** rstan answers some failures by
  declining to sample, printing its reason to the console and returning
  the shell of a `stanfit` with no error at all. `frm_sample()` has
  always refused that; `as_tmbstan()` passed it on, so a failed run
  looked like a successful one until an accessor came back empty. Both
  now raise the same message, naming the function that raised it. The
  wording distinguishes the two ways to get nothing out of this
  package, because they want different fixes: a model that cannot be
  sampled at all is refused by the compatibility pre-flight before any
  chain starts, and this one is a failure of the run in front of you.
  Nothing in the package or its vignettes relied on the empty object.

* **Fixed: a sampler control list handed to a fitted object was
  discarded without a word.** Up to and including 0.3.2,
  `frm_sample(fit, control = list(adapt_delta = 0.99))` on the fit
  route ran the chains with rstan's own defaults and said nothing:
  `control` was read only where the formula route assembles a model,
  and the fit route does not go there. A brms user porting a call that
  raised `adapt_delta` therefore got the untightened sampler. The
  argument now reaches rstan on both routes, and
  `ds$stanfit@stan_args[[1]]$control` records what was asked for.

* **A pre-flight against the compatibility registry**, on
  `frm_sample()` and on `as_tmbstan()`. Both read the rows that name
  `frm_sample` with the status `"refused"` before anything is taped,
  before the formula route calls `frm()`, and before the rstan and
  tmbstan namespaces are loaded, so a model the sampler cannot run
  costs the registry read and nothing else: 0.36 s to 0.42 s on the
  first call of a session against the 2.9 s failure it replaces, and
  that failure also left rstan unusable for the rest of the session.
  The refusal repeats the registered note and names the row.
  `frmtmb.ode` registers the first such row.

  What the check refuses on is deliberately narrow, because what
  identifies a feature in a model is WHERE its call sits in the
  formula. It decides three ways. A response FAMILY is matched by
  name, which is exact; it does not look inside `mixture()`. An
  ADDITION TERM is matched by POSITION, since it is written inside the
  bar on the left of the response formula and nowhere else, which
  separates `y | se(v) ~ x` from `y ~ a * se(x)` with a user's own
  `se()`. Everything else the vocabulary writes as a call is matched
  by NAME, outside addition-term position, and only when the name is
  not also a base R function; that route is justifiable from names but
  not sound, so a refusal registered on `mo()`, `ma()`, `cosy()` or
  `mm()` would also fire on a user's function of that name. Everything
  the check cannot decide warns and samples rather than refusing a
  model that may be correct. A registered refusal cannot be revoked by
  a later registration: display names are unique, core refuses a second
  registration of one, and a name the vocabulary writes without
  parentheses is declared not to be a call at all.

# frmtmb.sample 0.3.2

Requires frmtmb 0.55.0, for the hazard-container lint that now
runs in this package's own check.

# frmtmb.sample 0.3.1

* The sampling articles say why they are unevaluated, and the
  reason they used to give was false: four chains on the model in
  question take seconds, not the minutes claimed. They are
  unevaluated because the samplers are optional dependencies and
  because that page documents a Stan build which samples the wrong
  density silently. The two chunks that need no sampler now run,
  so the page shows output.

* `loo()` on a group-unit matrix says which unit it left out.
  `loo::loo.matrix()` never sees the attribute the log-likelihood
  matrix carries, so the printed estimate was indistinguishable
  from a per-observation one.

# frmtmb.sample 0.3.0

`conditional_effects()` on draws returns core's frame, grid and
display quantity; `hypothesis()` reports Savage-Dickey evidence ratios
validated against brms; `loo()` says which unit it leaves out. Requires
frmtmb 0.53.0 for the seam exports.

* `conditional_effects()` on a draws object now returns the frame
  core's fit method returns, which closes the divergence 0.2.0
  recorded. The columns and their order are core's `ce_frame()`: the
  varied predictor(s), then every other model variable at the value it
  is held at, then `cond__`, `cats__`, `effect1__` and `effect2__`,
  then the band. `cond__` is always present as a factor over every
  condition level, so brms's own `plot()` has something to facet on; a
  per-category ordinal display is keyed `"x:cats__"` rather than `"x"`,
  so brms's `plot()` reads the pair back out of the effects attribute;
  and `int_conditions =`, `categorical =` and `seed =` are formals here
  as they are on the fit method. A ported script that indexed the
  returned frame positionally, or that read `ce$x` on an ordinal fit,
  has to be read again.

* BUG FIX. `conditional_effects()` on draws drew the `mu` predictor
  where the fit method draws the EXPECTED RESPONSE, on every family
  whose mean is not the inverse link of `mu` and on every `trunc()`
  response. A zero-inflated fit was plotted `1 / (1 - zi)` times too
  high, a mixture showed one component's mean instead of the mixture's,
  and a response truncated below at zero reported a NEGATIVE mean and
  drew it. The column names, the grid and the keys were all correct
  throughout, which is why the first round's identity checks passed.
  The display quantity now comes from core's own `ce_pred_dpar()`, so
  the two surfaces cannot make different choices again. Measured after
  the fix, largest deviation from the fit curve over a five-point grid:
  1.4% zero-inflated, 1.0% truncated, 0.5% mixture, against 21.7% and a
  sign change before it. An explicitly named `dpar =` was always
  correct and is unchanged.

* `conditional_effects()` on draws accepts `allow_new_levels` and
  lme4's `allow.new.levels`, which the fit method has always accepted
  and this method warned about as unknown. It reads them with core's
  own `ce_dots()`, so the accepted set cannot drift again.

* `hypothesis()` refuses `abs(x) = 0` instead of reporting a Bayes
  factor twice too large. The affine check probed only non-negative
  points, where `abs()` is the identity; it now probes both signs. The
  prior of `|X|` for `X ~ N(0, 1)` is a half-normal with twice the
  density at zero, so the denominator was half what it should have
  been, silently.

* The `"evid_ratio_mcse"` attribute is measured over at least four
  blocks of draws. With two chains it was the standard deviation of two
  numbers over `sqrt(2)` - one degree of freedom - and understated a
  measured 0.46 gap between two runs of one model as 0.18. Chains are
  split in half when there are fewer than four, which keeps each block
  consecutive and so keeps the autocorrelation the estimator exists to
  respect.

* A point hypothesis on a natural-scale summary says so. `sigma = 0.9`
  used to refuse with "not a population-level coefficient", which is
  true of the name and false of the model: `sigma` IS a coefficient
  here, reported on a scale its prior is not written on. That is now
  the reason given, for every such quantity.

* BEHAVIOR CHANGE. `conditional_effects(re_formula = NULL)` on draws
  conditions on a NEW group, as it does on a fit, instead of silently
  taking the first observed level. The grid's grouping column reports
  `NA` and each POSTERIOR DRAW gets that group's effects drawn from its
  own covariance parameters, which is brms's construction; the fit
  method draws them per bootstrap replicate instead. The band is wider
  than the population band by the amount the new group's variance adds
  (1.92 to 3.93 in mean width on the model in `dev/sample-ce-findings.md`).
  `seed =` makes the draw reproducible.

* `conditional_effects()` on draws reports an argument it cannot use
  instead of discarding it in silence.

* `hypothesis()` on a draws object reports `evid_ratio` and
  `post_prob`, brms's `Evid.Ratio` and `Post.Prob`. A directional
  hypothesis gets the posterior odds of the claim, which needs no
  prior. A POINT hypothesis gets the Savage-Dickey density ratio: the
  posterior density of the tested quantity at the point, from a kernel
  density with `stats::density()`'s default `nrd0` bandwidth over 4096
  grid points read by spline, over the PRIOR density there, evaluated
  from the specification the model was sampled under rather than
  estimated from prior draws. brms estimates both halves by kernel
  density from `sample_prior = "yes"` draws; the two agree to the
  accuracy of that estimate, and the numbers are in
  `dev/sample-ce-findings.md`. The between-chain Monte Carlo error of
  each ratio rides on the returned object as the `"evid_ratio_mcse"`
  attribute, because a density ratio read without one is easy to
  over-believe.

* A point hypothesis whose parameter has no proper prior REFUSES BY
  NAME, with a warning that says which parameter and why, and leaves
  `evid_ratio` `NA`; the rest of the row is unaffected and the other
  hypotheses of the same call still report. `frm_sample()` leaves class
  `"b"` flat by default exactly as brms does, so this is the commonest
  answer until a prior is written. A hypothesis on a variance
  component, on a nonlinear function of the coefficients, or on the
  intercept (whose prior is about the intercept at the predictor means,
  not about the coefficient) refuses for its own stated reason.

* `log_lik()`, `loo()` and `waic()` accept a structured family that
  declares how its likelihood factorizes
  (`frmtmb::frmtmb_structure(loglik_group = )` or `(loglik_row = )`),
  instead of refusing every family whose likelihood is not rowwise.
  The columns are the pieces the family declares, at the coarsest
  granularity it gives, because a family that groups its likelihood is
  saying its rows are not independently droppable. The matrix carries
  `attr(x, "unit")` naming what a column is when it is not an
  observation. A family that supplies one number for the whole
  response is still refused, in the words it already used.

# frmtmb.sample 0.2.0

* BEHAVIOR CHANGE, following core. `ranef()` on a `frmtmb_draws`
  object is keyed by the GROUPING FACTOR rather than by the block, so
  `ranef(draws)[["1 | g"]]` returns `NULL` where it used to return an
  array and `ranef(draws)$g` returns it instead. The method delegates
  to core's `ranef()` once per draw, so it follows core's re-key
  exactly; the block label rides along in each array's `"term"`
  attribute, as it does in core. A model with two blocks on ONE factor
  gives two entries under one name, and the method now assembles them
  BY POSITION: keyed by name, `out[[tn]] <- st` wrote the same name
  twice, so the SECOND block was dropped entirely and the list came
  back with one entry instead of two.

* `conditional_effects()` on a draws object inherits four grid changes
  from core, because it builds its grids with core's
  `ce_grids_build()`: the two-variable grid now varies the FIRST effect
  slowest (brms's row order), a numeric moderator is held at the exact
  `mean +/- sd` rather than at `signif(mean +/- sd, 3)`, a `mo()`
  predictor gets one grid point per level instead of a 100-point
  continuous grid, and a `trials()` variable is held at 1 with a
  message that is new here. What it does NOT yet inherit is core's
  frame shape: the returned columns, the always-present `cond__`, the
  `"x:cats__"` key for a per-category display, `int_conditions =`,
  `categorical =` and the new-group meaning of `re_formula = NULL` are
  core-only for now, so the two surfaces disagree in both directions
  until the draws method is migrated.

* `frm_sample()` refuses parallel chains on Windows when the model's
  family comes from a namespace built by `pkgload::load_all()`, naming the
  package and the two remedies. A worker process cannot load such a
  namespace, and R silently substitutes the global environment for it,
  so the family failed at its first internal call and rstan reported
  only that the chains returned no draws.

* The defaults this package registers with `frmtmb::get_prior()` now
  answer `route = "sample"` only. Loading this package used to change
  the table `get_prior()` returned for every caller, an author asking
  about `frm()` included; it no longer does. Ask for
  `get_prior(..., route = "sample")` to see what `frm_sample()`
  applies, and with this package unloaded that route is refused rather
  than reported as flat. `frm_sample()` itself is unchanged: the same
  defaults under the same call-over-fit-over-defaults precedence.

* The frmtmb floor rises to 0.51.0. The registration in `.onLoad()`
  now passes `frmtmb_register_compat(expects = )`, which frmtmb gained
  in 0.51.0, and an older frmtmb stops the load outright with
  `unused argument (expects = c("hmm", "lca"))`. The old floor of
  0.46.0 permitted an installation that could not load, so this is a
  hard requirement rather than a preference.

* `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb
  for that package, so the tests here that sample one of them qualify
  the call and skip when it is absent. The two compatibility rules this
  package registers for those pairs are declared in
  `frmtmb_register_compat(expects = )`, which says they name a feature
  another package owns. Before frmtmb checked, those rules were dropped
  without a word in every session that did not also load
  `frmtmb.latent`; now `frm_compat()` reports them as unresolved until
  it is loaded.

# frmtmb.sample 0.1.0

First release, extracted from frmtmb 0.46.0.

* NUTS sampling for frmtmb models through tmbstan: `frm_sample()` on
  a fit or a formula, with brms 2.23's default priors on both routes,
  priors stacking most-explicit-first, and non-centered sampling of
  qualifying random-effect blocks. Draws carry frmtmb parameter
  names.
* The full posterior method surface on core's generics and on the
  brms, loo, posterior, bayesplot, rstantools, coda and
  bridgesampling generics: `log_lik()`, `loo()`, `waic()`,
  `posterior_epred()`, `posterior_predict()`, `hypothesis()`,
  `conditional_effects()` on draws, `bayes_R2()`, `mcmc_plot()`, and
  the rest.
* `check_laplace()` measures the Laplace and Wald approximations
  against NUTS on the same objective, centered and with no added
  priors.
* At load, the package registers its sampling default priors with
  core's `get_prior()` and its feature rows with the compatibility
  matrix, so both stay truthful whether or not this package is
  installed.
* Refuses a tmbstan built against StanHeaders 2.39 or later, whose
  code generator leaves one log-density overload unpatched so every
  chain silently samples a standard normal instead of the model.
* Parallel chains work on every platform, Windows included: tmbstan
  rebuilds the tape on each socket worker from the serialized
  objective closure, and frmtmb's generated closures are
  self-contained.
