# brms argument-name diff

Every frmtmb export that shares a name with a brms export, with its
formals diffed against brms's. Written for the release that finishes
the naming rule v0.40.1 started.

## The rule

A brms-NAMED function speaks brms's argument names. frmtmb's own fit
surface (`predict()`, `simulate()`, `frm_bootstrap()`,
`dharma_residuals()`) keeps lme4's, because that is the heritage the
name comes from and a reader should be able to tell which library a
call was written against.

Where a brms-named function has ALREADY shipped under the lme4
spelling, retiring it would break working code, so both spellings stay
live and resolve to one internal setting. brms does the same:
`posterior_epred.brmsfit()` carries `re_formula` and `re.form` side by
side. The one difference is that frmtmb REFUSES both at once, where
brms silently prefers `re.form`; two names for one setting supplied
together is a question about intent, not a preference to guess at.

## How this was produced

`intersect(getNamespaceExports("frmtmb"), getNamespaceExports("brms"))`
against brms 2.23.0 and frmtmb 0.40.1, then `formals()` on each pair.
For a generic, the S3 method is read rather than the generic (both are
`(x, ...)` otherwise): the frmtmb column prefers `.frmtmb_fit` and
falls back to `.frmtmb_draws`, the brms column reads `.brmsfit`. The
"signature read" column names exactly which function each row is a
claim about. `...` is excluded from all three set columns. Scripts:
`api-diff.R`, `api-probe3.R`, `api-table.R` in the session scratchpad.

Counts: 134 frmtmb exports, 306 brms exports, 92 shared names.

## The table

| function | frmtmb signature read | shared | brms-only (gap) | frmtmb-only (extra) |
| --- | --- | --- | --- | --- |
| `acat()` | `acat` | `link` | `link_disc`, `threshold` | - |
| `as.mcmc()` | `as.mcmc.frmtmb_draws` | `x`, `combine_chains` | `pars`, `fixed`, `inc_warmup` | - |
| `as_draws()` | `as_draws.frmtmb_draws` | `x` | `variable`, `regex`, `inc_warmup` | - |
| `as_draws_array()` | `as_draws_array.frmtmb_draws` | `x` | `variable`, `regex`, `inc_warmup` | - |
| `as_draws_df()` | `as_draws_df.frmtmb_draws` | `x` | `variable`, `regex`, `inc_warmup` | - |
| `as_draws_list()` | `as_draws_list.frmtmb_draws` | `x` | `variable`, `regex`, `inc_warmup` | - |
| `as_draws_matrix()` | `as_draws_matrix.frmtmb_draws` | `x` | `variable`, `regex`, `inc_warmup` | - |
| `as_draws_rvars()` | `as_draws_rvars.frmtmb_draws` | `x` | `variable`, `regex`, `inc_warmup` | - |
| `asym_laplace()` | `asym_laplace` | `link` | `link_sigma`, `link_quantile` | - |
| `bayes_factor()` | `bayes_factor.frmtmb_draws` | - | `x1`, `x2`, `log` | `x` |
| `bayes_R2()` | `bayes_R2.frmtmb_draws` | `object`, `resp`, `summary`, `probs` | `robust` | `ndraws` |
| `bernoulli()` | `bernoulli` | `link` | - | - |
| `Beta()` | `Beta` | `link` | `link_phi` | - |
| `beta_binomial()` | `beta_binomial` | `link` | `link_phi` | - |
| `bf()` | `bf` | `formula`, `family`, `nl` | `flist`, `autocor`, `loop`, `center`, `cmc`, `sparse`, `decomp` | - |
| `bridge_sampler()` | `bridge_sampler.frmtmb_draws` | - | `samples`, `recompile` | `x` |
| `categorical()` | `categorical` | `link` | `refcat` | `levels`, `K` |
| `conditional_effects()` | `conditional_effects.frmtmb_fit` | `x`, `effects`, `resolution`, `prob`, `method`, `re_formula`, `conditions`, `surface` | `int_conditions`, `robust`, `spaghetti`, `categorical`, `ordinal`, `transform`, `select_points`, `too_far`, `probs` | `resp`, `dpar`, `band`, `ndraws`, `boot`, `profile_points`, `seed`, `data` |
| `cox()` | `cox` | `link` | - | `df`, `degree`, `intercept` |
| `cratio()` | `cratio` | `link` | `link_disc`, `threshold` | - |
| `cumulative()` | `cumulative` | `link` | `link_disc`, `threshold` | - |
| `custom_family()` | `custom_family` | `dpars`, `links`, `type` | `name`, `lb`, `ub`, `vars`, `loop`, `specials`, `threshold`, `log_lik`, `posterior_predict`, `posterior_epred`, `predict`, `fitted`, `env` | `family`, `lpdf`, `valid_y`, `init_dpars`, `post`, `sim`, `sim_ctx`, `sim_refusal`, `primary_dpars`, `lcdf`, `extra_pars`, `drop_intercept` |
| `exgaussian()` | `exgaussian` | `link` | `link_sigma`, `link_beta` | - |
| `exponential()` | `exponential` | `link` | - | - |
| `expose_functions()` | `expose_functions.frmtmb_draws` | `x` | `vectorize`, `env` | - |
| `fixef()` | `fixef.frmtmb_fit` | `object` | `summary`, `robust`, `probs`, `pars` | - |
| `geometric()` | `geometric` | `link` | - | - |
| `get_prior()` | `get_prior` | `formula` | - | `data`, `family`, `data2` |
| `hurdle_gamma()` | `hurdle_gamma` | `link` | `link_shape`, `link_hu` | - |
| `hurdle_lognormal()` | `hurdle_lognormal` | `link` | `link_sigma`, `link_hu` | - |
| `hurdle_poisson()` | `hurdle_poisson` | `link` | `link_hu` | - |
| `hypothesis()` | `hypothesis.frmtmb_fit` | `x`, `hypothesis`, `alpha`, `seed`, `class`, `group` | `scope`, `robust` | `method`, `nsim`, `vcov` |
| `kfold()` | `kfold.frmtmb_draws` | `x` | `K`, `Ksub`, `folds`, `group`, `joint`, `compare`, `resp`, `model_names`, `save_fits`, `recompile`, `future_args` | - |
| `lf()` | `lf` | - | `flist`, `dpar`, `resp`, `center`, `cmc`, `sparse`, `decomp` | - |
| `log_lik()` | `log_lik.frmtmb_draws` | `object`, `ndraws`, `resp` | `newdata`, `re_formula`, `draw_ids`, `pointwise`, `combine`, `add_point_estimate`, `cores` | - |
| `log_posterior()` | `log_posterior.frmtmb_draws` | `object` | - | - |
| `lognormal()` | `lognormal` | `link` | `link_sigma` | - |
| `loo()` | `loo.frmtmb_draws` | `x`, `resp` | `compare`, `pointwise`, `moment_match`, `reloo`, `k_threshold`, `save_psis`, `moment_match_args`, `reloo_args`, `model_names` | `ndraws` |
| `LOO()` | `LOO.frmtmb_draws` | `x` | `compare`, `resp`, `pointwise`, `moment_match`, `reloo`, `k_threshold`, `save_psis`, `moment_match_args`, `reloo_args`, `model_names` | - |
| `loo_compare()` | `loo_compare.frmtmb_draws` | `x`, `criterion`, `model_names` | - | - |
| `loo_moment_match()` | `loo_moment_match.frmtmb_draws` | `x` | `loo`, `k_threshold`, `newdata`, `resp`, `check`, `recompile` | - |
| `loo_subsample()` | `loo_subsample.frmtmb_draws` | `x` | `compare`, `resp`, `model_names` | - |
| `mcmc_plot()` | `mcmc_plot.frmtmb_draws` | `object`, `type`, `variable` | `pars`, `regex`, `fixed` | - |
| `mixture()` | `mixture` | - | `flist`, `nmix`, `order` | `groups` |
| `multinomial()` | `multinomial` | - | `link`, `refcat` | `K` |
| `mvbf()` | `mvbf` | `rescor` | `flist` | - |
| `nchains()` | `nchains.frmtmb_draws` | `x` | - | - |
| `ndraws()` | `ndraws.frmtmb_draws` | `x` | - | - |
| `neff_ratio()` | `neff_ratio.frmtmb_draws` | `object` | `pars` | - |
| `negbinomial()` | `negbinomial` | `link` | `link_shape` | - |
| `ngrps()` | `ngrps.frmtmb_fit` | `object` | - | - |
| `niterations()` | `niterations.frmtmb_draws` | `x` | - | - |
| `nlf()` | `nlf` | `formula`, `loop` | `flist`, `dpar`, `resp` | - |
| `nsamples()` | `nsamples.frmtmb_draws` | `object` | `subset`, `incl_warmup` | - |
| `nuts_params()` | `nuts_params.frmtmb_draws` | `object` | `pars` | - |
| `nvariables()` | `nvariables.frmtmb_draws` | `x` | - | - |
| `parnames()` | `parnames.frmtmb_draws` | `x` | - | - |
| `post_prob()` | `post_prob.frmtmb_draws` | `x` | `prior_prob`, `model_names` | - |
| `posterior_epred()` | `posterior_epred.frmtmb_draws` | `object`, `newdata`, `resp`, `re_formula`, `re.form`, `ndraws` | `dpar`, `nlpar`, `draw_ids`, `sort` | - |
| `posterior_interval()` | `posterior_interval.frmtmb_draws` | `object`, `prob`, `variable` | `pars` | - |
| `posterior_linpred()` | `posterior_linpred.frmtmb_draws` | `object`, `transform`, `newdata`, `resp`, `re_formula`, `re.form`, `dpar`, `ndraws` | `nlpar`, `incl_thres`, `draw_ids`, `sort` | - |
| `posterior_predict()` | `posterior_predict.frmtmb_draws` | `object`, `newdata`, `resp`, `re_formula`, `re.form`, `ndraws` | `transform`, `negative_rt`, `draw_ids`, `sort`, `ntrys`, `cores` | - |
| `posterior_samples()` | `posterior_samples.frmtmb_draws` | `x` | `pars`, `fixed`, `add_chain`, `subset`, `as.matrix`, `as.array` | - |
| `posterior_summary()` | `posterior_summary.frmtmb_draws` | `probs`, `robust`, `variable` | `x`, `pars` | `object` |
| `pp_check()` | `pp_check.frmtmb_fit` | `object`, `type`, `ndraws` | `prefix`, `group`, `x`, `newdata`, `resp`, `draw_ids`, `nsamples`, `subset` | `re_formula`, `re.form` |
| `pp_mixture()` | `pp_mixture.frmtmb_draws` | `x`, `summary`, `ndraws` | `newdata`, `re_formula`, `resp`, `draw_ids`, `log`, `robust`, `probs` | - |
| `predictive_error()` | `predictive_error.frmtmb_draws` | `object`, `resp`, `re_formula`, `re.form`, `ndraws` | `newdata`, `method`, `draw_ids`, `sort` | - |
| `predictive_interval()` | `predictive_interval.frmtmb_draws` | `object`, `prob` | - | `newdata`, `resp`, `re_formula`, `re.form`, `ndraws` |
| `prior_summary()` | `prior_summary.frmtmb_fit` | `object` | `all` | - |
| `psis()` | `psis.frmtmb_draws` | `log_ratios`, `resp` | `newdata`, `model_name` | `ndraws` |
| `ranef()` | `ranef.frmtmb_fit` | `object` | `summary`, `robust`, `probs`, `pars`, `groups` | `condVar` |
| `reloo()` | `reloo.frmtmb_draws` | `x` | `loo`, `k_threshold`, `newdata`, `resp`, `check`, `recompile`, `future_args` | - |
| `restructure()` | `restructure.frmtmb_draws` | `x` | - | - |
| `rhat()` | `rhat.frmtmb_draws` | - | `x`, `pars` | `object` |
| `set_prior()` | `set_prior` | `prior`, `class`, `coef`, `dpar`, `group`, `lb`, `ub` | `resp`, `nlpar`, `tag`, `check` | - |
| `set_rescor()` | `set_rescor` | `rescor` | - | `rescor_value` |
| `shifted_lognormal()` | `shifted_lognormal` | `link` | `link_sigma`, `link_ndt` | - |
| `skew_normal()` | `skew_normal` | `link` | `link_sigma`, `link_alpha` | - |
| `sratio()` | `sratio` | `link` | `link_disc`, `threshold` | - |
| `stancode()` | `stancode.frmtmb_draws` | `object` | `version`, `regenerate`, `threads`, `backend` | - |
| `standata()` | `standata.frmtmb_draws` | `object` | `newdata`, `re_formula`, `newdata2`, `new_objects`, `incl_autocor` | - |
| `student()` | `student` | `link` | `link_sigma`, `link_nu` | - |
| `VarCorr()` | `VarCorr.frmtmb_fit` | `x` | `sigma`, `summary`, `robust`, `probs` | - |
| `variables()` | `variables.frmtmb_fit` | `x` | - | - |
| `von_mises()` | `von_mises` | `link` | `link_kappa` | - |
| `waic()` | `waic.frmtmb_draws` | `x`, `resp` | `compare`, `pointwise`, `model_names` | `ndraws` |
| `WAIC()` | `WAIC.frmtmb_draws` | `x` | `compare`, `resp`, `pointwise`, `model_names` | - |
| `weibull()` | `weibull` | `link` | `link_shape` | - |
| `zero_inflated_beta()` | `zero_inflated_beta` | `link` | `link_phi`, `link_zi` | - |
| `zero_inflated_binomial()` | `zero_inflated_binomial` | `link` | `link_zi` | - |
| `zero_inflated_negbinomial()` | `zero_inflated_negbinomial` | `link` | `link_shape`, `link_zi` | - |
| `zero_inflated_poisson()` | `zero_inflated_poisson` | `link` | `link_zi` | - |

The table rows for the seven functions changed in this release already
show the new signatures.

## Triage

The buckets overlap by design, so these are counts of rows carrying at
least one item of each kind, not a partition. A family constructor is
in (b) for keeping `link` alone on purpose AND in (c) for the `link_*`
arguments that would close the gap.

Structural counts, read straight off the table:

- 13 of the 92 rows match brms exactly, with no gap and no extra:
  `bernoulli()`, `exponential()`, `geometric()`, `log_posterior()`,
  `loo_compare()`, `nchains()`, `ndraws()`, `ngrps()`,
  `niterations()`, `nvariables()`, `parnames()`, `restructure()`,
  `variables()`.
- 8 methods across 7 exported functions are changed in this release
  (bucket a). Nothing else in the table is touched.
- 75 rows have at least one brms argument frmtmb lacks. The 22 family
  constructors and `conditional_effects()`, `hypothesis()`, `loo()`,
  `waic()`, `psis()`, `log_lik()`, `bayes_R2()`, `ranef()`,
  `custom_family()`, `multinomial()`, `categorical()`, `mixture()`,
  `bf()`, `lf()`, `nlf()` and `mvbf()` have a reason in (b) for the
  parts that are deliberate; every gap is logged in (c).
- 20 rows have a frmtmb-only argument. 15 are extensions with a reason
  in (b). The other 5 are the only pure NAME divergences in the whole
  diff: `posterior_summary(object)`, `rhat(object)`,
  `bayes_factor(x)`, `bridge_sampler(x)` and
  `set_rescor(rescor_value)`. `set_rescor()` is fixed in (a),
  `rhat()` is correct as it stands (see b), and the other three are
  in (c).

### (a) Mechanical alias fixes, made in this release

Same setting, two spellings, both accepted, both together refused.
Seven methods across six exported functions take the random-effect
switch, plus one boolean.

| function | primary (brms) | alias kept | default when neither is given |
| --- | --- | --- | --- |
| `posterior_epred.frmtmb_draws` | `re_formula` | `re.form` | `NULL` |
| `posterior_linpred.frmtmb_draws` | `re_formula` | `re.form` | `NULL` |
| `posterior_predict.frmtmb_draws` | `re_formula` | `re.form` | `NULL` |
| `predictive_interval.frmtmb_draws` | `re_formula` | `re.form` | `NULL` |
| `predictive_error.frmtmb_draws` | `re_formula` | `re.form` | `NULL` |
| `pp_check.frmtmb_fit` | `re_formula` | `re.form` | `NA` |
| `pp_check.frmtmb_draws` | `re_formula` | `re.form` | `NULL` (both new) |
| `set_rescor` | `rescor` | `rescor_value` | `TRUE` |

Notes.

- Every default is the one that was already there. `NULL` on the draws
  methods (a draw carries its own random effects, so conditioning on
  them is the point), `NA` on `pp_check()` for a fit (one point
  estimate, so new levels per replicate are what produces a spread at
  all). Only the spellings changed.
- `pp_check.frmtmb_draws` had NEITHER spelling: an `re.form` passed to
  it fell into the dots and reached the bayesplot `ppc_*` function as
  an unknown argument. It now takes both and forwards to
  `posterior_predict()`, which is what the fit method's `re.form` has
  always done through `simulate()`.
- `set_rescor(rescor_value = )` was never a considered choice: the
  argument shares an Rd page with `mvbf(rescor = )` and roxygen cannot
  carry two `@param rescor` entries on one page. brms spells it
  `rescor`; so does frmtmb now, with the old spelling kept as the
  alias. Every call site in the package, the tests and the vignettes is
  positional, so nothing had to change.
- The `probs` / `prob` hypothesis was checked and is EMPTY. Wherever
  frmtmb and brms both have the argument they agree on the spelling:
  `probs` on `bayes_R2()` and `posterior_summary()`, `prob` on
  `posterior_interval()`, `predictive_interval()` and
  `conditional_effects()`. The apparent divergences on `fixef()`,
  `ranef()`, `VarCorr()` and `pp_mixture()` are absences, not
  misspellings, and belong to bucket (c).

### (b) Deliberate divergences, kept

- **The fit surface keeps lme4's `re.form`.** `predict.frmtmb_fit()`,
  `simulate.frmtmb_fit()`, `frm_bootstrap()` and `dharma_residuals()`
  are frmtmb's own functions, not brms's, and their heritage is lme4
  and DHARMa. They do not gain `re_formula`, and a test asserts they
  do not.
- **`conditional_effects(method =, band =, boot =, profile_points =,
  seed =, ndraws =, data =, resp =, dpar =)`.** A frequentist band has
  choices a posterior band does not: `band` selects Wald, profile or
  bootstrap, `boot`/`profile_points`/`seed` drive them, `ndraws` sizes
  a resample. brms's `method` names a prediction function
  (`"posterior_epred"`); frmtmb's names the same idea for a fit.
- **`hypothesis(method =, nsim =, vcov =)`.** The frequentist test
  needs a covariance and, for a nonlinear contrast, a delta-method or
  simulation choice. brms's `robust` and `scope` are posterior
  summaries with no analog on a fit.
- **`cox(df =, degree =, intercept =)`.** brms's `cox()` has a fixed
  baseline-hazard basis; frmtmb's exposes the M-spline that
  approximates it.
- **`custom_family(family =, lpdf =, ...)` vs brms's `name =`.**
  Different objects: brms names a Stan function it will compile,
  frmtmb takes the R closure itself. `lb`, `ub`, `vars`, `loop`,
  `specials`, `env` are all Stan-code concerns.
- **`multinomial(K)`, `categorical(levels =, K =)`.** frmtmb needs the
  category count at family-construction time because the likelihood is
  built from it; brms reads it off the response and offers `refcat`
  instead.
- **`mixture(groups =)`.** A group-level latent class, which brms has
  no equivalent of. brms's `order` and `nmix` are Stan-side
  identification devices.
- **`get_prior(data =, family =, data2 =)`.** brms's `get_prior()` is
  a generic over `...`; frmtmb names the arguments it needs, which is
  the same call for every user who writes them by name.
- **`predictive_interval(newdata =, resp =, re_formula =, re.form =,
  ndraws =)`.** brms passes these through `...` to
  `posterior_predict()`. Naming them is compatible with every brms
  call and better documented.
- **`pp_check(re_formula =, re.form =)`.** brms's `pp_check.brmsfit()`
  has neither, because a `brmsfit` conditions on its draws. A
  frequentist check has to choose whether replicates get new random
  effects, so the argument must exist; it takes brms's spelling
  because `pp_check()` is a brms function.
- **`ranef(condVar =)`.** lme4/glmmTMB's conditional SDs of the modes,
  the fit-side analog of brms's `summary`/`probs` columns.
- **`bayes_R2(ndraws =)`, `loo(ndraws =)`, `waic(ndraws =)`,
  `psis(ndraws =)`, `log_lik(ndraws =)`.** Thinning matters here in a
  way it does not for brms: every draw re-runs the whole model
  evaluation rather than a vectorized Stan expression.
- **`rhat(object)` vs brms's `rhat(x)`.** frmtmb registers this method
  on `bayesplot::rhat`, whose first formal IS `object`. Following the
  generic it dispatches on is correct; the divergence is between brms
  and bayesplot, not between brms and frmtmb.
- **The 22 family constructors that take `link` only** (`acat()`,
  `asym_laplace()`, `Beta()`, `beta_binomial()`, `cratio()`,
  `cumulative()`, `exgaussian()`, `hurdle_gamma()`,
  `hurdle_lognormal()`, `hurdle_poisson()`, `lognormal()`,
  `negbinomial()`, `shifted_lognormal()`, `skew_normal()`, `sratio()`,
  `student()`, `von_mises()`, `weibull()`, `zero_inflated_beta()`,
  `zero_inflated_binomial()`, `zero_inflated_negbinomial()`,
  `zero_inflated_poisson()`). See (c) for what adding `link_*` would
  take; the divergence is deliberate for now, not accidental.

### (c) Gaps worth closing later, not implemented here

Effort notes are the implementation cost, not the doc cost.

- **`link_*` on 22 family constructors** (`link_sigma`, `link_shape`,
  `link_phi`, `link_zi`, `link_hu`, `link_disc`, `link_nu`,
  `link_alpha`, `link_beta`, `link_kappa`, `link_ndt`,
  `link_quantile`). Medium: the secondary dpar links are baked into
  each `fam_*()` definition, so this is a signature change plus a link
  lookup per dpar, times 22, plus the inverse-link plumbing in
  `predict()` and every simulator. Highest-value single item in this
  bucket, because a ported brms call with `link_sigma = "identity"`
  currently gets "unused argument".
- **`threshold =` on the four ordinal families** (`acat()`,
  `cratio()`, `cumulative()`, `sratio()`). Medium: `"flexible"` is
  what frmtmb does; `"equidistant"` and `"sum_to_zero"` are real
  reparameterizations of the cutpoints.
- **`summary = FALSE` / `robust` / `probs` on `fixef()`, `ranef()`,
  `VarCorr()`, `pp_mixture()`, `bayes_R2()`.** Small on the draws
  methods (the summarizing already happens in
  `draws_summarize_coef()`; `robust` is a median/MAD swap), meaningless
  on the fit methods, which have no draws to summarize. Do the draws
  half only.
- **`draw_ids =` and `sort =` across the predictive family**
  (`posterior_epred()`, `posterior_linpred()`, `posterior_predict()`,
  `predictive_error()`, `log_lik()`, `pp_mixture()`). Small:
  `draws_subsample()` already picks rows, so `draw_ids` is a
  pass-through of explicit indices; `sort` is a no-op here because
  frmtmb never reorders `newdata`.
- **`newdata =` on `predictive_error()`, `pp_mixture()`, `psis()`,
  `log_lik()`.** Small to medium: `posterior_predict()` already takes
  it, so the wiring exists, but each caller has to decide what the
  observed side of the comparison is on new rows.
- **`nlpar =` on `posterior_epred()` / `posterior_linpred()`.** Small
  given `dpar` already resolves; it is the nonlinear-parameter analog.
- **`pars` / `variable` / `regex` selection on `as_draws*()`,
  `mcmc_plot()`, `nuts_params()`, `neff_ratio()`, `posterior_samples()`.**
  Small: `draws_columns()` is the existing selector and needs a regex
  branch. `pars` is deprecated in brms; add `variable` and `regex`
  only.
- **`inc_warmup` on `as_draws*()` and `as.mcmc()`.** Blocked, not
  small: `frm_sample()` discards warmup before the draws matrix is
  built, so this is a change to what is stored, not to an accessor.
- **`compare =`, `model_names =`, `pointwise =` on `loo()` / `waic()`
  / `LOO()` / `WAIC()` / `loo_subsample()`.** Small: `loo_compare()`
  already has `model_names`, and `compare` is the "run
  `loo_compare()` on the extra arguments" convenience.
- **`posterior_summary(object)` should be `posterior_summary(x)`.**
  Small but BREAKING, so it needs a release that can carry it. This is
  a live defect, not just a style difference: the method is registered
  on `brms::posterior_summary`, whose generic dispatches on `x`, so
  `brms::posterior_summary(x = ds)` fails (the generic dispatches on `x`, our method signature starts at `object`, so the call errors with `argument "object" is missing, with no default`)
  while `posterior_summary(ds)` works. The fix is to rename the first
  formal on the generic, the `.default` method and the `.frmtmb_draws`
  method together. Found by the S3 first-formal audit in
  `api-probe3.R`, which reports exactly three mismatches package-wide.
- **`bayes_factor(x1, x2)` and `bridge_sampler(samples)` first
  formals.** The other two mismatches from that audit. Both methods
  are refusal stubs, so the only consequence is that
  `bridgesampling::bridge_sampler(samples = ds)` gives "unused
  argument" instead of the refusal that explains why marginal
  likelihoods are unavailable. Small, and worth doing the next time
  those stubs are touched.
- **`bf()` structural arguments** (`center`, `cmc`, `sparse`,
  `decomp`, `loop`, `autocor`, `flist`) and `lf()` / `nlf()`'s
  (`dpar`, `resp`, `flist`). Mixed: `autocor` has a frmtmb spelling
  already (terms inside the formula), `sparse` exists as a fit-level
  option, `center`/`cmc` are contrast-coding switches worth having,
  `decomp` is a QR reparameterization that would touch the objective.
- **`set_prior(resp =, nlpar =, tag =, check =)`.** Medium: `resp` is
  the multivariate selector and is the one that matters; the priors
  frame currently keys on class/coef/dpar/group only.
- **`hypothesis(scope =)`.** Medium: `"ranef"` and `"coef"` test
  group-level quantities, which needs the per-group covariance rather
  than the fixed-effect block.
- **`conditional_effects(int_conditions =, categorical =, ordinal =,
  transform =, select_points =, too_far =, spaghetti =)`.**
  `int_conditions` and `too_far` are small and useful (custom
  interaction levels; masking extrapolation on a surface).
  `categorical`/`ordinal` are partly covered by the automatic
  category handling. `spaghetti` needs per-draw curves and only makes
  sense with `band = "boot"`.
- **`prior_summary(all =)`, `VarCorr(sigma =)`, `nsamples(subset =)`,
  `posterior_interval(pars =)`, `expose_functions(vectorize =, env =)`,
  `stancode(version =, backend =, threads =)`,
  `standata(newdata2 =, incl_autocor =)`, `restructure()`,
  `posterior_samples()` extras.** Zero or near-zero value: the last
  five are on refusal stubs whose whole job is to explain that frmtmb
  has no Stan program, and adding arguments to a function that always
  stops buys nothing.
- **`kfold()`, `reloo()`, `loo_moment_match()` argument sets.** These
  are large features, not argument gaps; tracked in
  `dev/feature-gaps.md` rather than here.
