# Changelog

## frmtmb.sample 0.11.0

- The `pp_check(type = "error_binned")` refusal on draws of a
  multinomial fit calls the response a set of counts over categories,
  not a category, as the fit method now does.

- Two test files no longer leave `Rplots.pdf` in the tests directory:
  `test-draws-methods.R`, whose device guard covered one block and not
  the [`pairs()`](https://rdrr.io/r/graphics/pairs.html) call, and the
  generated `test-brms-suite-methods.R`, whose generator now opens a
  null device for the file.

- **[`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  on draws now draws a `cs()` term.** On an
  [`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html),
  [`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
  or
  [`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
  model with a category-specific effect it drew every row as if the term
  were absent, in sample and at `newdata`, with no warning.
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  on draws go through it and change with it.
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  and
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  were right and are unchanged. Needs the frmtmb that exports
  [`cs_offsets_add()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html).

- On draws, a `newdata` with a grouping level the fit did not see is
  refused by the draws methods themselves, with or without
  `allow_new_levels`, and the refusal says that
  `allow_new_levels = TRUE` is refused too. Without the flag, core’s
  message reached the caller and recommended “Use allow_new_levels =
  TRUE”, which the draws methods then refused. The same holds for
  `sample_new_levels` alone and for a `newdata` without the grouping
  column.

- `re_formula = NA` with `allow_new_levels = TRUE` and an unseen level
  now answers, as brms does, instead of being refused: the check asked
  about the level without the `re_formula` that drops it.

- The default priors of a categorical, multinomial, mixture or
  latent-class model carry the `dpar` of each location (`mub`, `mu1`),
  as brms writes them. They used to be written without it, which frmtmb
  now refuses, so this is what keeps those models sampling. The resolved
  defaults are the same parameters with the same densities. With no
  prior of your own, a categorical or mixture fit samples exactly as
  before. With a `dpar` prior of your own on a mixture, the default for
  that slot now steps aside instead of being written and then
  overridden, which changes the order the priors are applied in: the
  draws at a fixed seed differ (by up to 8.74 on one measured fit), and
  [`prior_summary()`](https://mc-stan.org/rstantools/reference/prior_summary.html)
  no longer lists the overridden defaults. The resp-, dpar- and
  nlpar-qualified defaults of the other families with several location
  dpars (`lba()`, `rdm()`, `hmm()`, `lca()`,
  [`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html))
  follow the same rule.

- The default-prior announcement names each slot with its `resp`, `dpar`
  and `nlpar`, so every line is a spelling
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
  takes. On a multivariate model it used to print `Intercept` with no
  response.

- [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)’s
  advice for a coefficient with no proper prior gives a full
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.html)
  call, with `resp`, `dpar`, `nlpar` and `coef` as the model needs them.
  It used to say `set_prior(class = "b")`, which a multivariate model
  refuses.

## frmtmb.sample 0.10.0

- **The default priors are brms’s on a multivariate model**: each
  response gets its own `Intercept`, `sigma` and `sd` defaults, read off
  that response, and `set_rescor(TRUE)` gets `lkj(1)` on the residual
  correlation. Until now a multivariate model was sampled with no
  default priors at all. This moves the draws of every multivariate
  model sampled without `prior = "flat"`.

- **The default intercept location takes off the mean
  [`offset()`](https://rdrr.io/r/stats/offset.html)**, as brms does:
  `y ~ 1 + offset(off)` with `y` near 2 and `off = 10` gets
  `student_t(3, -8, 2.5)`, where it got `student_t(3, 2, 2.5)`. A
  location is now written to 15 significant digits, as brms writes it.
  This moves the draws of a model with an offset in its location
  formula.

- **The `sd` default is written once per response, distributional
  parameter and nonlinear parameter**, to follow frmtmb’s brms scoping
  of class `"sd"`. With no `prior =`, every standard deviation gets the
  density it got before. With `prior = set_prior(..., class = "sd")`,
  that prior now replaces the location parameter’s default only, and a
  `phi` or `sigma` block keeps its own default, as in brms. A SMOOTH’s
  smoothing standard deviation is keyed the same way, so `sigma ~ s(z)`
  gets its own default and its own `sd` row. This moves the draws of a
  distributional model sampled with a class-wide `sd` prior.

- **BREAKING:
  [`predictive_error()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
  and [`residuals()`](https://rdrr.io/r/stats/residuals.html) on draws
  refuse an ordinal, categorical or multinomial model**, as brms refuses
  them.

- [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  on draws refuses `type = "error_binned"` for a category response in
  the same words the fit method uses, and ignores `resp` on a model with
  one response, which is what brms does with it. A multivariate model
  still selects with `resp` and still refuses a name none of its
  responses has.

- **A one-sided `re_formula` on draws keeps the terms it names**, in
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html). Each draw is
  evaluated through
  [`frmtmb::frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.html),
  which kept every term for any formula; the fix is in frmtmb’s
  development version, so it needs that frmtmb. On a
  `(1 + x | g) + (1 | h)` fit, `re_formula = ~ (1 | g)` used to be off
  the draw’s own intercept-only construction by up to 2.98; it now
  equals it. A term the fit does not have is an error that names it. The
  draws already carry each draw’s sampled group effects, so the
  uncertainty at a known level was in them before and nothing is added.

## frmtmb.sample 0.9.0

- **[`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  is frmtmb’s generic now, re-exported here**, rather than a second
  generic of the same name defined in this package. frmtmb owns the name
  because it has a method for it: `log_lik(fit)` on a maximum-likelihood
  fit refuses and names the route to draws, as
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`waic()`](https://mc-stan.org/loo/reference/waic.html) do, where it
  used to be “could not find function”. `log_lik(draws)` is unchanged,
  and the `frmtmb_draws` method is still registered on rstantools’
  generic as well. The method page is
  [`?log_lik.frmtmb_draws`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-log_lik.md);
  [`?frmtmb::log_lik`](https://aforren1.github.io/frmtmb/reference/log_lik.html)
  documents the generic.

- [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  and
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  on draws refuse every name in `...` they do not read, where they used
  to ignore it, as the other methods here already do. The new-level
  arguments pass through to the refusal described below.

- Requires frmtmb 0.61.0, for
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.html),
  [`fixef_by_dpar()`](https://aforren1.github.io/frmtmb/reference/fixef_by_dpar.html),
  [`brms_fixef_rows()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html)
  and the brms shape helpers.

- **[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) on draws
  exist** (item 2.6f). They are brms’s three summarizing methods, each a
  summary of the draws method beside it:
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  and
  [`predictive_error()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md).
  Before this there was no method at all, so `fitted(ds)` reached
  [`stats::fitted.default()`](https://rdrr.io/r/stats/fitted.values.html),
  read `ds$fitted.values` and returned `NULL`, a wrong answer with
  nothing said. `summary = FALSE` gives the draws, `robust` and `probs`
  are brms’s, and [`predict()`](https://rdrr.io/r/stats/predict.html) on
  an ordinal or categorical response returns brms’s `P(Y = k)`
  proportions.

- **BREAKING:
  [`nsamples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md),
  [`posterior_samples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  and
  [`parnames()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  answer instead of refusing.** All three were refused as “the
  deprecated brms spelling”. brms keeps all three live with a
  deprecation warning, and the refusal stopped a ported script dead
  where brms’s warning does not. Each now does what brms’s does and
  warns as brms warns:

  - `posterior_samples(x)` is `as.data.frame(x)`, with brms’s `pars` (a
    regular expression unless `fixed = TRUE`), `add_chain`, `subset`,
    `as.matrix` and `as.array`;
  - `nsamples(x)` is `ndraws(x)`, with brms’s `subset`;
  - `parnames(x)` is `variables(x)`.

  `nsamples(incl_warmup = TRUE)` is still refused:
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  discards the warmup rather than storing it, so there is nothing to
  count.

- **[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) on an
  ordinal draws object reports the thresholds**, `Intercept[1]`,
  `Intercept[2]`, then the coefficients and any `cs()` terms, in brms’s
  order and on the model’s own scale, which is what the fit’s
  [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) reports.
  It reported the slopes alone, because the sampler stores the
  thresholds as `tau_raw`.

- [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) on draws
  refuse `allow_new_levels = TRUE` when `newdata` holds a grouping level
  the fit did not see, with a message that names the function called and
  says why. That case used to fail in the design builder with “Use
  allow_new_levels = TRUE”, which the caller had just done. Predicting
  an unseen level from draws is not implemented;
  `predict(fit, allow_new_levels = TRUE)` on the maximum-likelihood fit
  does it. Every other use answers as before:
  `allow_new_levels = FALSE`, brms’s default, and `TRUE` with no
  `newdata` or with levels the fit saw.

- **BREAKING:
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  on draws no longer returns an object of class `brmshypothesis`.** The
  class is `"frmtmb_hypothesis"` alone; the shape and every element are
  unchanged, and what stops working is `is(x, "brmshypothesis")` and
  `inherits(x, "brmshypothesis")` in a ported script. frmtmb owns
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for its own
  class, so neither changes.

- **`point_estimate` and `ndraws_point_estimate` are honored** by
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  and
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html).
  They were accepted and ignored, so
  `posterior_epred(ds, point_estimate = "median", ndraws_point_estimate = 2)`
  returned all 25 draws where brms returns

  2.  As in brms it is a PARAMETER-space operation: the draws are
      collapsed to their mean or median first and the method runs once
      at that one vector, repeated `ndraws_point_estimate` times.

- [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  and
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  build each draw’s prediction with
  [`frmtmb::frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.html)
  rather than [`predict()`](https://rdrr.io/r/stats/predict.html), which
  is brms’s predictive summary in frmtmb 0.61.0.

## frmtmb.sample 0.8.0

- Requires frmtmb 0.60.0, for
  [`frm_stop()`](https://aforren1.github.io/frmtmb/reference/frmtmb-conditions.html)
  and the other condition helpers.

- `tmbstan (>= 1.2.1)` in Suggests, and `check_tmbstan_build()` now
  names that version as the remedy. The CI step that pinned a dated
  snapshot of the Stan trio is removed.

- Requires the frmtmb release that exports
  [`frm_stop()`](https://aforren1.github.io/frmtmb/reference/frmtmb-conditions.html);
  the `frmtmb (>= 0.59.0)` floor must move to it.

- **BREAKING:** every error, warning and message that frmtmb.sample
  raises is classed. An error has the class
  `c("frmtmb_sample_error", "frmtmb_error", "error", "condition")`, and
  warnings and messages follow the same pattern, so
  `tryCatch(frmtmb_error = )` catches any refusal. The class vector no
  longer contains `simpleError`, `simpleWarning` or `simpleMessage`. See
  `?frmtmb::frmtmb-conditions`.

- **BREAKING:** a value that matches none of the choices of
  `loo_compare(criterion =)`, `hypothesis(scope =)`,
  `pp_check(prefix =)` and `predictive_error(method =)` is refused with
  a `frmtmb_sample_error` that names the argument, the value and the
  choices. The old text was `'arg' should be one of ...`. A partial
  value still matches.

- [`posterior_interval()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md),
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
  [`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  and the other accessors refuse an exact variable name that the draws
  do not carry as a `frmtmb_sample_error`. The posterior package used to
  raise it as a plain error. The text keeps posterior’s words,
  `missing in the draws object`.

## frmtmb.sample 0.7.0

- Requires frmtmb 0.59.0, for
  [`brms_par_labels()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  [`brms_coef_table()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  [`brms_stan_name()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  [`expand_b()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  [`varcorr_layout()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  the
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  helpers
  ([`hyp_eval_in()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  [`hyp_expr_vars()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html))
  and
  [`check_prior_slots()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html)
  on `?frmtmb::frmtmb-sampling-api`.

- **BREAKING: the draws carry brms’s names.** `b_Intercept`, `b_x`,
  `b_sigma_Intercept` for the coefficients and `r_g[1,Intercept]` for a
  group-level coefficient (`r_gs[lev.1,Intercept]` for a level with a
  space, `r_g:h[1_p,Intercept]` for an interaction group), in
  `ds$draws`,
  [`variables()`](https://mc-stan.org/posterior/reference/variables.html),
  every `as_draws_*()`,
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
  [`as.array()`](https://rdrr.io/r/base/array.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html),
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
  [`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html).
  They were `Intercept`, `x`, `sigma_Intercept` and `b[1]`. A covariance
  parameter brms does not sample keeps the name
  [`confint()`](https://rdrr.io/r/stats/confint.html) gives it on the
  fit, `theta_1`; a block with no brms counterpart (reduced rank,
  smooth, GP, CAR, SPDE) keeps `b[i]`. What stops working:
  `ds$draws[, "x"]`, and anything that selected `b[` columns.

- **BREAKING: a distributional parameter nobody wrote a formula for is
  stored on its natural scale, under brms’s name.** The column is
  `sigma` (or `shape`, `nu`, and `sigma_ya` for response `y_a` of a
  multivariate model) and holds sigma, as brms’s does. It used to be
  `b_sigma_Intercept` holding log sigma, which brms reports only when
  `sigma ~ 1` is written out, so `fixef(ds)` had a row, `coef(ds)$g` a
  slice, and `posterior_summary(ds)` and `as_draws_df(ds)` a column that
  brms does not have. They are gone. What stops working:
  `ds$draws[, "b_sigma_Intercept"]` on such a model; use `sigma`, or
  `log(ds$draws[, "sigma"])` for the old values.

  A mixture’s weights are stored as brms’s simplex:
  `theta1 ... theta<K-1>` hold the mixing probabilities, and brms’s
  `thetaK` is a derived column before `lp__`. They used to be the log
  ratio under the name `theta1_Intercept`. Every reader that hands a
  draw to the model maps all `K` back to the log ratios together.

- **Draws labels a level repeats are suffixed as brms suffixes them.**
  Levels `lvl 1` and `lvl.1` are both `r_gd[lvl.1,Intercept]` in brms,
  which names the later one `r_gd[lvl.1,Intercept]__1`. Without the
  suffix every accessor failed on “Duplicate variable names”.
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  on draws refuses a name that matches two columns instead of reading
  the first.

- **BREAKING: the accessors return brms’s objects.**
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
  [`as.array()`](https://rdrr.io/r/base/array.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
  brms’s unclassed draws objects and take brms’s `pars`, `variable`,
  `draw` and `subset`, with brms’s deprecation warnings; `as_draws_*()`
  take `variable`, `regex` and `inc_warmup`;
  [`as_draws()`](https://mc-stan.org/posterior/reference/draws.html) is
  brms’s
  [`as_draws_list()`](https://mc-stan.org/posterior/reference/draws_list.html).
  `posterior_summary(x, pars, variable, probs, robust)` and
  [`posterior_interval()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
  cover every variable by default, as brms’s do. Compared against brms’s
  own installed methods on the same draws: 88 of 90 calls
  [`identical()`](https://rdrr.io/r/base/identical.html), the other 2
  `rvars` objects that differ only by a cache environment
  (`dev/brmsnames-findings.md`).

- **BREAKING:
  [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) are brms’s.**
  `summary = FALSE` returns the raw draws, as in brms, where it used to
  return a summary; `robust`, `probs`, `pars` and `groups` are brms’s.
  [`coef()`](https://rdrr.io/r/stats/coef.html) broadcasts every
  population-level coefficient, as brms’s does, and
  [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) has brms’s
  `sd`/`cor`/`cov`/`residual__` structure, with one `residual__` row per
  response on a multivariate model.
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) on a
  model with no group-level effects is refused with brms’s message.
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) and
  [`coef()`](https://rdrr.io/r/stats/coef.html) on a reduced-rank block,
  whose draws are factor scores rather than coefficients, compute the
  coefficients per draw; on `frm_sample(laplace = TRUE)` draws, which
  carry no group-level draws, they are refused by name rather than
  returned as `NA`.

- **BREAKING:
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  returns brms’s `brmshypothesis` object**, takes `class = "b"` by
  default and brms’s argument order, supports `robust`, `seed` and
  `scope = "ranef"` or `"coef"`, and reports `Evid.Ratio` and
  `Post.Prob` in brms’s columns. `evid_ratio`, `post_prob`, `z`, `p` and
  `attr(h, "draws")` are gone; the draws are `h$samples`. The string is
  read with brms’s renaming, so `x:fe` is the interaction coefficient
  and not R’s `:` operator, and a stored column such as
  `r_g[1,Intercept]` can be named with `class = NULL`, as in brms.

- **BREAKING:
  [`bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html)
  on a multivariate model returns one row per response, `R2ya` and
  `R2y2`, as brms does.** It returned the first response alone, named
  `R2y_a`. `resp` takes brms’s spelling of a response (`ya`).

- **BREAKING: brms’s positional slots on six more methods.**
  `bayes_R2(ds, NULL, TRUE, TRUE)` now asks for brms’s `robust` summary;
  the fourth slot used to be `probs` and returned a `Q100` column.
  [`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html)’s
  second slot is `pars`,
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)’s
  third is `class`, [`pairs()`](https://rdrr.io/r/graphics/pairs.html)’s
  second is `pars`,
  [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)’s
  third is `conditions` and
  [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)’s
  fourth is `prefix`. The formals audit of `dev/brmsmatch-beyond.R` now
  finds 0 of 68 methods diverging, from 6.

- **BREAKING: [`summary()`](https://rdrr.io/r/base/summary.html) has
  brms’s columns and slots**: `Estimate`, `Est.Error`, `l-95% CI`,
  `u-95% CI`, `Rhat`, `Bulk_ESS`, `Tail_ESS`, with `prob`, `robust` and
  `mc_se`.

- **BREAKING:
  [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  defaults to brms’s `robust = TRUE`**: `estimate__` is the median of
  the drawn curves and `se__` their MAD. `robust = FALSE` gives the mean
  and SD it used to.

- [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  takes brms’s `prefix`, `group`, `x`, `newdata`, `resp` and `draw_ids`,
  and brms’s default draw count with brms’s message.
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  takes brms’s `pointwise`, `combine`, `add_point_estimate` and `cores`,
  and refuses an argument it does not have, which it used to accept and
  ignore.

## frmtmb.sample 0.6.0

- **BREAKING, and a published number moves.
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  and
  [`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  now answer brms’s question.** They called
  [`bayesplot::rhat()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html)
  and
  [`bayesplot::neff_ratio()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html)
  on the `stanfit`, which is `rstan::summary()`’s classic split-R-hat
  and rstan’s `n_eff`. brms computes neither. `rhat.brmsfit` is
  `posterior::summarise_draws(rhat = posterior::rhat)`, the
  rank-normalized split-R-hat, and `neff_ratio.brmsfit` is
  `min(ess_bulk, ess_tail) / ndraws`. Both now do the same, computed on
  the draws array rather than on the `stanfit`.

  Measured on 4 chains of 500 (`dev/brmsmatch-findings.md`): the old
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  differed from brms’s by 0.00563108 at most, against a whole signal of
  `max|rhat - 1| = 0.0100318`, so by 1.15 times the entire excess over 1
  that the diagnostic reports.
  [`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  differed by 0.43165, which is 77 times larger and is the number a user
  reads to decide whether to sample longer.

  These two no longer agree with `ds$stanfit`, which carries a different
  definition of the same idea. `rstan::summary(ds$stanfit)$summary`
  still has the sampler’s own numbers for anyone who wants them.

- **BREAKING. [`summary()`](https://rdrr.io/r/base/summary.html)’s
  convergence columns move with them.** It reported rstan’s `n_eff` and
  classic split-R-hat; it now reports `Rhat`, `Bulk_ESS` and `Tail_ESS`
  from posterior, which is what `brms:::summary.brmsfit()` reports,
  under brms’s own column names.

  **The `n_eff` column is gone.** `summary(ds)[, "Rhat"]` is now exactly
  `rhat(ds)` on those rows and `neff_ratio(ds)` is exactly
  `pmin(Bulk_ESS, Tail_ESS) / ndraws(ds)`, so the package no longer
  disagrees with itself. It did: on the fit measured in
  `dev/reviews/20260915-brmsmatch.md`, `sigma_Intercept` read 0.99990396
  in the `Rhat` column against 1.0037623 from `rhat(ds)`, BELOW 1 in one
  place and above it in the other, and `x` read 1600.37 effective draws
  against 1028.60, a 55.6 percent overstatement of the number a user
  reads to decide whether to sample longer.

- **[`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  and
  [`neff_ratio()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  now report frmtmb’s parameter names.** They reported Stan’s, because
  they read the `stanfit`: 4 of 11 entries on the measured fit
  (`beta[1]`, `beta[2]`, `betad`, `theta`) were not addressable by the
  names `variables(ds)` lists, so `rhat(ds)["x"]` was `NA`. It is now a
  number. This falls out of the change above rather than being separate
  work.

- Both take brms’s `pars`, and it is brms’s `pars` FOR THESE TWO, which
  is not the `pars` that
  [`as.mcmc()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-as_draws.md),
  [`mcmc_plot()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  and
  [`posterior_interval()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
  take. `brms:::rhat.brmsfit()` passes `variable = pars` to
  [`as_draws_array()`](https://mc-stan.org/posterior/reference/draws_array.html)
  rather than through `extract_pars()`, so `NULL` is the default and
  means every variable, a string is an EXACT variable name,
  `regex = TRUE` makes it a regular expression, and a name that is not
  there is an error. The other three keep the `extract_pars` rule, where
  `NA` means every variable and a string is a regular expression.
  `?draws-diagnostics` documents both under *Two different `pars` rules,
  both brms’s*.

- **BREAKING. Ten methods take brms’s arguments in brms’s POSITIONS.** A
  positional call ported from brms used to answer a different question,
  twice without saying so. Every method below now matches brms’s method
  name for name, as far as this package’s own arguments go, and
  `tests/testthat/test-draws-spellings.R` asserts that against the
  installed brms.

  - `as.mcmc(x, pars, fixed, combine_chains, inc_warmup, ...)`, was
    `(x, combine_chains, ...)`. `as.mcmc(ds, TRUE)` used to return the
    pooled draws; brms reads that slot as `pars` and refuses a `pars`
    that is neither `NA` nor character, and so does this now.
  - `posterior_interval(object, pars, variable, prob, regex, fixed, ...)`,
    was `(object, prob, variable, ...)`. `posterior_interval(ds, 0.9)`
    used to return a 90% interval; it is now the same refusal brms
    gives.
  - `log_lik(object, newdata, re_formula, resp, ndraws, draw_ids, ...)`,
    was `(object, ndraws, resp, ...)`.
  - `mcmc_plot(object, pars, type, variable, regex, fixed, ...)`, was
    `(object, type, variable, ...)`.
  - `posterior_epred(object, newdata, re_formula, re.form, resp, dpar, nlpar, ndraws, draw_ids, ...)`,
    was `(object, newdata, resp, re_formula, re.form, ndraws, ...)`.
  - `posterior_linpred(object, transform, newdata, re_formula, re.form, resp, dpar, nlpar, incl_thres, ndraws, draw_ids, ...)`.
  - `posterior_predict(object, newdata, re_formula, re.form, transform, resp, negative_rt, ndraws, draw_ids, ...)`.
  - `pp_mixture(x, newdata, re_formula, resp, ndraws, draw_ids, log, summary, robust, probs, ...)`,
    was `(x, summary, ndraws, ...)`.
  - `predictive_error(object, newdata, re_formula, re.form, method, resp, ndraws, draw_ids, ...)`,
    was `(object, resp, re_formula, re.form, ndraws, ...)`.
  - `psis(log_ratios, newdata, resp, model_name, ndraws, ...)`, was
    `(log_ratios, ndraws, resp, ...)`.

  **A call that named its arguments is unaffected. A call that passed
  them positionally must be re-read.**

- New along the way, because brms puts them in those positions:
  `draw_ids` on the predictive methods, on
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  and on
  [`pp_mixture()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/pp_mixture.md)
  names the draws to use by row index; `dpar` on
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md);
  `log`, `robust` and `probs` on
  [`pp_mixture()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/pp_mixture.md);
  `method` on
  [`predictive_error()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md),
  which also takes `newdata` now and re-evaluates the response term on
  it.

- Eight brms arguments are carried in brms’s position and refused, each
  with the reason and the replacement: `log_lik(newdata =)`,
  `log_lik(re_formula =)`, `psis(newdata =)`, `pp_mixture(newdata =)`,
  `pp_mixture(re_formula =)`, `as.mcmc(inc_warmup = TRUE)`,
  `posterior_linpred(incl_thres = TRUE)` and
  `posterior_predict(negative_rt = TRUE)`. `re_formula = NULL` is brms’s
  own default and is refused nowhere: it is the no-op these methods
  already do.

- Two spellings of the random-effect switch survive on the five methods
  where brms itself ACCEPTS both, which is not the same as the four that
  declare both.
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  and
  [`predictive_error()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
  declare `re_formula` and lme4’s `re.form` on `brmsfit`.
  [`predictive_interval()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
  declares neither and honors both anyway, because its whole brms body
  is `posterior_predict(object, ...)`, so the alias reaches a formal one
  frame down. `tests/testthat/test-draws-spellings.R` derives “brms
  accepts it” from brms’s own code rather than restating a list.

- **BREAKING, one method.**
  [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  no longer takes `re.form`. Note that this is the one place the package
  is knowingly narrower than brms: `pp_check.brmsfit()` forwards its
  dots to
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  through `do_call()`, so brms accepts the alias there too.
  `dev/argspell-brms-accepts.R` in the monorepo has the measurement, and
  the decision is open rather than settled.

- **BREAKING.** An argument that lands in a method’s `...` is an error
  naming it, rather than being swallowed. Every method of this package
  that took a `...` it never read now refuses one;
  `dev/argspell-report.R` in the monorepo reports 0 still swallowing.
  See the frmtmb NEWS entry for why the refusal is the part that
  prevents the class.

## frmtmb.sample 0.5.0

- **frmtmb.sample no longer breaks brms, rstantools, loo,
  bridgesampling, bayesplot, posterior, coda or gratia.** This package
  defined its own generic for 28 names those packages own, such as
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`psis()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-loo.md),
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  and
  [`as.mcmc()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-as_draws.md).
  [`UseMethod()`](https://rdrr.io/r/base/UseMethod.html) reads the
  method table of the namespace where the generic it reached was
  defined, so after
  [`library(brms); library(frmtmb.sample)`](https://github.com/paul-buerkner/brms),
  `log_lik(fit)` on a `brmsfit` found no method. Measured on the
  previous release: all 28 lost after brms, all 28 with brms only
  loaded, all 28 inside a package that imports frmtmb.sample, and every
  gratia method for
  [`posterior_samples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  after gratia. It is now 0 in each of those orders. The package now
  calls
  [`frmtmb::frm_install_generics()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html)
  from its `.onLoad()`, so each name resolves to the owner’s generic
  while the owner is loaded and to this package’s own while it is not.

- Two methods reached a foreign `.default` in silence and now reach
  their own. With posterior loaded after this package,
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  on draws went to posterior’s `rhat.default`; with gratia loaded after
  it,
  [`posterior_samples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  went to gratia’s. Both are now also registered on those owners’
  generics.

- **BREAKING.** Three generics take their owner’s first argument, so a
  call that NAMED the old one fails:

  - `bridge_sampler(samples, ...)`, was `(x, ...)`, as bridgesampling.
  - `bayes_factor(x1, x2, log = FALSE, ...)`, was `(x, ...)`, as
    bridgesampling.
  - `rhat(x, ...)`, was `(object, ...)`. posterior and bayesplot both
    define `rhat` with different first arguments; this follows
    posterior, which is the one brms uses.

  Two only GAIN arguments and break nothing:
  [`post_prob()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-loo-refusals.md)
  gains `prior_prob` and `model_names` from bridgesampling, and
  [`posterior_samples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  gains `pars` from brms.
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)’s
  generic gains `transform`, which its method already had.

  **If you write an S3 method on any of these names, match the owner’s
  formals.**

- While brms is loaded,
  [`parnames()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  and
  [`posterior_samples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md)
  on draws give brms’s own deprecation warning before this package’s
  refusal, because the generic that runs is brms’s.

- `gratia (>= 0.9.0)` joins Suggests, the first gratia with
  [`posterior_samples()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frmtmb-draws-deprecated.md),
  and `posterior` gains the floor `(>= 1.0.0)`, whose
  [`rhat()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md)
  this package now registers a method on. As in frmtmb, a partial
  install of either on the library path that lacks the generic, loaded
  before this package, stops it loading.

- The `frmtmb` floor moves to the release that exports
  [`frm_install_generics()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html).

## frmtmb.sample 0.4.2

- **Requires frmtmb 0.56.0, and the requirement is a hard one.**
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  on a `frmtmb_draws` object now arms the note about a covariate that
  shadows a distributional parameter itself, through
  `frmtmb::hyp_shadow_arm()` and `hyp_shadow_disarm()`, which no earlier
  frmtmb exports. Before this release the draws method relied on core’s
  generic to arm that note; once frmtmb 0.56.0 stopped owning the
  `hypothesis` generic, the note fired on a fitted model and never on
  draws, in every session and whether or not brms was loaded.

- The generics this package re-exports from frmtmb no longer break
  brms’s methods when brms is loaded. With brms loaded but not attached,
  26 of 26 re-exports lost their brms method in 0.4.1 and 0 of 26 do
  now. That repair comes from frmtmb 0.56.0 and reaches this package for
  free.

- Three method signatures now carry their generics’ full formals,
  because an S3 method must accept every argument its generic declares
  and the generics are now the owners’ rather than frmtmb’s.

- This package still DEFINES 28 generics of its own that break brms’s
  methods in the same way frmtmb’s did. That is not fixed here.

## frmtmb.sample 0.4.1

- **A broken `tmbstan` now skips the sampler tests instead of erroring
  them**, and `skip_sampler()` asks the package’s own detector rather
  than only whether `tmbstan` is installed. 111 test blocks reach a
  sampler; 47 of them had been using a file-local copy of
  `skip_sampler()`, so searching for the helper’s name looked like
  coverage. All 111 go through it now.

- **The `\donttest` examples carried the same weak guard.** Fifteen
  topics called a sampler behind
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html), which a
  broken build satisfies. An example cannot skip, so it must not run.

- **A check that does not depend on the marker string.** On an affected
  build
  [`rstan::grad_log_prob()`](https://mc-stan.org/rstan/reference/stanfit-method-logprob.html)
  disagrees with `-fit$obj$gr()`, since the patched log-density overload
  is `require_not_st_var` and the missed one is `require_st_var`. With
  infinite bounds the two are the same number bitwise, 0 ulp over four
  model shapes and five points, so this catches an upstream that renames
  the placeholder, which the marker cannot.

- The check workflow pins the Stan trio to a dated snapshot and fails if
  `model.hpp` is still unpatched, reporting the state it measured rather
  than a boolean. Public RSPM rebuilt `tmbstan` the day after
  StanHeaders 2.39.1 landed, and four consecutive green runs of this
  package’s checks sampled a standard normal instead of the model.

## frmtmb.sample 0.4.0

BREAKING. `frm_sample(control =)` is the SAMPLER’s control list now,
brms’s spelling and brms’s meaning, passed to tmbstan unchanged. The
fit-time options move to `fit_control =`.
[`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
already meant `control` that way, so the two entry points agree rather
than offering a third convention, and brms code that tightens
`adapt_delta` ports across unchanged.

- That closes a silent discard. On the fitted-object route
  `control = list(adapt_delta = 0.99)` was IGNORED:
  `stan_args[[1]]$control` came back NULL where 0.97 was asked for, and
  the mean `accept_stat__` was 0.9266 against 0.9836 once the setting
  arrives. `fit_control`, `start`, `data2`, `na.action` and `REML` were
  silently ignored on that route too, and are refused by name. An
  abbreviation is refused as well: `contro =` fell into `...`, reached
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)
  and partial matched onto ITS `control`, bypassing every check.

- [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
  no longer returns an empty `stanfit` without saying so. That was worse
  than a missing refusal: it bricked the session, and every later model
  failed with `empty_nested() must be true`. Both doors share one check,
  and a run that produced nothing is worded apart from a model that
  cannot be sampled at all, because those want different fixes.

- [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  reads the compatibility registry before taping and refuses a fit the
  sampler cannot run, in about half a second with rstan never loaded,
  against the roughly three seconds the failure used to take. Where the
  formula position settles it the match is made on position; where it
  cannot, the match is made on names, and
  [`?frm_sample`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  says that is the largest set justifiable that way rather than a sound
  one.

- **`frm_sample(control =)` is the SAMPLER’s control list now, which is
  what brms means by the name.**
  `frm_sample(fit, control = list(adapt_delta = 0.99, max_treedepth = 12))`
  ports straight from a `brm()` call and reaches rstan through tmbstan.
  The fit-time options
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) takes
  are `fit_control =`, and they still apply only on the formula route.
  This is a breaking rename, and the old shape is refused rather than
  reinterpreted: a
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.html)
  list arriving in `control` is caught by the field names it carries and
  the refusal names `fit_control`. The two vocabularies share no name,
  which is what makes the two lists decidable, and a test asserts they
  stay disjoint. An option rstan does not have is also refused by name,
  because rstan answers an unknown one by declining to sample and
  returning an empty fit, which used to be reported as a solver failure.
  So is an ABBREVIATION of the name: `control` follows `...` in the
  formals, so only the exact spelling binds to it, and a shorter one
  would land in `...` and then partial-match rstan’s own `control`,
  reaching the sampler unchecked.
  [`vignette("brms-posterior")`](https://aforren1.github.io/frmtmb/frmtmb.sample/articles/brms-posterior.md)
  carried a row saying `adapt_delta` had no route through this function;
  it now says the opposite.

- **[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  refuses an assembly argument on a fitted object instead of accepting
  and discarding it.** `fit_control`, `start`, `data2`, `na.action` and
  `REML` build a model, and a fitted object is already built, so on that
  route they were read by nothing and said nothing:
  `frm_sample(fit, REML = TRUE)` sampled the model it had. `data` and
  `family` were already refused there; these five now are too, and the
  message points at `control` for the caller who meant the sampler.
  `fit_control` matters most, because it is the name the `control`
  refusal sends a confused user to.

- **Breaking, and a return-contract change:
  [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
  never hands back an empty `stanfit` now.** rstan answers some failures
  by declining to sample, printing its reason to the console and
  returning the shell of a `stanfit` with no error at all.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  has always refused that;
  [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md)
  passed it on, so a failed run looked like a successful one until an
  accessor came back empty. Both now raise the same message, naming the
  function that raised it. The wording distinguishes the two ways to get
  nothing out of this package, because they want different fixes: a
  model that cannot be sampled at all is refused by the compatibility
  pre-flight before any chain starts, and this one is a failure of the
  run in front of you. Nothing in the package or its vignettes relied on
  the empty object.

- **Fixed: a sampler control list handed to a fitted object was
  discarded without a word.** Up to and including 0.3.2,
  `frm_sample(fit, control = list(adapt_delta = 0.99))` on the fit route
  ran the chains with rstan’s own defaults and said nothing: `control`
  was read only where the formula route assembles a model, and the fit
  route does not go there. A brms user porting a call that raised
  `adapt_delta` therefore got the untightened sampler. The argument now
  reaches rstan on both routes, and `ds$stanfit@stan_args[[1]]$control`
  records what was asked for.

- **A pre-flight against the compatibility registry**, on
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  and on
  [`as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.md).
  Both read the rows that name `frm_sample` with the status `"refused"`
  before anything is taped, before the formula route calls
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html), and
  before the rstan and tmbstan namespaces are loaded, so a model the
  sampler cannot run costs the registry read and nothing else: 0.36 s to
  0.42 s on the first call of a session against the 2.9 s failure it
  replaces, and that failure also left rstan unusable for the rest of
  the session. The refusal repeats the registered note and names the
  row. `frmtmb.ode` registers the first such row.

  What the check refuses on is deliberately narrow, because what
  identifies a feature in a model is WHERE its call sits in the formula.
  It decides three ways. A response FAMILY is matched by name, which is
  exact; it does not look inside
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html).
  An ADDITION TERM is matched by POSITION, since it is written inside
  the bar on the left of the response formula and nowhere else, which
  separates `y | se(v) ~ x` from `y ~ a * se(x)` with a user’s own
  `se()`. Everything else the vocabulary writes as a call is matched by
  NAME, outside addition-term position, and only when the name is not
  also a base R function; that route is justifiable from names but not
  sound, so a refusal registered on `mo()`, `ma()`, `cosy()` or `mm()`
  would also fire on a user’s function of that name. Everything the
  check cannot decide warns and samples rather than refusing a model
  that may be correct. A registered refusal cannot be revoked by a later
  registration: display names are unique, core refuses a second
  registration of one, and a name the vocabulary writes without
  parentheses is declared not to be a call at all.

## frmtmb.sample 0.3.2

Requires frmtmb 0.55.0, for the hazard-container lint that now runs in
this package’s own check.

## frmtmb.sample 0.3.1

- The sampling articles say why they are unevaluated, and the reason
  they used to give was false: four chains on the model in question take
  seconds, not the minutes claimed. They are unevaluated because the
  samplers are optional dependencies and because that page documents a
  Stan build which samples the wrong density silently. The two chunks
  that need no sampler now run, so the page shows output.

- [`loo()`](https://mc-stan.org/loo/reference/loo.html) on a group-unit
  matrix says which unit it left out.
  [`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html)
  never sees the attribute the log-likelihood matrix carries, so the
  printed estimate was indistinguishable from a per-observation one.

## frmtmb.sample 0.3.0

[`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
on draws returns core’s frame, grid and display quantity;
[`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
reports Savage-Dickey evidence ratios validated against brms;
[`loo()`](https://mc-stan.org/loo/reference/loo.html) says which unit it
leaves out. Requires frmtmb 0.53.0 for the seam exports.

- [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  on a draws object now returns the frame core’s fit method returns,
  which closes the divergence 0.2.0 recorded. The columns and their
  order are core’s
  [`ce_frame()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html):
  the varied predictor(s), then every other model variable at the value
  it is held at, then `cond__`, `cats__`, `effect1__` and `effect2__`,
  then the band. `cond__` is always present as a factor over every
  condition level, so brms’s own
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) has something
  to facet on; a per-category ordinal display is keyed `"x:cats__"`
  rather than `"x"`, so brms’s
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) reads the
  pair back out of the effects attribute; and `int_conditions =`,
  `categorical =` and `seed =` are formals here as they are on the fit
  method. A ported script that indexed the returned frame positionally,
  or that read `ce$x` on an ordinal fit, has to be read again.

- BUG FIX.
  [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  on draws drew the `mu` predictor where the fit method draws the
  EXPECTED RESPONSE, on every family whose mean is not the inverse link
  of `mu` and on every [`trunc()`](https://rdrr.io/r/base/Round.html)
  response. A zero-inflated fit was plotted `1 / (1 - zi)` times too
  high, a mixture showed one component’s mean instead of the mixture’s,
  and a response truncated below at zero reported a NEGATIVE mean and
  drew it. The column names, the grid and the keys were all correct
  throughout, which is why the first round’s identity checks passed. The
  display quantity now comes from core’s own
  [`ce_pred_dpar()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  so the two surfaces cannot make different choices again. Measured
  after the fix, largest deviation from the fit curve over a five-point
  grid: 1.4% zero-inflated, 1.0% truncated, 0.5% mixture, against 21.7%
  and a sign change before it. An explicitly named `dpar =` was always
  correct and is unchanged.

- [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  on draws accepts `allow_new_levels` and lme4’s `allow.new.levels`,
  which the fit method has always accepted and this method warned about
  as unknown. It reads them with core’s own
  [`ce_dots()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html),
  so the accepted set cannot drift again.

- [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  refuses `abs(x) = 0` instead of reporting a Bayes factor twice too
  large. The affine check probed only non-negative points, where
  [`abs()`](https://rdrr.io/r/base/MathFun.html) is the identity; it now
  probes both signs. The prior of `|X|` for `X ~ N(0, 1)` is a
  half-normal with twice the density at zero, so the denominator was
  half what it should have been, silently.

- The `"evid_ratio_mcse"` attribute is measured over at least four
  blocks of draws. With two chains it was the standard deviation of two
  numbers over `sqrt(2)` - one degree of freedom - and understated a
  measured 0.46 gap between two runs of one model as 0.18. Chains are
  split in half when there are fewer than four, which keeps each block
  consecutive and so keeps the autocorrelation the estimator exists to
  respect.

- A point hypothesis on a natural-scale summary says so. `sigma = 0.9`
  used to refuse with “not a population-level coefficient”, which is
  true of the name and false of the model: `sigma` IS a coefficient
  here, reported on a scale its prior is not written on. That is now the
  reason given, for every such quantity.

- BEHAVIOR CHANGE. `conditional_effects(re_formula = NULL)` on draws
  conditions on a NEW group, as it does on a fit, instead of silently
  taking the first observed level. The grid’s grouping column reports
  `NA` and each POSTERIOR DRAW gets that group’s effects drawn from its
  own covariance parameters, which is brms’s construction; the fit
  method draws them per bootstrap replicate instead. The band is wider
  than the population band by the amount the new group’s variance adds
  (1.92 to 3.93 in mean width on the model in
  `dev/sample-ce-findings.md`). `seed =` makes the draw reproducible.

- [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  on draws reports an argument it cannot use instead of discarding it in
  silence.

- [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  on a draws object reports `evid_ratio` and `post_prob`, brms’s
  `Evid.Ratio` and `Post.Prob`. A directional hypothesis gets the
  posterior odds of the claim, which needs no prior. A POINT hypothesis
  gets the Savage-Dickey density ratio: the posterior density of the
  tested quantity at the point, from a kernel density with
  [`stats::density()`](https://rdrr.io/r/stats/density.html)’s default
  `nrd0` bandwidth over 4096 grid points read by spline, over the PRIOR
  density there, evaluated from the specification the model was sampled
  under rather than estimated from prior draws. brms estimates both
  halves by kernel density from `sample_prior = "yes"` draws; the two
  agree to the accuracy of that estimate, and the numbers are in
  `dev/sample-ce-findings.md`. The between-chain Monte Carlo error of
  each ratio rides on the returned object as the `"evid_ratio_mcse"`
  attribute, because a density ratio read without one is easy to
  over-believe.

- A point hypothesis whose parameter has no proper prior REFUSES BY
  NAME, with a warning that says which parameter and why, and leaves
  `evid_ratio` `NA`; the rest of the row is unaffected and the other
  hypotheses of the same call still report.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  leaves class `"b"` flat by default exactly as brms does, so this is
  the commonest answer until a prior is written. A hypothesis on a
  variance component, on a nonlinear function of the coefficients, or on
  the intercept (whose prior is about the intercept at the predictor
  means, not about the coefficient) refuses for its own stated reason.

- [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) and
  [`waic()`](https://mc-stan.org/loo/reference/waic.html) accept a
  structured family that declares how its likelihood factorizes
  (`frmtmb::frmtmb_structure(loglik_group = )` or `(loglik_row = )`),
  instead of refusing every family whose likelihood is not rowwise. The
  columns are the pieces the family declares, at the coarsest
  granularity it gives, because a family that groups its likelihood is
  saying its rows are not independently droppable. The matrix carries
  `attr(x, "unit")` naming what a column is when it is not an
  observation. A family that supplies one number for the whole response
  is still refused, in the words it already used.

## frmtmb.sample 0.2.0

- BEHAVIOR CHANGE, following core.
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) on a
  `frmtmb_draws` object is keyed by the GROUPING FACTOR rather than by
  the block, so `ranef(draws)[["1 | g"]]` returns `NULL` where it used
  to return an array and `ranef(draws)$g` returns it instead. The method
  delegates to core’s
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) once per
  draw, so it follows core’s re-key exactly; the block label rides along
  in each array’s `"term"` attribute, as it does in core. A model with
  two blocks on ONE factor gives two entries under one name, and the
  method now assembles them BY POSITION: keyed by name,
  `out[[tn]] <- st` wrote the same name twice, so the SECOND block was
  dropped entirely and the list came back with one entry instead of two.

- [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  on a draws object inherits four grid changes from core, because it
  builds its grids with core’s
  [`ce_grids_build()`](https://aforren1.github.io/frmtmb/reference/frmtmb-sampling-api.html):
  the two-variable grid now varies the FIRST effect slowest (brms’s row
  order), a numeric moderator is held at the exact `mean +/- sd` rather
  than at `signif(mean +/- sd, 3)`, a `mo()` predictor gets one grid
  point per level instead of a 100-point continuous grid, and a
  `trials()` variable is held at 1 with a message that is new here. What
  it does NOT yet inherit is core’s frame shape: the returned columns,
  the always-present `cond__`, the `"x:cats__"` key for a per-category
  display, `int_conditions =`, `categorical =` and the new-group meaning
  of `re_formula = NULL` are core-only for now, so the two surfaces
  disagree in both directions until the draws method is migrated.

- [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  refuses parallel chains on Windows when the model’s family comes from
  a namespace built by
  [`pkgload::load_all()`](https://pkgload.r-lib.org/reference/load_all.html),
  naming the package and the two remedies. A worker process cannot load
  such a namespace, and R silently substitutes the global environment
  for it, so the family failed at its first internal call and rstan
  reported only that the chains returned no draws.

- The defaults this package registers with
  [`frmtmb::get_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.html)
  now answer `route = "sample"` only. Loading this package used to
  change the table
  [`get_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.html)
  returned for every caller, an author asking about
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
  included; it no longer does. Ask for
  `get_prior(..., route = "sample")` to see what
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  applies, and with this package unloaded that route is refused rather
  than reported as flat.
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  itself is unchanged: the same defaults under the same
  call-over-fit-over-defaults precedence.

- The frmtmb floor rises to 0.51.0. The registration in `.onLoad()` now
  passes `frmtmb_register_compat(expects = )`, which frmtmb gained in
  0.51.0, and an older frmtmb stops the load outright with
  `unused argument (expects = c("hmm", "lca"))`. The old floor of 0.46.0
  permitted an installation that could not load, so this is a hard
  requirement rather than a preference.

- `frmtmb.latent` joins `Suggests`. `hmm()` and `lca()` left frmtmb for
  that package, so the tests here that sample one of them qualify the
  call and skip when it is absent. The two compatibility rules this
  package registers for those pairs are declared in
  `frmtmb_register_compat(expects = )`, which says they name a feature
  another package owns. Before frmtmb checked, those rules were dropped
  without a word in every session that did not also load
  `frmtmb.latent`; now
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html)
  reports them as unresolved until it is loaded.

## frmtmb.sample 0.1.0

First release, extracted from frmtmb 0.46.0.

- NUTS sampling for frmtmb models through tmbstan:
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  on a fit or a formula, with brms 2.23’s default priors on both routes,
  priors stacking most-explicit-first, and non-centered sampling of
  qualifying random-effect blocks. Draws carry frmtmb parameter names.
- The full posterior method surface on core’s generics and on the brms,
  loo, posterior, bayesplot, rstantools, coda and bridgesampling
  generics:
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`loo()`](https://mc-stan.org/loo/reference/loo.html),
  [`waic()`](https://mc-stan.org/loo/reference/waic.html),
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  [`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html),
  [`conditional_effects()`](https://paulbuerkner.com/brms/reference/conditional_effects.brmsfit.html)
  on draws,
  [`bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html),
  [`mcmc_plot()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/draws-diagnostics.md),
  and the rest.
- [`check_laplace()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/check_laplace.md)
  measures the Laplace and Wald approximations against NUTS on the same
  objective, centered and with no added priors.
- At load, the package registers its sampling default priors with core’s
  [`get_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.html)
  and its feature rows with the compatibility matrix, so both stay
  truthful whether or not this package is installed.
- Refuses a tmbstan built against StanHeaders 2.39 or later, whose code
  generator leaves one log-density overload unpatched so every chain
  silently samples a standard normal instead of the model.
- Parallel chains work on every platform, Windows included: tmbstan
  rebuilds the tape on each socket worker from the serialized objective
  closure, and frmtmb’s generated closures are self-contained.
