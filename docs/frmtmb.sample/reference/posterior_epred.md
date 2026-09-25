# Expected-value and predictive draws from sampled parameters

`posterior_epred()` evaluates the response-scale expectation per draw;
`posterior_predict()` additionally simulates responses from the family,
giving the posterior predictive distribution. Both condition on each
draw's own random effects (`re_formula = NA` drops them; `re.form` is an
accepted alias here, see *Argument spellings*).

## Usage

``` r
posterior_epred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_epred(
  object,
  newdata = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  point_estimate = NULL,
  ndraws_point_estimate = 1,
  ...
)

posterior_linpred(object, transform = FALSE, ...)

# S3 method for class 'frmtmb_draws'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  incl_thres = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  point_estimate = NULL,
  ndraws_point_estimate = 1,
  ...
)

posterior_predict(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_predict(
  object,
  newdata = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  transform = NULL,
  resp = NULL,
  negative_rt = FALSE,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ntrys = NULL,
  cores = NULL,
  point_estimate = NULL,
  ndraws_point_estimate = 1,
  ...
)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than a silently ignored name. The exceptions are brms's
  `allow_new_levels` (and `allow.new.levels`) and `sample_new_levels`.
  `allow_new_levels = FALSE`, and `TRUE` with levels the fit saw, answer
  as the call without it does. A level the fit did not see, in a term
  `re_formula` keeps, is refused with or without the flag, and so is a
  `newdata` that leaves the grouping column out under `TRUE`: brms draws
  that level's effect from each posterior draw, which is not built here,
  and predicting it at the population level would drop the group
  variance from every draw.
  [`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html)
  predicts unseen levels from the maximum-likelihood fit, and
  `re_formula = NA` predicts here at the population level.

- newdata, resp:

  As in
  [`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html).

- re_formula:

  The random-effect switch, in brms's spelling: `NULL` (the default)
  conditions on each draw's own random effects, `NA` or `~0` gives the
  population-level quantity. Its meaning is
  [`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html)'s
  `re_formula`; see *Argument spellings*.

- re.form:

  lme4's spelling of `re_formula`, accepted on the four methods where
  brms accepts it too. Pass one or the other, not both.

- dpar:

  Which distributional parameter to evaluate: its linear predictor for
  `posterior_linpred()`, its response-scale value for
  `posterior_epred()`. The default is the family's `mu`.

- nlpar:

  The parameter an `nlf()` body names. brms keeps it apart from `dpar`;
  frmtmb asks for either by the `dpar` name, so this is the same setting
  and the slot is here for brms's position.

- ndraws:

  Number of draws to use (default: all).

- draw_ids:

  The draws to use, by row index, instead of the evenly spaced subsample
  `ndraws` takes. Give one or the other.

- sort:

  brms's argument. Rows come back in the order of the data here, always,
  so `sort = TRUE` is refused by name.

- point_estimate, ndraws_point_estimate:

  brms's arguments, which collapse the draws to their `"mean"` or
  `"median"` FIRST and run the method once at that one parameter vector,
  repeated `ndraws_point_estimate` times. It is a parameter-space
  operation and not a summary of the output, as in brms.

- transform:

  For `posterior_linpred()`: if `TRUE`, apply the inverse link (the
  value of the `mu` dpar on its natural scale, brms's convention; unlike
  `posterior_epred()` this is not the response mean for zero-inflated
  and similar families).

- incl_thres:

  For `posterior_linpred()`: refused. brms subtracts a cumulative
  family's thresholds from the predictor; frmtmb returns the latent
  predictor itself.

- negative_rt:

  For `posterior_predict()`: refused. It is brms's sign convention for
  its own wiener family.

- ntrys, cores:

  brms's arguments, carried so that a positional brms call lands where
  brms lands it. Both are refused by name: the rejection limit of a
  [`trunc()`](https://rdrr.io/r/base/Round.html)ed draw is the family
  simulator's own, and the draws are replayed in one process.

## Value

A draws-by-observations matrix; for a categorical outcome
`posterior_epred()` returns a draws-by-observations-by-categories array
(see the section below).

## Categorical outcomes

An ordinal family predicts a DISTRIBUTION per observation, not one
number: each draw's `frm_linpred(type = "response")` is an `n x K`
matrix of category probabilities. Those stack into a 3-D
`draws x observations x categories` array. `dimnames` are
`list(NULL, <observation names or NULL>, <category levels>)`, so
`ep[, , "high"]` is the draws-by-observations matrix for one category
and `ep[k, , ]` is draw `k`'s own `n x K` prediction, the matrix
`frm_linpred(type = "response")` returns. Every `ep[k, i, ]` sums to 1
for an ordinal family.

This is brms's convention:
[`?brms::posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.html)
documents "an S x N x C array" for categorical and ordinal models and an
S x N matrix otherwise, and frmtmb follows brms spelling for brms-origin
functions. Any family whose per-draw response-scale prediction is a
matrix takes the array shape; every family that predicts one number per
observation keeps the plain `draws x observations` matrix.

`posterior_predict()` is unaffected for an ordinal or categorical family
(it draws one category per observation), and so is
`posterior_linpred()`, which is a statement about one distributional
parameter and stays an `n`-column matrix of the latent predictor. What
does take the array shape in `posterior_predict()` is a matrix-valued
RESPONSE:
[`frmtmb::multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
counts,
[`frmtmb::mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html)
draws and
[`frmtmb.latent::lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.html)
item codes give one row per observation, so the draws stack into
`draws x observations x columns`.

## Structured draws

`posterior_predict()` uses the same simulator
[`simulate()`](https://rdrr.io/r/stats/simulate.html) does, including
the structured families
([`frmtmb.latent::hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.html),
`mixture(groups = )`,
[`frmtmb::mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html))
and residual correlation terms; see the Structured draws section of
[`frmtmb::simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.html).
Those draws index the rows the model was fitted on, so `newdata` is
refused for them.

## Argument spellings

One rule decides every name in this package: where lme4 or glmmTMB and
brms disagree, brms wins. The random-effect switch is therefore
`re_formula` everywhere, on the draws methods here and on
[`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html)
and
[`frmtmb::simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.html)
alike; lme4's `re.form` was dropped from the fit surface and is refused
there by name.

Five methods take BOTH spellings, and they are exactly the five where
brms ITSELF accepts both. Four declare them: `posterior_epred()`,
`posterior_linpred()`, `posterior_predict()` and
[`predictive_error()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
carry `re_formula` and `re.form` side by side on `brmsfit`. The fifth
does not declare them and accepts them anyway:
`predictive_interval.brmsfit()`'s whole body is
`posterior_predict(object, ...)`, so the alias reaches a formal one
frame down. What brms ACCEPTS is the test, not what it declares.

[`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html) is
the one method that lost the alias, and the same test is why it stays
lost. brms DOES honor `re.form` there, through the same dots forwarding,
but it warns "unrecognized and ignored" while doing it. Matching brms
means matching what brms decided, and a warn-then-honor path is a leak
rather than a decision: the argument changes the answer and the message
says it did not.
[`predictive_interval()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
is the contrast, where brms honors the alias silently and this package
follows.

Giving both at once is refused rather than resolved. Two names for one
setting supplied together is a question about what was meant, and
guessing at it would silently ignore one of them.

The argument ORDER is brms's too, so a positional brms call means the
same thing here: `newdata` then `re_formula` then `re.form` then `resp`,
after `transform` in `posterior_linpred()` and before it in
`posterior_predict()`.

The literal default of both formals is an internal "not supplied" marker
rather than a value, because `NULL` (keep the random effects) and `NA`
(drop them) are both real settings here and neither can double as
"unset". The behavior when neither is given is unchanged: `NULL` on
every draws method, `NA` on
[`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
for a fit.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rpois(80, exp(0.3 + 0.4 * dd$x + rnorm(8, 0, 0.5)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

nd <- data.frame(x = c(-1, 0, 1),
                 g = factor(1, levels = levels(dd$g)))

# the expected response per draw: uncertainty in the mean
ep <- posterior_epred(ds, newdata = nd)
apply(ep, 2, quantile, c(0.025, 0.5, 0.975))

# the predictive distribution adds the family's own noise, so its
# intervals are wider
pp <- posterior_predict(ds, newdata = nd)
apply(pp, 2, quantile, c(0.025, 0.5, 0.975))

# the linear predictor itself, on the link scale by default
head(posterior_linpred(ds, newdata = nd, ndraws = 5))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0, 2.5)
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: The largest R-hat is 1.13, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>             [,1]      [,2]      [,3]
#> [1,]  0.48954070 0.6555560 0.8215714
#> [2,]  0.33465244 0.7663111 1.1979698
#> [3,] -0.34015303 0.1554798 0.6511126
#> [4,]  0.08794828 0.3320401 0.5761320
#> [5,]  0.28325834 0.5987430 0.9142276
# }
```
