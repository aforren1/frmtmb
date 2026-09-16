# Summaries and intervals of draws

[`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html)
reduces draws to estimate, error and quantiles in brms's column layout
(`Estimate`, `Est.Error`, `Q2.5`, `Q97.5`); `posterior_interval()` gives
the central interval alone, in rstantools' layout. Both work on a
`frmtmb_draws` object and on any matrix of draws, which is what makes
`posterior_summary(bayes_R2(ds, summary = FALSE))` work.

## Usage

``` r
# S3 method for class 'frmtmb_draws'
posterior_summary(
  x,
  probs = c(0.025, 0.975),
  robust = FALSE,
  variable = NULL,
  ...
)

posterior_interval(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_interval(
  object,
  pars = NA,
  variable = NULL,
  prob = 0.95,
  regex = FALSE,
  fixed = FALSE,
  ...
)

predictive_interval(object, ...)

# S3 method for class 'frmtmb_draws'
predictive_interval(
  object,
  prob = 0.9,
  newdata = NULL,
  resp = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  ndraws = NULL,
  ...
)

predictive_error(object, ...)

# S3 method for class 'frmtmb_draws'
predictive_error(
  object,
  newdata = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  method = "posterior_predict",
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)
```

## Arguments

- x:

  The same, for
  [`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html),
  whose generic is brms's and names its first argument `x`.

- probs:

  Quantiles for
  [`posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html).

- robust:

  If `TRUE`, median and MAD instead of mean and SD.

- variable:

  Optional subset of variables, by name.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

- object:

  A `frmtmb_draws`, or a matrix of draws (variables in columns).

- pars:

  brms's alias of `variable`, in brms's own second position on
  `posterior_interval()`: `NA` (the default) for every variable,
  otherwise a character vector matched as a regular expression unless
  `fixed = TRUE`. brms refuses a `pars` that is neither `NA` nor
  character, and so does this, which is why `posterior_interval(x, 0.9)`
  is a refusal and not an interval.

- prob:

  Central interval width for `posterior_interval()` and
  `predictive_interval()`.

- regex:

  If `TRUE`, `variable` is a regular expression.

- fixed:

  If `TRUE`, `pars` is matched by exact name.

- re_formula, re.form:

  Passed to
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md),
  which takes brms's `re_formula` and accepts lme4's `re.form` as an
  alias of it. Pass one or the other; see the *Argument spellings*
  section of
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md).

- ndraws, draw_ids, newdata, resp:

  Passed to
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md).
  `predictive_error(newdata =)` re-evaluates the response term on
  `newdata`, so `newdata` must carry the response.

- method:

  For `predictive_error()`, which predictive draws the error is taken
  against: `"posterior_predict"` (the default) or `"posterior_epred"`.

## Value

A matrix with one row per variable (or per observation, for the
predictive functions), except `predictive_error()`, which returns a
draws-by-observations matrix.

## Details

`predictive_interval()` is the same central interval of
[`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
draws, and `predictive_error()` is the matrix of predictive residuals
`y - yrep`, one row per draw.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  posterior_summary(ds, variable = c("Intercept", "x"))
  posterior_interval(ds, prob = 0.9, variable = "x")
  head(predictive_interval(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>              5%      95%
#> [1,] -1.0121251 1.978241
#> [2,] -1.3895099 1.570833
#> [3,] -0.5332048 2.593338
#> [4,] -1.0027390 2.427115
#> [5,] -0.5427595 2.690614
#> [6,] -1.7351044 1.568193
# }
```
