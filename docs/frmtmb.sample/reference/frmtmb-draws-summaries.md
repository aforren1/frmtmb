# Summaries of the posterior predictive quantities

brms's three summarizing methods, each a summary of the draws method
beside it:

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) summarizes
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  (or
  [`posterior_linpred()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md)
  under `scale = "linear"`);

- [`predict()`](https://rdrr.io/r/stats/predict.html) summarizes
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/posterior_epred.md);

- [`residuals()`](https://rdrr.io/r/stats/residuals.html) summarizes
  [`predictive_error()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md).

`summary = FALSE` returns the draws themselves, which is what the method
it wraps returns.

## Usage

``` r
# S3 method for class 'frmtmb_draws'
fitted(
  object,
  newdata = NULL,
  re_formula = arg_unset(),
  scale = c("response", "linear"),
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)

# S3 method for class 'frmtmb_draws'
predict(
  object,
  newdata = NULL,
  re_formula = arg_unset(),
  transform = NULL,
  resp = NULL,
  negative_rt = FALSE,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ntrys = NULL,
  cores = NULL,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)

# S3 method for class 'frmtmb_draws'
residuals(
  object,
  newdata = NULL,
  re_formula = arg_unset(),
  method = "posterior_predict",
  type = c("ordinary", "pearson"),
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- newdata, re_formula, resp, dpar, nlpar, ndraws, draw_ids:

  Passed to the draws method, where they are documented.

- scale:

  `"response"` for the expected response, `"linear"` for the linear
  predictor.

- summary:

  If `FALSE`, the draws instead of their summary.

- robust:

  If `TRUE`, the median and MAD instead of the mean and standard
  deviation.

- probs:

  Probabilities of the quantile columns.

- ...:

  Passed to the draws method, which is where brms's `point_estimate` and
  `ndraws_point_estimate` are answered.

- transform:

  Applied to the draws before they are summarized.

- negative_rt, ntrys, cores, sort:

  brms's remaining arguments, carried so that a positional brms call
  asks the same question; each is passed on or refused by the method it
  belongs to.

- method:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html), which
  predictive distribution the error is taken against.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  `"ordinary"` (brms's name for the raw error) or `"pearson"`.

## Value

With `summary = TRUE` an observations-by-four matrix, or an
observations-by-four-by-K array for a category distribution. With
`summary = FALSE` the draws.
