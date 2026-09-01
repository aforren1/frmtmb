# Predictions from a frmtmb fit

Predictions from a frmtmb fit

## Usage

``` r
# S3 method for class 'frmtmb_fit'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "conditional", "zprob", "zlink", "disp"),
  dpar = NULL,
  resp = NULL,
  re.form = NULL,
  se.fit = FALSE,
  allow_new_levels = FALSE,
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- newdata:

  Optional data frame to predict on. Defaults to the training data.

- type:

  `"link"` for the linear predictor, `"response"` for the expected
  response (which equals
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on the
  training data; for zero-inflated, hurdle, and similar families this is
  the response mean, not the `mu` dpar). When `dpar` is given,
  `"response"` is that dpar on its natural scale. The glmmTMB spellings
  `"conditional"` (the `mu` dpar on its natural scale),
  `"zprob"`/`"zlink"` (the zero-inflation/hurdle probability on the
  response/link scale), and `"disp"` (the dispersion dpar) are accepted
  as aliases.

- dpar:

  Which distributional parameter to predict; defaults to the family's
  first location parameter (`"mu"` for most families).

- resp:

  For multivariate fits: which response to predict (defaults to the
  first).

- re.form:

  `NULL` (default) includes random effects; `NA` or `~0` gives
  population-level predictions.

- se.fit:

  If `TRUE`, return a list with elements `fit` and `se.fit`
  (delta-method standard errors accounting for fixed-effect and
  random-effect uncertainty). Exact `gp()` terms predict unseen
  positions by kriging: the conditional mean at the fitted kernel, with
  the GP conditional variance added to the standard errors.

- allow_new_levels:

  Predict unseen grouping-factor levels at the population level instead
  of erroring.

- ...:

  Unused.

## Value

A numeric vector, or a list when `se.fit = TRUE`.
