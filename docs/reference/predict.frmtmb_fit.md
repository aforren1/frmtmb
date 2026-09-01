# Predictions from a frmtmb fit

Predictions from a frmtmb fit

## Usage

``` r
# S3 method for class 'frmtmb_fit'
predict(
  object,
  newdata = NULL,
  type = c("link", "response"),
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

  `"link"` for the linear predictor, `"response"` for the dpar on its
  natural scale.

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
