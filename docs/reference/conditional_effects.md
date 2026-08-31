# Conditional effects of predictors

For each requested effect, predicts over a grid of that predictor with
every other predictor held at a reference value (numeric: mean; factor:
first level; matrix covariate: column means) and random effects excluded
(`re.form = NA`). Confidence bands are Wald intervals computed on the
link scale and back-transformed. Smooth terms are included, so this also
covers what brms calls `conditional_smooths()`.

## Usage

``` r
conditional_effects(x, ...)

# S3 method for class 'frmtmb_fit'
conditional_effects(
  x,
  effects = NULL,
  resp = NULL,
  dpar = NULL,
  resolution = 100,
  prob = 0.95,
  conditions = list(),
  data = NULL,
  ...
)
```

## Arguments

- x:

  A `frmtmb_fit`.

- ...:

  Passed to [`predict.frmtmb_fit()`](predict.frmtmb_fit.md).

- effects:

  Character vector of variable names, or `"x:z"` pairs; for a pair, the
  first variable is varied over its range while the second is held at
  its levels (factors) or at mean and mean plus or minus one SD
  (numeric). Default: every fixed-effect and smooth variable of the
  selected linear predictor.

- resp, dpar:

  Response and distributional parameter, as in
  [`predict.frmtmb_fit()`](predict.frmtmb_fit.md).

- resolution:

  Number of grid points for a varied numeric predictor.

- prob:

  Coverage of the confidence bands (brms spelling).

- conditions:

  Named list overriding reference values, e.g. `list(x2 = 1, g = "b")`.

- data:

  The original model data. Only needed when the model frame does not
  store a raw variable (e.g. a variable used only inside
  [`poly()`](https://rdrr.io/r/stats/poly.html)).

## Value

A named list of data frames (one per effect) with the varied variable(s)
plus `estimate__`, `se__` (link scale), `lower__`, and `upper__`;
printing it draws the plots.
