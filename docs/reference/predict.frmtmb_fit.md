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

## Details

When the fixed-effect design was rank deficient, the aliased columns
were dropped at fit time and some coefficient combinations are not
estimable. Rows of `newdata` that load on a dropped direction get `NA`
(and `NA` standard errors), with one warning naming the dropped columns;
every other row is unaffected. The test is the one
[`stats::predict.lm()`](https://rdrr.io/r/stats/predict.lm.html) uses: a
row is non-estimable when it is not orthogonal to the null space of the
fitted design, up to a relative tolerance of `1e-8`. Two limits follow.
It is a numerical test, so near-aliased designs sit on a threshold
rather than a clean yes/no. And it covers the parametric fixed-effect
block only: smooth null-space, `gp()`, `mo()` and `mi()` columns are
appended after the rank check and are never dropped.

## Truncated responses

For a response with [`trunc()`](https://rdrr.io/r/base/Round.html)
bounds, `type = "response"` (and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)) report the
truncated mean `E[Y | lb <= Y <= ub]`, matching the likelihood the model
was fitted with. Predictions of a distributional parameter
(`type = "link"`, `dpar = `, or `type = "conditional"`) stay
**untruncated**: they are statements about the latent parameter, not
about the observed, truncated response. Bounds are re-evaluated on
`newdata` the same way `trials()` and `se()` are: a literal bound
carries over unchanged, and a bound given as a variable must be a column
of `newdata` of the right length.

## Standard errors of the expected response

For a family whose mean is the `mu` dpar, `se.fit` on
`type = "response"` is the usual one-predictor delta method:
`|dmu/deta| * se(eta)`.

When the mean is a function of several dpars (zero-inflated and hurdle
families, `lognormal`, a `trials()` binomial, or any
[`trunc()`](https://rdrr.io/r/base/Round.html)ed response), the delta
method runs jointly over every dpar's linear predictor: `se^2 = g' V g`,
where row `i` of `g` stacks `dm_i/deta_k` times the design row of
predictor `k`, and `V` is the joint covariance of all the coefficients
([`vcov()`](https://rdrr.io/r/stats/vcov.html)'s `jointPrecision` block,
so the cross-predictor covariances and the shared random-effect block
are included). The gradients `dm/deta_k` are central differences of the
family mean, taken one predictor at a time with a relative step.

Random effects enter conditional on their modes, the same convention
`se.fit` uses for the linear predictor. Unseen grouping levels
(`allow_new_levels = TRUE`) add their block's marginal variance,
propagated through the same gradients.
