# Expected-value and predictive draws from sampled parameters

`posterior_epred()` evaluates the response-scale expectation per draw;
`posterior_predict()` additionally simulates responses from the family,
giving the posterior predictive distribution. Both condition on each
draw's own random effects (`re.form = NA` drops them).

## Usage

``` r
posterior_epred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_epred(
  object,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  ndraws = NULL,
  ...
)

posterior_linpred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  dpar = NULL,
  ndraws = NULL,
  ...
)

posterior_predict(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_predict(
  object,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Unused.

- newdata, resp, re.form:

  As in
  [`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md).

- ndraws:

  Number of draws to use (default: all).

- transform:

  For `posterior_linpred()`: if `TRUE`, apply the inverse link (the
  value of the `mu` dpar on its natural scale, brms's convention; unlike
  `posterior_epred()` this is not the response mean for zero-inflated
  and similar families).

- dpar:

  For `posterior_linpred()`: which distributional parameter's linear
  predictor to evaluate.

## Value

A draws-by-observations matrix.
