# Residuals from a frmtmb fit

`"osa"` gives one-step-ahead (conditional quantile) residuals via
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html):
standard-normal under a correctly specified model, valid under
correlated observations where pearson residuals mislead.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
residuals(
  object,
  type = c("response", "pearson", "osa"),
  osa_method = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- type:

  `"response"`, `"pearson"`, or `"osa"`.

- osa_method:

  Method for
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html);
  defaults to `"fullGaussian"` for gaussian models and
  `"oneStepGeneric"` otherwise. A truncated, censored or ordinal
  response always uses `"oneStepGeneric"` (a truncated gaussian is not
  gaussian) with the integration domain and discrete support taken from
  the [`trunc()`](https://rdrr.io/r/base/Round.html) bounds or the
  censoring window, which must then be the same for every row.

- ...:

  For `type = "osa"`: passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

## Value

A numeric vector, `NA` on censored rows.

## Details

On a [`trunc()`](https://rdrr.io/r/base/Round.html)ed response,
`"response"` residuals are taken against the truncated mean
`E[Y | lb <= Y <= ub]`. `"pearson"` divides by the untruncated family
variance, so it is conservative there. `"osa"` builds its conditional
CDF on `[lb, ub]` (see `osa_method`).

On a `cens()`ed response, `"osa"` returns `NA` for every censored row:
what is observed there is an event (`Y > c`), not a value, and an event
has no one-step CDF. The uncensored rows get residuals conditional on
the censoring events, which needs one censoring point per side (type-I
censoring); row-varying censoring times and interval censoring are
refused. [`dharma_residuals()`](dharma_residuals.md) is not a substitute
on a censored fit, because
[`simulate.frmtmb_fit()`](simulate.frmtmb_fit.md) draws the latent
uncensored response and those draws are not comparable with the observed
censored values.

Ordinal responses use `"oneStepGeneric"` over the discrete support
`1..K`, which makes the residuals randomized quantile residuals.
