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
  `"oneStepGeneric"` otherwise.

- ...:

  For `type = "osa"`: passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

## Value

A numeric vector.
