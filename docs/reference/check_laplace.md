# Check the Laplace/Wald approximation against NUTS

Samples the fitted objective (see
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md))
and compares the ML estimates and sdreport standard errors against
posterior means and SDs. Close agreement supports the Laplace
approximation and Wald intervals; a posterior SD much larger than the
Wald SE, or a shifted mean, flags parameters where they are unreliable
(typically variance components with few groups).

## Usage

``` r
check_laplace(fit, chains = 2, iter = 1000, ...)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- chains, iter:

  Passed to
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Passed to
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

## Value

A data frame (one row per outer parameter) with columns `ml`,
`post_mean`, `wald_se`, `post_sd`, `z_shift` ((post_mean - ml)/post_sd)
and `sd_ratio` (post_sd/wald_se).
