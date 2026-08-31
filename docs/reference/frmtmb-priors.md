# Prior specifications for frm_sample

Priors apply on the INTERNAL parameter scale: coefficients are on their
link scale, and covariance parameters (`theta_*`) are the unconstrained
parameterization (log-SDs, scaled-Cholesky terms), so
`prior_normal(0, 1)` on `theta_1` is a lognormal prior on that SD.

## Usage

``` r
prior_normal(location = 0, scale = 1)

prior_t(df = 3, location = 0, scale = 1)
```

## Arguments

- location, scale, df:

  Prior parameters.

## Value

A `frmtmb_prior` object.
