# Model structure behind a set of draws

[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`formula()`](https://rdrr.io/r/stats/formula.html),
[`family()`](https://rdrr.io/r/stats/family.html),
[`getCall()`](https://rdrr.io/r/stats/update.html) and
[`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.html) on a
`frmtmb_draws` report the model the sampler ran, by delegating to the
fit stored inside it. They read structure only, so they work on draws
from the formula route, which has no maximum-likelihood estimate.

## Usage

``` r
# S3 method for class 'frmtmb_draws'
nobs(object, ...)

# S3 method for class 'frmtmb_draws'
formula(x, ...)

# S3 method for class 'frmtmb_draws'
family(object, ...)

# S3 method for class 'frmtmb_draws'
getCall(x, ...)

# S3 method for class 'frmtmb_draws'
coef(object, ...)
```

## Arguments

- object, x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- ...:

  Unused.

## Value

As for the corresponding `frmtmb_fit` method.

## Details

[`coef()`](https://rdrr.io/r/stats/coef.html) is a posterior quantity,
not a structural one: it summarizes the per-group coefficients (fixed
effects plus that group's own random effects) over the draws, in the
same nested shape
[`frmtmb::coef.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/coef.frmtmb_fit.html)
returns, with a `levels x statistics x coefficients` array in place of
each data frame. That is brms's `coef.brmsfit` layout.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  nobs(ds)
  ngrps(ds)
  coef(ds)$g[1:3, , "(Intercept)"]
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#>    Estimate Est.Error        Q2.5    Q97.5
#> 1 0.8430866 0.2615370  0.37389391 1.354788
#> 2 0.5154862 0.2922578 -0.06618975 0.965672
#> 3 0.9697659 0.2697711  0.54900417 1.565481
# }
```
