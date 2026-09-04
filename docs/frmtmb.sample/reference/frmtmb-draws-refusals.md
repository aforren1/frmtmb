# Methods a ported brms script may call that frmtmb does not have

These `brmsfit` methods either describe machinery frmtmb does not use
(Stan code and Stan data) or are brms spellings that have been renamed.
They are defined so that a ported script gets the reason and the
replacement rather than "could not find function", which is what the
vignette-port audit measured most of its post-processing failures as.

## Usage

``` r
stancode(object, ...)

# S3 method for class 'frmtmb_draws'
stancode(object, ...)

standata(object, ...)

# S3 method for class 'frmtmb_draws'
standata(object, ...)

# S3 method for class 'frmtmb_draws'
expose_functions(x, ...)

# S3 method for class 'frmtmb_draws'
plot(x, ...)

# S3 method for class 'frmtmb_draws'
update(object, ...)

restructure(x, ...)

# S3 method for class 'frmtmb_draws'
restructure(x, ...)

posterior_samples(x, ...)

# S3 method for class 'frmtmb_draws'
posterior_samples(x, ...)

nsamples(object, ...)

# S3 method for class 'frmtmb_draws'
nsamples(object, ...)

parnames(x, ...)

# S3 method for class 'frmtmb_draws'
parnames(x, ...)
```

## Arguments

- object, x, ...:

  Ignored; these functions always stop.

## Value

These functions never return; they signal an error.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  # each refusal names its reason and the replacement
  try(stancode(ds))
  try(nsamples(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1.1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Error : stancode() has no meaning for frmtmb: there is no Stan program. The model is an R closure built by build_objective() from the assembled frame and differentiated by RTMB, and the closure IS the source: print `ds$fit$obj$fn` for the evaluator and `ds$fit$frame` for everything baked into it
#> Error : nsamples() is the deprecated brms spelling. Use ndraws(x) for the pooled draw count, niterations(x) for the per-chain count
# }
```
