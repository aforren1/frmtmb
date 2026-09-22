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
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(1)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  # each refusal names its reason and the replacement
  try(stancode(ds))
  try(standata(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1.1, 2.5)
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
#>   b                  (flat), as brms leaves slopes
#> Error : stancode() has no meaning for frmtmb: there is no Stan program. The model is an R closure built by build_objective() from the assembled frame and differentiated by RTMB, and the closure IS the source: print `ds$fit$obj$fn` for the evaluator and `ds$fit$frame` for everything baked into it
#> Error : standata() has no meaning for frmtmb: nothing is exported to a Stan data list. The assembled frame `ds$fit$frame` holds the same content (the response, the design matrices, the sparse Z, the addition terms), and model.matrix(), getME() and model.frame() read the pieces of it individually
# }
```
