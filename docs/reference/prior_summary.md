# Priors used in a fit

Priors used in a fit

## Usage

``` r
prior_summary(object, ...)

# S3 method for class 'frmtmb_fit'
prior_summary(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

The `frmtmb_priorlist` the fit was penalized with, or (invisibly) `NULL`
for a plain maximum-likelihood fit.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# what was actually applied, after set_prior() was matched to the
# coefficients of this design
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           priors = set_prior("normal(0, 1)", class = "b") +
                    set_prior("exponential(1)", class = "sd"))
prior_summary(fit)
#> normal(0, 1) class=b
#> exponential(1) class=sd

# a plain maximum-likelihood fit reports that it had none
prior_summary(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd))
#> No priors were set (plain maximum likelihood).
```
