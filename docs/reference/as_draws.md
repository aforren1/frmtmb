# Convert to a posterior draws object

The `as_draws` family converts an object holding posterior draws into
one of the posterior package's draws formats. frmtmb defines the
generics so that they work whether or not posterior is attached, and
registers methods with posterior so that its own spellings dispatch too.

## Usage

``` r
as_draws(x, ...)

as_draws_matrix(x, ...)

as_draws_array(x, ...)

as_draws_df(x, ...)

as_draws_list(x, ...)

as_draws_rvars(x, ...)

# S3 method for class 'frmtmb_fit'
as_draws(x, ...)

# S3 method for class 'frmtmb_fit'
as_draws_matrix(x, ...)

# S3 method for class 'frmtmb_fit'
as_draws_array(x, ...)

# S3 method for class 'frmtmb_fit'
as_draws_df(x, ...)

# S3 method for class 'frmtmb_fit'
as_draws_list(x, ...)

# S3 method for class 'frmtmb_fit'
as_draws_rvars(x, ...)

# S3 method for class 'frmtmb_multiple'
as_draws(x, ...)
```

## Arguments

- x:

  An object holding draws.

- ...:

  Passed to methods.

## Value

A posterior draws object of the requested format.

## Details

Core has no object that carries draws:
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) is maximum
likelihood and
[`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
pools point estimates, so its methods explain that rather than inventing
a draws matrix. Install `frmtmb.sample` and sample with
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
to get an object these convert.

## Examples

``` r
# a maximum-likelihood fit carries no draws, and says so by name
dd <- data.frame(y = rnorm(40), x = rnorm(40))
try(as_draws_df(frm(bf(y ~ x) + gaussian(), data = dd)))
#> Error : as_draws_df() needs posterior draws and a frmtmb_fit has none: frm() is maximum likelihood, so it carries one parameter vector rather than chains. Install frmtmb.sample and sample first, with as_draws_df(frmtmb.sample::frm_sample(fit)); for the point estimates and their covariance use fixef(), vcov() or confint() on the fit itself

# frm_multiple() pools estimates rather than carrying draws, so it
# answers with the reason rather than a matrix
fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
try(as_draws(fits))
#> Error : as_draws() needs draws, and a frm_multiple() result has none: it is m maximum-likelihood fits pooled by Rubin's rules, with no chains. Read the pooled tables from `x$pooled` and `x$pooled_varcorr` or test with hypothesis(), and use frm_sample() on one imputation's fit (`x$fits[[1]]`) for draws
```
