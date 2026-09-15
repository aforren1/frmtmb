# Size of a draws object

`ndraws()` counts the post-warmup draws (all chains pooled),
`niterations()` the draws per chain, `nchains()` the chains and
`nvariables()` the sampled parameters. The names and meanings are
posterior's; frmtmb registers methods with posterior so that the
generics work whether or not that package is attached.

## Usage

``` r
ndraws(x)

nchains(x)

niterations(x)

nvariables(x, ...)

# S3 method for class 'frmtmb_fit'
ndraws(x)

# S3 method for class 'frmtmb_fit'
nchains(x)

# S3 method for class 'frmtmb_fit'
niterations(x)

# S3 method for class 'frmtmb_fit'
nvariables(x, ...)
```

## Arguments

- x:

  An object holding draws.

- ...:

  Unused. Carried because posterior's `nvariables()` generic has it and
  frmtmb hands these generics back to posterior when posterior is
  loaded.

## Value

A single integer.

## Details

As with the `as_draws` family, the objects that answer these with a
number come from
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html).

## Examples

``` r
dd <- data.frame(y = rnorm(40), x = rnorm(40))
fits <- frm_multiple(bf(y ~ x) + gaussian(), data = list(dd, dd))
try(ndraws(fits))
#> Error : ndraws() needs draws, and a frm_multiple() result has none: it is m maximum-likelihood fits pooled by Rubin's rules, with no chains. Read the pooled tables from `x$pooled` and `x$pooled_varcorr` or test with hypothesis(), and use frm_sample() on one imputation's fit (`x$fits[[1]]`) for draws
```
