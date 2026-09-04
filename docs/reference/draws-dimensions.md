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

nvariables(x)
```

## Arguments

- x:

  An object holding draws.

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
#> Error in UseMethod("ndraws") : 
#>   no applicable method for 'ndraws' applied to an object of class "frmtmb_multiple"
```
