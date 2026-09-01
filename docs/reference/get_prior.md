# Enumerate the targetable prior slots

The
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
counterpart of brms's `get_prior()`: one row per slot a prior can
target, with the class/coef/dpar/group values to pass to
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).
The default in every slot is flat (this is maximum likelihood until
priors are set). Class `"sd"` is targeted by `group` only; class
`"theta"` rows name the raw internal covariance parameters (escape
hatch, including correlations).

## Usage

``` r
get_prior(formula, data = NULL, family = NULL)
```

## Arguments

- formula:

  A [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula
  (with family), a plain formula, or an already fitted `frmtmb_fit`.

- data:

  Model data (ignored when `formula` is a fit).

- family:

  Family, when `formula` does not carry one.

## Value

A data frame with columns `prior`, `class`, `coef`, `group`, `dpar`,
`resp`, `lb`, `ub`.

## Examples

``` r
dd <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(1:6, 10)))
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#>    prior     class    coef group  dpar resp lb ub
#> 1 (flat) Intercept                          NA NA
#> 2 (flat)         b                          NA NA
#> 3 (flat)         b       x                  NA NA
#> 4 (flat) Intercept               sigma      NA NA
#> 5 (flat)        sd                          NA NA
#> 6 (flat)        sd             g            NA NA
#> 7 (flat)     theta                          NA NA
#> 8 (flat)     theta theta_1                  NA NA
```
