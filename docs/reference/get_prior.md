# Enumerate the targetable prior slots

The
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
counterpart of brms's `get_prior()`: one row per slot a prior can
target, with the class/coef/dpar/group values to pass to
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).
The default in every slot is flat (this is maximum likelihood until
priors are set; the formula route of
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
has its own brms defaults, which
[`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md)
reports). Classes `"sd"` and `"cor"` are targeted by `group` and
`nlpar`; class `"theta"` rows name the raw internal covariance
parameters (escape hatch, including correlations one at a time).

## Usage

``` r
get_prior(formula, data = NULL, family = NULL, data2 = list())
```

## Arguments

- formula:

  A [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula
  (with family), a plain formula, or an already fitted `frmtmb_fit`.

- data:

  A data frame of model data (ignored when `formula` is a fit).

- family:

  Family, when `formula` does not carry one.

- data2:

  Structural objects, as in
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) (ignored
  when `formula` is a fit, which carries its own).

## Value

A data frame with columns `prior`, `class`, `coef`, `group`, `dpar`,
`nlpar`, `resp`, `lb`, `ub`.

## Details

A nonlinear parameter's coefficients are listed under class `"b"` with
its name in the `nlpar` column, the intercept among them, which is how
brms lists them and what
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
addresses (see its Nonlinear parameters section).

## Examples

``` r
dd <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(1:6, 10)))
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#>    prior     class    coef group  dpar nlpar resp lb ub
#> 1 (flat) Intercept                                NA NA
#> 2 (flat)         b                                NA NA
#> 3 (flat)         b       x                        NA NA
#> 4 (flat) Intercept               sigma            NA NA
#> 5 (flat)        sd                                NA NA
#> 6 (flat)        sd             g                  NA NA
#> 7 (flat)     theta                                NA NA
#> 8 (flat)     theta theta_1                        NA NA
```
