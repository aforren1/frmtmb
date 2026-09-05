# Enumerate the targetable prior slots

The
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
counterpart of brms's `get_prior()`: one row per slot a prior can
target, with the class/coef/dpar/group values to pass to
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).
Classes `"sd"` and `"cor"` are targeted by `group` and `nlpar`; the
residual-correlation classes (`"ar"`, `"ma"`, `"cosy"`, `"cortime"`,
`"rescor"`) by `resp`; and class `"theta"` rows name the raw internal
covariance parameters (escape hatch, including correlations one at a
time, across all three covariance components).

## Usage

``` r
get_prior(
  formula,
  data = NULL,
  family = NULL,
  data2 = list(),
  route = c("fit", "sample")
)
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

- route:

  Which route's defaults the `prior` column reports. `"fit"` (the
  default) reports the defaults
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) applies,
  which are flat in every slot; it consults no registry, so its answer
  does not depend on which packages are attached. `"sample"` reports the
  defaults
  [`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
  applies, and refuses when no package has registered any. The returned
  object records the route and
  [`print()`](https://rdrr.io/r/base/print.html) names it on its first
  line.

## Value

A data frame of class `frmtmb_prior_rows` with columns `prior`, `class`,
`coef`, `group`, `dpar`, `nlpar`, `resp`, `lb`, `ub`, and a `route`
attribute.

## Details

A nonlinear parameter's coefficients are listed under class `"b"` with
its name in the `nlpar` column, the intercept among them, which is how
brms lists them and what
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
addresses (see its Nonlinear parameters section).

## Which route the defaults describe

Every column but `prior` is a property of the design, and the design
does not change with what is attached. The `prior` column is a property
of a ROUTE, and frmtmb has two of them with different defaults:
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) is maximum
likelihood and is flat in every slot until a prior is set, while
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
applies brms's weakly-informative defaults on both of its routes.
`route` makes the caller say which one is being asked about, so that the
answer is a property of the question rather than of the search path.

`route = "fit"`, the default, reads no registry at all. Its table is
identical whatever extension packages are loaded.

`route = "sample"` reads the defaults `frm_sample()` would apply, which
only frmtmb.sample can state. Without that package loaded the call is
refused rather than answered `(flat)`, because a flat table would be a
wrong answer about the sampling route and not a missing one.

brms's `get_prior()` describes what `brm()` would use, so the brms
reading of this function is `route = "sample"`. `route = "fit"` has no
brms counterpart: it describes
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

With frmtmb.sample attached,
`get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd, route = "sample")`
returns the same rows as the example below with brms's densities in the
`prior` column instead of `(flat)`: a Student-t on the intercept
centered on the response, a half-Student-t on `sigma` and on each
standard deviation, and `lkj(1)` on each correlation. Population-level
slopes stay `(flat)`, as they are in brms.

## Examples

``` r
dd <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(1:6, 10)))
# what frm() applies: flat, whatever else is loaded. For what
# frm_sample() applies, see "Which route the defaults describe"
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#> route = "fit": the prior defaults frm() applies
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
