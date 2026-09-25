# Default priors: the slots a prior can target

brms's `default_prior()`, and `get_prior()`, which brms has kept as its
alias since 2.20.14: one row per slot a prior can target, with the
class/coef/dpar/group values to pass to
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).
Classes `"sd"` and `"cor"` are targeted by `group` and `nlpar`; the
residual-correlation classes (`"ar"`, `"ma"`, `"cosy"`, `"cortime"`) by
`resp`; class `"rescor"` has one row and takes no `resp`; and class
`"theta"` rows name the raw internal covariance parameters (escape
hatch, including correlations one at a time, across all three covariance
components). Where the location is several distributional parameters (a
categorical or mixture model), each one's `b`, `Intercept` and `sd` rows
carry its `dpar`, as in brms.

## Usage

``` r
default_prior(
  object,
  data = NULL,
  family = NULL,
  data2 = list(),
  route = c("fit", "sample")
)

get_prior(formula, ...)
```

## Arguments

- object:

  A [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula
  (with family), a plain formula, or an already fitted `frmtmb_fit`.

- data:

  A data frame of model data (ignored when `object` is a fit).

- family:

  Family, when `object` does not carry one.

- data2:

  Structural objects, as in
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) (ignored
  when `object` is a fit, which carries its own).

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

- formula:

  The same as `object`, under the name brms's `get_prior()` gives it.

- ...:

  For `get_prior()`, the arguments of `default_prior()`.

## Value

A data frame of class `frmtmb_prior_rows` with brms's columns `prior`,
`class`, `coef`, `group`, `resp`, `dpar`, `nlpar`, `lb`, `ub` and
`source`, and a `route` attribute. `source` is `"default"` on every row;
[`validate_prior()`](https://aforren1.github.io/frmtmb/reference/validate_prior.md)
marks the rows a user set.

## Details

A nonlinear parameter's coefficients are listed under class `"b"` with
its name in the `nlpar` column, the intercept among them, which is how
brms lists them and what
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
addresses (see its Nonlinear parameters section).

A distributional parameter is listed under whichever of its two
spellings the model offers: its OWN class where the model gives it no
predictor, and class `"Intercept"` with its name in the `dpar` column
where the model gives it a formula. brms lists the same two, the same
way round. See A distributional parameter's own class in
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).

An ordinal family has no intercept column, so its class `"Intercept"`
row names the THRESHOLD vector, which is what the same row means in
brms. See the Ordinal thresholds section of
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).

## Where the table differs from brms's

The class `"theta"` rows are frmtmb's own. They name the internal
covariance parameters, which are real parameters of the fit, and for a
`gp()` length-scale, a `car()` dependence parameter or one entry of a
structured covariance they are the only spelling that reaches one
parameter. brms has no such class, so a brms table for the same model
has fewer rows.

brms lists a class `"sd"` row per coefficient of each block
(`sd_patient__Intercept`). frmtmb's class `"sd"` addresses a whole block
and refuses a `coef`, so those rows are not listed; class `"theta"`
reaches one standard deviation on its own.

A flat slot reads `"(flat)"` in the `prior` column, where brms stores an
empty string and prints `(flat)`.

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

brms's `default_prior()` describes what `brm()` would use, so the brms
reading of this function is `route = "sample"`. `route = "fit"` has no
brms counterpart: it describes
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

With frmtmb.sample attached, the example below with `route = "sample"`
added returns the same rows as the example below with brms's densities
in the `prior` column instead of `(flat)`: a Student-t on the intercept
centered on the response, a half-Student-t on `sigma` and on each
standard deviation, and `lkj(1)` on each correlation. Population-level
slopes stay `(flat)`, as they are in brms.

## See also

[`validate_prior()`](https://aforren1.github.io/frmtmb/reference/validate_prior.md)
for the table with a prior filled in.

## Examples

``` r
dd <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(1:6, 10)))
# what frm() applies: flat, whatever else is loaded. For what
# frm_sample() applies, see "Which route the defaults describe"
default_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#> route = "fit": the prior defaults frm() applies
#>    prior     class    coef group resp dpar nlpar lb ub  source
#> 1 (flat) Intercept                               NA NA default
#> 2 (flat)         b                               NA NA default
#> 3 (flat)         b       x                       NA NA default
#> 4 (flat)     sigma                               NA NA default
#> 5 (flat)        sd                               NA NA default
#> 6 (flat)        sd             g                 NA NA default
#> 7 (flat)     theta                               NA NA default
#> 8 (flat)     theta theta_1                       NA NA default
# the same table under brms's older name
get_prior(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#> route = "fit": the prior defaults frm() applies
#>    prior     class    coef group resp dpar nlpar lb ub  source
#> 1 (flat) Intercept                               NA NA default
#> 2 (flat)         b                               NA NA default
#> 3 (flat)         b       x                       NA NA default
#> 4 (flat)     sigma                               NA NA default
#> 5 (flat)        sd                               NA NA default
#> 6 (flat)        sd             g                 NA NA default
#> 7 (flat)     theta                               NA NA default
#> 8 (flat)     theta theta_1                       NA NA default
```
