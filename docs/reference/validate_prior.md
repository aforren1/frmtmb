# Check a prior against a model

brms's `validate_prior()`: resolve a prior specification against the
model it is meant for, refuse it if any part of it addresses nothing the
model has, and return the whole
[`default_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.md)
table with the prior filled in. Fitting runs the same check, so a prior
this function accepts is one
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) accepts,
and a prior it refuses is refused here with the message the fit would
give, before any fitting work.

## Usage

``` r
validate_prior(
  prior,
  formula,
  data,
  family = NULL,
  data2 = list(),
  route = c("fit", "sample")
)
```

## Arguments

- prior:

  A `frmtmb_priorlist` from
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md),
  a `brmsprior`, or `NULL` for none.

- formula:

  A [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula
  (with family) or a plain formula.

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

## Value

A `frmtmb_prior_rows` table, as
[`default_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.md)
returns.

## Details

A row the prior sets reads `source = "user"`. A row that takes its
density from a class row above it, the coefficients under a class `"b"`
prior or the blocks under a class `"sd"` prior, reads
`source = "(vectorized)"` and repeats that density, as brms prints it.
Every other row keeps its default.

Two specifications for the same slot are refused, as
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) and brms
refuse them. Where specifications for different slots reach one
parameter, the table shows the more specific one, which is the one
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) applies;
see
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).

The table can be passed back as `prior =`:
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) reads a
`frmtmb_prior_rows` table the way
[`as.brmsprior()`](https://aforren1.github.io/frmtmb/reference/as.brmsprior.md)
does, applying the rows with a density or a bound and skipping the flat
and vectorized ones.

## Examples

``` r
dd <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60),
                 g = factor(rep(1:6, 10)))
validate_prior(prior(normal(0, 10), class = b) +
                 prior(cauchy(0, 2), class = sd),
               y ~ x + z + (1 | g), data = dd)
#> route = "fit": the prior defaults frm() applies
#>           prior     class    coef group resp dpar nlpar lb ub       source
#> 1        (flat) Intercept                               NA NA      default
#> 2 normal(0, 10)         b                               NA NA         user
#> 3 normal(0, 10)         b       x                       NA NA (vectorized)
#> 4 normal(0, 10)         b       z                       NA NA (vectorized)
#> 5        (flat)     sigma                               NA NA      default
#> 6  cauchy(0, 2)        sd                               NA NA         user
#> 7  cauchy(0, 2)        sd             g                 NA NA (vectorized)
#> 8        (flat)     theta                               NA NA      default
#> 9        (flat)     theta theta_1                       NA NA      default
# a prior on a coefficient the model does not have is refused
try(validate_prior(prior(normal(0, 1), coef = w), y ~ x, data = dd))
#> Error : Prior target not found (class=b, coef=w)
```
