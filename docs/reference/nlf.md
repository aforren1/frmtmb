# Add a nonlinear parameter formula to a model formula

brms's `nlf()`: it declares that one parameter is computed by a
NONLINEAR expression rather than by a linear predictor. The names in
that expression are either other model parameters, each of which needs
its own formula, or columns of the data.

## Usage

``` r
nlf(formula, ..., loop = NULL)
```

## Arguments

- formula:

  A two-sided formula naming the parameter on the left and its nonlinear
  body on the right, e.g. `sigma ~ a * exp(b * x)`.

- ...:

  Further two-sided formulas, treated as LINEAR parameter formulas
  exactly as if passed to
  [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md) - the brms
  convention.

- loop:

  Accepted for brms source compatibility and ignored. frmtmb evaluates a
  nonlinear body once over whole vectors, which is brms's
  `loop = FALSE`; a body built from elementwise operations has the same
  value either way.

## Value

An object of class `frmtmb_nlf`, to be added to a
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md).

## Details

`bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)` is the same model as
`bf(y ~ exp(b * x), b ~ 1, nl = TRUE)`, written the other way round.
Where `nl = TRUE` makes the response formula the nonlinear body, `nlf()`
names the parameter it belongs to, so any parameter can have one:
`bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1)` is a nonlinear
model for the residual standard deviation with a linear `mu`, which
`nl = TRUE` cannot spell.

Bodies may be chained: an `nlf()` body can name a parameter that another
`nlf()` defines, to any depth. The parameters are evaluated in
dependency order, and a cycle is refused by name. Each parameter's link
is applied to the body's value, so `nlf(sigma ~ a + b * z)` gives
`sigma = exp(a + b * z)` under the default log link, exactly as in brms.

A body may also name another distributional parameter of the same
response, which reads that parameter's per-row VALUE. This is the one
place frmtmb goes beyond brms, where such a name is a data column and
the model is refused when no column has it. It buys the variance
function of the model's own mean that nlme writes as
`varPower(form = ~ fitted(.))`: with the default log link on `sigma`,
`nlf(sigma ~ ls + th * log(abs(mu)))` is `sd = exp(ls) * |mu|^th`. A
column of the data still wins over the parameter name, so a body ported
from brms keeps its meaning.

Like [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md), an
`nlf()` in a multivariate model must be added to the
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) of the
response it belongs to, before the responses are combined.

## Examples

``` r
# the composed spelling of a nonlinear model
bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)
#> y ~ a
#> a ~ exp(b * x) (nonlinear)
#> b ~ 1 

# a nonlinear sigma with a linear mu
bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1)
#> y ~ x
#> sigma ~ a + b * z (nonlinear)
#> a ~ 1 
#> b ~ 1 

# bodies chain: cc feeds a, a feeds mu
bf(y ~ a, nl = TRUE) + nlf(a ~ cc * x) + nlf(cc ~ exp(b)) + lf(b ~ 1)
#> y ~ a (nonlinear)
#> a ~ cc * x (nonlinear)
#> cc ~ exp(b) (nonlinear)
#> b ~ 1 

# a variance function of the fitted mean: sd = exp(ls) * |mu|^th
bf(y ~ x) + nlf(sigma ~ ls + th * log(abs(mu))) + lf(ls ~ 1, th ~ 1)
#> y ~ x
#> sigma ~ ls + th * log(abs(mu)) (nonlinear)
#> ls ~ 1 
#> th ~ 1 
```
