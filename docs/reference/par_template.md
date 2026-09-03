# Parameter names and starting values

The parameter vector a model optimizes, as an editable named list: one
component per parameter block (`beta`, `betad`, `theta`, ...), each a
named numeric vector. It answers "what do I call these?" before there is
a fit to ask, which is what `frm(start =)`,
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)'s
`newparams =`, and `frm(lower =)`/`upper =` all need.

## Usage

``` r
par_template(object, ...)

# S3 method for class 'frmtmb_fit'
par_template(object, ...)

# Default S3 method
par_template(
  object,
  data,
  family = NULL,
  start = NULL,
  prior = NULL,
  na.action = stats::na.omit,
  data2 = list(),
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`, or a
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula or
  plain formula.

- ...:

  Unused.

- data:

  Model data. Required when `object` is a formula.

- family:

  Family, when `object` does not carry one.

- start:

  Optional `start` list, applied exactly as
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) would,
  so the result shows the starting values a fit with that `start` uses.
  Formula method only.

- prior:

  Optional
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  specification, so the result shows the nonlinear starting values its
  locations place. Formula method only.

- na.action, data2:

  As in [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).
  Formula method only.

## Value

A `frmtmb_par_template`: a named list of named numeric vectors, accepted
directly as `frm(start =)` and as `frm_simulate(newparams =)`.

## Details

On a fitted model the values are the estimates. On a formula and data
the frame is assembled but nothing is fitted, and the values are the
ones [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) would
start from, `start =` and `prior =` included. Either way the result can
be edited and passed straight back:

    st <- par_template(bf(y ~ x) + gaussian(), data = dd)
    st$beta["x"] <- 2
    frm(bf(y ~ x) + gaussian(), dd, start = st)

## The names

Fixed-effect components (`beta`, `betad`) carry the design-matrix column
names, so they agree with
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) and
[`summary()`](https://rdrr.io/r/base/summary.html), and the intercept is
spelled `"(Intercept)"`. Covariance components (`theta`, `thetaac`,
`thetar`) have no design of their own and are spelled `theta_1`,
`theta_2`, ... - the same names `confint(parm =)` reports. When you
supply names BACK, the parentheses may be dropped: `"Intercept"` and
`"(Intercept)"` address the same coefficient.

## Two vocabularies, one model

This is the flat internal parameterization.
[`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
describes the same model in the structured
`class`/`coef`/`group`/`dpar`/`nlpar` vocabulary that
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
addresses, and
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
also accepts a natural-scale spelling (`sigma`, `sd_g__Intercept`). Use
[`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
to write priors and `par_template()` to write starting values and box
constraints.

brms has no counterpart for starting values: `brm(init =)` takes a list
of Stan program names and its own documentation calls it mainly an
internal testing facility. brms does not need one, because it
initializes at random in a bounded range and its priors carry location
information.
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) optimizes
rather than samples, so it evaluates the objective AT the starting
values; this is the discovery step that makes that survivable.

## See also

[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) for
`start =`,
[`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
for the prior addressing vocabulary,
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
for `newparams =`.

## Examples

``` r
dd <- data.frame(x = rnorm(30), g = factor(rep(1:5, 6)))
dd$y <- 1 + 2 * dd$x + rnorm(30)

# before fitting: the names and the cold starting values
par_template(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#> <frmtmb parameter template> starting values
#> $beta  (2)
#> (Intercept)           x 
#>    1.013626    0.000000 
#> $betad  (1)
#> sigma_(Intercept) 
#>         0.9645618 
#> $b  (5)
#> b_1 b_2 b_3 b_4 b_5 
#>   0   0   0   0   0 
#> $theta  (1)
#> theta_1 
#>       0 
#> Edit and pass back as frm(start = ) or frm_simulate(newparams = ). Parentheses are optional in names you supply.

# after fitting: the same layout, holding the estimates
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd)
par_template(fit)
#> <frmtmb parameter template> estimates
#> $beta  (2)
#> (Intercept)           x 
#>   0.9122575   2.4075210 
#> $betad  (1)
#> sigma_(Intercept) 
#>        -0.1233936 
#> $b  (5)
#>           b_1           b_2           b_3           b_4           b_5 
#> -2.860313e-10  1.661331e-09 -4.098257e-10 -1.698764e-09  7.332903e-10 
#> $theta  (1)
#>   theta_1 
#> -10.64264 
#> Edit and pass back as frm(start = ) or frm_simulate(newparams = ). Parentheses are optional in names you supply.

# discover, edit, fit
st <- par_template(bf(y ~ x) + gaussian(), data = dd)
st$beta["x"] <- 2
frm(bf(y ~ x) + gaussian(), dd, start = st)
#> frmtmb fit: y ~ x 
#> Family: gaussian   Method: ML 
#> logLik: -38.8663  AIC: 83.7327  nobs: 30 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      0.9123      2.4075 
#>  sigma:
#> (Intercept) 
#>     -0.1234 
```
