# Set up priors brms-style

Builds prior specifications with brms spelling:
`set_prior("normal(0, 5)", class = "b")`. Combine several with `+` or
[`c()`](https://rdrr.io/r/base/c.html). Distributions: `normal(mu, sd)`,
`student_t(df, mu, sd)`, `cauchy(mu, sd)`, `exponential(rate)`; an empty
string sets bounds only.

## Usage

``` r
set_prior(
  prior = "",
  class = "b",
  coef = "",
  dpar = "",
  group = "",
  lb = NA,
  ub = NA
)
```

## Arguments

- prior:

  Distribution string, e.g. `"normal(0, 5)"`, or a
  [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)/[`prior_t()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  object, or `""` for bounds only.

- class:

  `"b"`, `"Intercept"`, `"sd"`, or `"theta"`.

- coef:

  Restrict to one coefficient (classes `"b"`/`"Intercept"`).

- dpar:

  Distributional parameter (default: the location parameters).

- group:

  Restrict class `"sd"` to one grouping factor.

- lb, ub:

  Optional hard bounds.

## Value

A `frmtmb_priorlist`.

## Details

Classes and their scales:

- `"b"`: population-level coefficients of `dpar` (default: the location
  parameters), excluding the intercept; narrow to one coefficient with
  `coef`. Link scale.

- `"Intercept"`: the intercept of `dpar`. Link scale.

- `"sd"`: random-effect standard deviations (and smoothing SDs), on the
  NATURAL sd scale with the log-Jacobian applied, so
  `set_prior("exponential(1)", class = "sd")` means what it says; narrow
  with `group`. Correlation parameters are not covered (use class
  `"theta"` on the internal scale if you must).

- `"theta"`: raw internal covariance parameters (escape hatch).

When priors overlap, later specifications override earlier ones, so put
class-wide priors first and coefficient-specific ones after. `lb`/`ub`
become hard bounds (via Stan's constrained transforms in
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md));
for class `"sd"` they apply on the sd scale.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), z = rnorm(100),
                 g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# `+` combines specifications; the class-wide one goes first so the
# coefficient-specific one can override it
pr <- set_prior("normal(0, 1)", class = "b") +
  set_prior("normal(0, 0.2)", class = "b", coef = "z") +
  set_prior("exponential(1)", class = "sd", group = "g")
pr
#> normal(0, 1) class=b
#> normal(0, 0.2) class=b coef=z
#> exponential(1) class=sd group=g

# the priors penalize the likelihood: the fit is a MAP estimate
fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, priors = pr)
fixef(fit)$mu
#> (Intercept)           x           z 
#>  1.43922922  0.56382992  0.05323151 
# the tight prior on z shrinks it toward zero
fixef(frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd))$mu
#> (Intercept)           x           z 
#>  1.43901021  0.57150937  0.06947063 

# an empty distribution string sets a hard bound only
set_prior("", class = "b", coef = "x", lb = 0)
#> (bounds only) class=b coef=x lb=0

# get_prior() shows which rows a design offers
get_prior(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#>    prior     class    coef group  dpar resp lb ub
#> 1 (flat) Intercept                          NA NA
#> 2 (flat)         b                          NA NA
#> 3 (flat)         b       x                  NA NA
#> 4 (flat)         b       z                  NA NA
#> 5 (flat) Intercept               sigma      NA NA
#> 6 (flat)        sd                          NA NA
#> 7 (flat)        sd             g            NA NA
#> 8 (flat)     theta                          NA NA
#> 9 (flat)     theta theta_1                  NA NA
```
