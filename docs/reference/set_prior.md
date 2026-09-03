# Set up priors brms-style

Builds prior specifications with brms spelling:
`set_prior("normal(0, 5)", class = "b")`. Combine several with `+` or
[`c()`](https://rdrr.io/r/base/c.html). Distributions: `normal(mu, sd)`,
`student_t(df, mu, sd)`, `cauchy(mu, sd)`, `exponential(rate)`,
`lkj(eta)`; an empty string sets bounds only.

## Usage

``` r
set_prior(
  prior = "",
  class = "b",
  coef = "",
  group = "",
  resp = "",
  dpar = "",
  nlpar = "",
  lb = NA,
  ub = NA
)
```

## Arguments

- prior:

  Distribution string, e.g. `"normal(0, 5)"`, or a
  [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)/[`prior_t()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)/[`prior_lkj()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  object, or `""` for bounds only.

- class:

  `"b"`, `"Intercept"`, `"sd"`, `"cor"`, or `"theta"`.

- coef:

  Restrict to one coefficient (classes `"b"`/`"Intercept"`).

- group:

  Restrict class `"sd"` or `"cor"` to one grouping factor.

- resp:

  Response of a multivariate model.

- dpar:

  Distributional parameter (default: the location parameters).

- nlpar:

  Nonlinear parameter of an `nl = TRUE` formula. See Nonlinear
  parameters.

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
  with `group`.

- `"cor"`: the CORRELATION of a random-effect block, as a whole.
  `lkj(eta)` only, and it addresses a BLOCK the way class `"sd"` does,
  by `group`; `set_prior("lkj(2)", class = "cor")` covers every
  correlated block of the model, which is brms's spelling. See The LKJ
  prior below.

- `"theta"`: raw internal covariance parameters (escape hatch).

When priors overlap, later specifications override earlier ones, so put
class-wide priors first and coefficient-specific ones after. A class
`"theta"` prior on a position an earlier `"cor"` prior covers replaces
that whole LKJ term, and the other way round, so "later wins" holds
between the two spellings as well. `lb`/`ub` become hard bounds (via
Stan's constrained transforms in
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md));
for class `"sd"` they apply on the sd scale.

## The LKJ prior

`lkj(eta)` is the density `det(C)^(eta - 1)` over a block's correlation
matrix `C`, normalized: `eta = 1` is uniform over correlation matrices,
larger `eta` concentrates toward the identity. frmtmb holds a
correlation as an unconstrained row-normalized Cholesky parameter rather
than as `C`, so the density is carried onto those parameters with the
exact Jacobian of that map (the derivation is in the source of
`R/priors.R`; `tests/testthat/test-lkj.R` checks the sampled
correlations against the closed-form LKJ marginals). The prior a FLAT
correlation parameter carries instead is `(1 - rho^2)^(-3/2)`, which is
improper.

It fits `us()` and `gr(cov = )` blocks of two or more terms, which hold
a whole correlation matrix, and the one-parameter structures `cs()`,
`ar1()` and `hetar1()`, whose single bounded correlation takes the LKJ
marginal `(1 - rho^2)^(eta - 1)` with that structure's own Jacobian. A
`cs()` correlation is bounded below at `-1/(d - 1)`, where a
compound-symmetric matrix stops being positive definite, and the density
is renormalized over that window. `toep()` is refused: its
parameterization is not positive definite everywhere, so it has no
correlation matrix to put a density on.

## Nonlinear parameters

`nlpar` addresses one parameter of an `nl = TRUE` formula, brms's
spelling: `set_prior("normal(5000, 1000)", nlpar = "ult")`, or
`prior(normal(5000, 1000), nlpar = "ult")`. Class `"b"` there covers
EVERY coefficient of that parameter, its intercept included, because a
nonlinear parameter's sub-formula is not centered and brms holds its
intercept in the same coefficient vector as its slopes. That is why the
vignette spelling above lands on `ult_(Intercept)` rather than on
nothing. Narrow to one column with `coef` (`"Intercept"` and
`"(Intercept)"` both name the intercept), or write
`class = "Intercept", nlpar = "ult"`, which is frmtmb's spelling of the
same slot. `nlpar` narrows classes `"sd"` and `"cor"` to the
random-effect blocks of that parameter as well.

An identification prior does NOT stand in for
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)'s `start`.
brms places its sampler with the priors;
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) optimizes,
and it evaluates the objective at the starting values before any penalty
can steer it, which for a nonlinear body usually means an undefined
likelihood at zero. The prior means read across as starting values, and
the refusal names `start` when they are missing.

`resp` picks one response of a multivariate model; the default priors of
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
still stay off there (see its Default priors section), so a multivariate
model's priors are the ones written by hand.

brms's `tag` and `check` have no counterpart: `tag` names a prior for
reuse inside a Stan program, and `check` passes an unchecked string
through to one. frmtmb compiles no Stan program, so both are omitted
rather than accepted and ignored.

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
fit <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, prior = pr)
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

# class "cor" addresses a correlated block as a whole, brms's
# spelling; eta > 1 pulls the correlation toward zero
dd$z <- rnorm(100)
dd$y2 <- dd$y + rnorm(10, 0, 0.6)[dd$g] * dd$z
fitc <- frm(bf(y2 ~ x + z + (z | g)) + gaussian(), data = dd,
            prior = set_prior("lkj(4)", class = "cor"))
VarCorr(fitc)
#>   z | g 
#>         Name Std.Dev. (Intercept)
#>  (Intercept)  0.84786            
#>            z  0.76010       0.334

# get_prior() shows which rows a design offers
get_prior(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
#>    prior     class    coef group  dpar nlpar resp lb ub
#> 1 (flat) Intercept                                NA NA
#> 2 (flat)         b                                NA NA
#> 3 (flat)         b       x                        NA NA
#> 4 (flat)         b       z                        NA NA
#> 5 (flat) Intercept               sigma            NA NA
#> 6 (flat)        sd                                NA NA
#> 7 (flat)        sd             g                  NA NA
#> 8 (flat)     theta                                NA NA
#> 9 (flat)     theta theta_1                        NA NA

# prior() quotes its first argument, brms's spelling, and reaches
# the same machinery
prior(normal(0, 1), class = "b")
#> normal(0, 1) class=b
```
