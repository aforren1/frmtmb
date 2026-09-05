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

  Optional hard bounds, on the scale described in Hard bounds.

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

- `"ar"`, `"ma"`, `"cosy"`, `"cortime"`: the R-side residual correlation
  of an [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`,
  `cosy()` or `unstr()` term, under brms's own class names. `"ar"`,
  `"ma"` and `"cosy"` take an ordinary density on the NATURAL
  coefficient with that map's Jacobian applied, as class `"sd"` does;
  `"cortime"` takes `lkj(eta)` on an `unstr()` time correlation as a
  whole. Narrow to one response with `resp`. See Residual correlation
  below.

- `"rescor"`: the residual correlation BETWEEN responses of a
  multivariate model (`set_rescor(TRUE)`), as a whole. `lkj(eta)` only,
  as brms spells it.

- `"theta"`: raw internal covariance parameters (escape hatch). `coef`
  names one by its internal name and spans all three covariance
  components: `"theta_2"` for a random-effect block, `"thetaac_1"` for a
  residual autocorrelation, `"thetar_1"` for a residual correlation.
  This is the one spelling that reaches a single parameter of a
  structure whose natural coefficients are not free of one another.

When priors overlap, later specifications override earlier ones, so put
class-wide priors first and coefficient-specific ones after. A class
`"theta"` prior on a position an earlier `"cor"` prior covers replaces
that whole LKJ term, and the other way round, so "later wins" holds
between the two spellings as well. `lb`/`ub` become hard bounds. See
Hard bounds.

## Hard bounds

`lb`/`ub` are how a box constraint is written. A specification may carry
bounds alone (`prior = ""`), a distribution alone, or both, and a later
bounds-only specification tightens an entry an earlier distribution
created rather than replacing it.

A bound is addressed exactly like the distribution beside it, so
`set_prior("", nlpar = "guess", lb = 0, ub = 1)` bounds the nonlinear
parameter `guess`, and `dpar`, `resp`, `group` and `coef` narrow a bound
the same way they narrow a density. As with a distribution, class `"b"`
with `nlpar` covers EVERY coefficient of that parameter; `coef` picks
out one.

The scale is the parameter's own: class `"sd"` bounds a standard
deviation on the sd scale (frmtmb stores its log), classes `"ar"`,
`"ma"` and `"cosy"` bound the natural coefficient of a residual
structure, and class `"theta"` bounds an internal covariance parameter
itself, one position at a time with `coef = "theta_2"`, `"thetaac_1"` or
`"thetar_1"`. Everything else is bounded on the internal (link) scale,
so a bound on a log-linked dispersion is a bound on its logarithm.

In [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) a bound
is a box constraint handed to the optimizer; in
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
it becomes one of Stan's constrained transforms. Both take this spelling
and no other: the `lower`/`upper` arguments of releases before 0.49 are
gone rather than aliased, and a call still using them fails as an unused
argument. Every outer parameter they could reach has a class here, down
to a single internal covariance parameter.

Where the two spellings differed, this one broadcasts: a bound carried
by `nlpar =` covers every coefficient of that parameter, the way a prior
does, and `coef` narrows it to one. When two specifications bound the
same parameter the later one wins, so a bounds-only specification after
a wide one tightens it.

## Residual correlation

frmtmb holds an [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`,
`arma()`, `cosy()` or `unstr()` residual block in one unconstrained
vector, chosen so the optimizer cannot step outside the stationary and
invertible region. A prior is still written about the parameter brms
names, and carried onto that vector with the log Jacobian of the map,
exactly as class `"sd"` carries a density on a standard deviation onto
its logarithm. So `set_prior("normal(0, 0.5)", class = "ar")` is a
density on the AR coefficient itself, and
[`summary()`](https://rdrr.io/r/base/summary.html) reports the parameter
the prior was written about.

Bounds behave the same way where the map allows it. A first-order `ar`,
`ma` or `cosy` coefficient is a monotone function of one internal
parameter, so `lb`/`ub` map exactly onto a box. At order two and above
they do not: `ar[1]` is a function of every internal parameter of the
block at once, so no box in internal space is the box asked for, and
`lb`/`ub` are refused rather than approximated. Little is lost, because
the parameterization already guarantees stationarity and invertibility,
which is what such a bound is usually for; where a hard box really is
wanted, `class = "theta"` with `coef = "thetaac_1"` bounds one internal
parameter.

`cosy` is bounded below at `-1/(d - 1)` for `d` time points, where a
compound-symmetric matrix stops being positive definite, and a bound
outside that window is refused rather than clamped. brms bounds `cosy`
on `[0, 1]` instead, so a negative estimate here has no brms
counterpart.

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

A prior with a location places
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)'s `start`
for a nonlinear parameter. `normal()`, `student_t()` and `cauchy()` all
carry one, and where `start` does not set a nonlinear coefficient, that
coefficient begins at the prior's location, reported in a message. Other
parameters keep their usual starts: a prior is a penalty, not a claim
about where to begin. Without a located prior a nonlinear model still
needs `start`, because
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) evaluates
the objective AT the starting values;
[`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md)
names them.

`resp` picks one response of a multivariate model; the default priors of
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
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

# bounds address a nonlinear parameter the way a distribution does,
# so a guessing rate is held in [0, 1]
set_prior("", nlpar = "guess", lb = 0, ub = 1)
#> (bounds only) class=b nlpar=guess lb=0 ub=1

# the residual-correlation classes are brms's own names
set_prior("normal(0, 0.5)", class = "ar")
#> normal(0, 0.5) class=ar
set_prior("lkj(2)", class = "rescor")
#> lkj(2) class=rescor
# and one internal covariance parameter, by its name
set_prior("", class = "theta", coef = "thetaac_1", lb = -2, ub = 2)
#> (bounds only) class=theta coef=thetaac_1 lb=-2 ub=2

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
#> route = "fit": the prior defaults frm() applies
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
