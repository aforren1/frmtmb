# Conditional effects of predictors

For each requested effect, predicts over a grid of that predictor with
every other predictor held at a reference value (numeric: mean; factor:
first level; matrix covariate: column means) and random effects excluded
(`re.form = NA`). Confidence bands are Wald intervals computed on the
link scale and back-transformed. Smooth terms are included, so this also
covers what brms calls `conditional_smooths()`.

## Usage

``` r
conditional_effects(x, ...)

# S3 method for class 'frmtmb_fit'
conditional_effects(
  x,
  effects = NULL,
  resp = NULL,
  dpar = NULL,
  resolution = 100,
  prob = 0.95,
  method = c("epred", "predict"),
  band = c("wald", "profile", "boot"),
  re_formula = NA,
  ndraws = 400,
  boot = NULL,
  profile_points = 25,
  seed = NULL,
  conditions = list(),
  surface = FALSE,
  data = NULL,
  ...
)

# S3 method for class 'frmtmb_draws'
conditional_effects(
  x,
  effects = NULL,
  resp = NULL,
  dpar = NULL,
  resolution = 100,
  prob = 0.95,
  ndraws = NULL,
  re_formula = NA,
  conditions = list(),
  data = NULL,
  ...
)
```

## Arguments

- x:

  A `frmtmb_fit`.

- ...:

  Passed to
  [`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md).

- effects:

  Character vector of variable names, or `"x:z"` pairs; for a pair, the
  first variable is varied over its range while the second is held at
  its levels (factors) or at mean and mean plus or minus one SD
  (numeric). Default: every fixed-effect and smooth variable of the
  selected linear predictor, plus one `"a:b"` pair per fitted
  interaction (brms's default); a term of order three or more
  contributes its leading pair.

- resp, dpar:

  Response and distributional parameter, as in
  [`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md).

- resolution:

  Number of grid points for a varied numeric predictor.

- prob:

  Coverage of the confidence bands (brms spelling).

- method:

  `"epred"` (default): Wald bands for the expected response.
  `"predict"`: prediction intervals - quantile bands from `ndraws`
  responses simulated from the family at each grid point (observation
  noise; random effects stay excluded, as in brms with
  `re_formula = NA`), around the expected response on the same scale as
  the draws (a count under `trials()`, the truncated mean under
  [`trunc()`](https://rdrr.io/r/base/Round.html)). The draws respect the
  response's addition terms: literal
  [`trunc()`](https://rdrr.io/r/base/Round.html) bounds apply, and
  `trials()`, `se()` or variable
  [`trunc()`](https://rdrr.io/r/base/Round.html) bounds must be pinned
  in `conditions` (a grid row is an artificial observation, so a
  reference value for those is meaningless and is an error rather than a
  silent default).

- band:

  How the confidence band is built: `"wald"` (default, the delta method
  on the link scale), `"profile"` (likelihood-root inversion per grid
  point) or `"boot"` (parametric-bootstrap percentiles). See the band
  section. Only for `method = "epred"`.

- re_formula:

  The population switch, in brms's spelling: `NA` (the default) draws
  the population-level curve, `NULL` conditions on the random effects of
  the grid's reference group levels (set them with `conditions =`), and
  a one-sided formula keeps the named terms. The fit surface's
  [`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
  spells the same setting `re.form` after lme4; `conditional_effects()`
  takes brms's name because it is brms's function, and says so if handed
  the other spelling. `band = "profile"` exists only for the
  population-level curve.

- ndraws:

  Simulated responses per grid point for `method = "predict"`. For the
  draws method: how many evenly spaced posterior draws the curves are
  computed over (default all).

- boot:

  For `band = "boot"`: `NULL` (default) runs one
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  of 200 refits and says so, a number runs that many, and a
  `frmtmb_boot` object from an earlier identical call
  (`attr(ce, "boot")`) is reused without refitting anything.

- profile_points:

  For `band = "profile"`: how many points of a numeric grid are
  profiled, the band being interpolated between them on the link scale.

- seed:

  Seed for `band = "boot"`, passed to
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md).

- conditions:

  Named list overriding reference values, e.g. `list(x2 = 1, g = "b")`;
  or a data frame whose rows define multiple condition sets (brms
  style), labeled by a `cond__` column from its row names.

- surface:

  Accepted for brms compatibility. `TRUE` (a fitted surface over two
  predictors) is refused: ask for the two-variable effect `"x1:x2"`
  instead, which varies the first predictor at three values of the
  second.

- data:

  The original model data. Only needed when the model frame does not
  store a raw variable (e.g. a variable used only inside
  [`poly()`](https://rdrr.io/r/stats/poly.html)).

## Value

A named list of data frames (one per effect) with the varied variable(s)
plus `estimate__`, `se__` (link scale), `lower__`, and `upper__`;
printing it draws the plots. An ordinal fit adds a `cats__` column and
one block of rows per response category. `plot(ce, points = TRUE)`
overlays the raw observations (the brms argument), each panel showing
only the observations that belong to its own condition; see the faceting
section. No points are drawn for a per-category ordinal display, a
non-mean `dpar`, or a matrix response (a message says so).

## Several conditions become one faceted page

A `conditions` data frame of several rows gives the effect one panel per
row, laid out as small multiples on a SINGLE page with a shared scale,
the way brms's `facet_wrap("cond__")` does. `plot(ce, ncol = )` sets the
number of columns; the default lays the panels out roughly square, as
brms's `ncol = NULL` does.

The panels are drawn with tinyplot when it is installed, which supplies
the shared axes and a single outer legend. Without it the fallback is a
grid of ordinary base-graphics panels, still one page and still honoring
`ncol`, labeled by condition on the y axis. The per-category ordinal
display always takes the fallback grid: its grouping slot already
carries the response category.

`points = TRUE` draws only the observations belonging to each condition,
matching brms's `make_point_frame()`. A condition claims the rows of the
data that MATCH it on the variables it sets, so with one condition per
level of a factor each observation appears once, in its own panel. Two
rules bound that, both as in brms:

- A NUMERIC condition variable is dropped from the match, because a
  reference value such as `list(x2 = 0.37)` is a point on a continuum
  that names no observation; matching on it would empty every panel. A
  grouping factor stored as a number is exempt, being a label rather
  than a continuum. A condition left with nothing to match on therefore
  claims EVERY observation, and its panel differs from the others only
  in its curve.

- A condition variable that names no column of the model data cannot
  select anything either, so the points are left unsplit and every panel
  draws all of them.

Unlike brms, a variable the `conditions` data frame does not set is not
silently pinned to its first level for this purpose, so an unmentioned
factor does not drop observations from every panel.

## Ordinal responses

[`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
and
[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
have no mean, so the display is per CATEGORY, as brms's
`categorical = TRUE` is: each effect data frame gains a `cats__` factor
of the response's own levels and carries the fitted category probability
in `estimate__`, with one curve per category in the plot (a second
predictor gets a panel of its own).

`se__` is then on the probability scale, and the band is a Wald interval
on the logit of the probability so it cannot leave `[0, 1]`. The
standard errors are the delta method over the joint covariance of the
coefficients, the thresholds AND the `cs()` coefficients: a category
probability depends on all of them, and holding the thresholds fixed
would understate every band.

`method = "predict"` is refused there (the category probabilities are
already the whole predictive distribution). Naming a distributional
parameter, `dpar = "mu"`, opts back into the ordinary display of the
latent linear predictor.

## Confidence bands

`band` picks how `lower__` and `upper__` are found. The estimate is the
same curve in all three cases; only the band changes.

- `"wald"` (default, and free): the delta-method standard error on the
  link scale, back-transformed. Symmetric on the link scale by
  construction.

- `"profile"`: one likelihood-root search
  ([`TMB::tmbroot()`](https://rdrr.io/pkg/TMB/man/tmbroot.html)) per
  grid point. A grid point's linear predictor is a linear combination of
  the coefficients, so the search inverts the likelihood ratio along
  that combination; the link is monotone, so the two endpoints map
  through it. The band is not symmetric and does not assume a quadratic
  log-likelihood, which is what makes it worth its cost near a boundary
  or at a small sample size.

- `"boot"`: percentiles of the grid predictions across the refits of ONE
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  (`re.form = NA`, the same population-level grid). One bootstrap serves
  every effect, condition set and ordinal category of the call;
  `attr(ce, "boot")` returns it, and passing it back as `boot =` costs
  no refits at all.

What `se__` means follows the band: the Wald standard error (link scale,
or the probability scale on the ordinal display) for `"wald"` and
`"profile"` - the profile changes the endpoints, not the standard
error - and the standard deviation of the bootstrap draws, on the
displayed scale, for `"boot"`.

Cost, and how it is capped. A root search is two constrained
optimizations, so `band = "profile"` profiles at most `profile_points`
(default 25) of a numeric grid, spread over its range and including both
ends, and interpolates the endpoints linearly between them on the link
scale. `resolution` still governs the estimate curve. A factor grid is
profiled at every level. Points whose search does not converge become
`NA` and one warning names how many. `band = "boot"` costs its refits
once, however many effects are asked for.

Refusals. `band` other than `"wald"` needs `method = "epred"`: a
prediction interval is already a simulation quantile. `band = "profile"`
additionally needs the displayed quantity to BE a linear combination of
the parameters, so it is refused (naming `"boot"`) for an ordinal
category probability, for an expected response that runs through several
distributional parameters or through truncation bounds (`dpar =` names
one predictor and opts back in), for a nonlinear predictor, for a
predictor carrying `s()`, `gp()` or `hsgp()` (basis coefficients are
inner parameters), and for a REML fit or
`frmtmb_control(profile = TRUE)` (the coefficients are not outer
parameters there).

A nonlinear predictor (`nl = TRUE`) has no delta-method standard error
at all, so `band = "boot"` is its only band and the other two are
refused. The same holds for the per-category display of a nominal
family, whose probabilities have no threshold Jacobian to differentiate.

## Draws objects

On a `frmtmb_draws` object from
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
the same grids are evaluated once per posterior draw, and `estimate__`,
`lower__`, `upper__` and `se__` are the pointwise mean, quantiles and
standard deviation of the drawn curves. There is no `band =` or
`method =` to choose (the band IS the posterior quantile band), a
nonlinear predictor and a nominal per-category display work without a
delta method, and the method runs on formula-route draws that have no
maximum-likelihood fit behind them. `ndraws` thins the draws evenly for
a cheaper curve. Draws from `frm_sample(laplace = TRUE)` are refused:
the sampled vector no longer aligns with the model's parameter template.

## Which predictors are plotted by default

Every variable of the selected linear predictor that the display can
vary: its fixed-effect terms, its smooth terms and its `mo()` terms
(whose design columns are placeholders, so the variable is read from the
term itself). On a nonlinear fit they live one level down - the
covariates the nonlinear formula reads, plus the terms of each nonlinear
parameter's own predictor - and all of those are collected too.
Matrix-valued columns are excluded, a grid over a matrix covariate not
being a curve. Every fitted interaction whose components are single
plottable variables adds an `"a:b"` display, as in brms: a fitted
interaction hidden by default invites reading the main-effect curves as
the whole story. The display takes two variables, so a three-way or
deeper term contributes its leading pair, with the remaining variables
at their reference values until `conditions =` pins them. Naming
`effects =` overrides the search.

## Examples

``` r
set.seed(5)
dd <- data.frame(x = rnorm(120), f = factor(rep(c("a", "b"), 60)))
dd$y <- rnorm(120, 1 + 0.5 * dd$x + (dd$f == "b"), 1)
fit <- frm(bf(y ~ x * f), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = c("x", "x:f"))
plot(ce, ask = FALSE)


# prediction intervals instead of epred bands
ce_p <- conditional_effects(fit, effects = "x", method = "predict")
# \donttest{
# a likelihood-profile band: asymmetric, and no quadratic assumption
ce_pr <- conditional_effects(fit, effects = "x", band = "profile",
                             resolution = 20, profile_points = 5)

# a bootstrap band, reused for a second effect without refitting
ce_b <- conditional_effects(fit, effects = "x", band = "boot",
                            boot = 25, seed = 1)
ce_b2 <- conditional_effects(fit, effects = "x", band = "boot",
                             boot = attr(ce_b, "boot"))
# }
```
