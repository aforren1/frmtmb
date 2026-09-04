# Hypothesis tests on parameter expressions

The frequentist analog of brms's `hypothesis()`: evaluates expressions
of the model parameters at the estimates and tests them against zero. A
hypothesis is `"expr"` (tested against 0), `"expr = rhs"`, e.g.
`"x1 - x2 = 0"` or `"exp(Intercept) = 1"`, or brms's directional
`"lhs > rhs"` / `"lhs < rhs"`.

## Usage

``` r
hypothesis(x, ...)

# S3 method for class 'frmtmb_fit'
hypothesis(
  x,
  hypothesis,
  alpha = 0.05,
  method = c("wald", "profile", "boot"),
  nsim = 500,
  seed = NULL,
  class = NULL,
  group = NULL,
  vcov = NULL,
  ...
)

# S3 method for class 'frmtmb_multiple'
hypothesis(x, hypothesis, alpha = 0.05, class = NULL, group = NULL, ...)
```

## Arguments

- x:

  A `frmtmb_fit`, or a `frmtmb_multiple` for pooled tests.

- ...:

  Backend controls: passed to
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html) for
  `method = "profile"` (e.g. `ytol`, `ystep`, `maxit`, `parm.range`) and
  to
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  for `method = "boot"` (e.g. `re.form = NULL` for a conditional
  bootstrap). Unused for `"wald"` (a warning).

- hypothesis:

  Character vector of hypotheses.

- alpha:

  Test level; the reported interval covers `1 - alpha` (brms spelling).

- method:

  `"wald"`, `"profile"`, or `"boot"`.

- nsim:

  Bootstrap draws for `method = "boot"`; all hypotheses share one
  bootstrap run.

- seed:

  Optional seed for `method = "boot"`.

- class, group:

  brms shorthand for the parameter names: the hypothesis is written with
  bare names and `class` (and `group`, for the `sd_`/`cor_` summaries)
  supplies the prefix. The default `NULL` (like brms's `class = "b"`)
  takes the names as written.

- vcov:

  `method = "wald"` only: a covariance matrix over the whole outer
  parameter vector to use in place of the model-based one -
  [`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md)
  with `full = TRUE`, or a function of the fit returning such a matrix.
  The delta-method standard error is then the cluster-robust one, and a
  matrix carrying reference degrees of freedom switches the test to a
  `t` reference.

## Value

A `frmtmb_hypothesis` object: a data frame with one row per hypothesis
(`estimate`, `se`, `lwr`, `upr`, `z`, `p`) carrying the method payload
in attributes - the bootstrap draws matrix (`attr(., "draws")`) or the
profile curves (`attr(., "profiles")`).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
bootstrap distribution, the profile curve, or the implied Wald normal
density, one panel per hypothesis. Subsetting with `[` drops the
attributes; keep the full object for plotting.

## Directional hypotheses

`"lhs > rhs"` and `"lhs < rhs"` test the same difference `(lhs) - (rhs)`
against zero with a one-sided alternative, so the reported `p` is the
one-sided tail probability and the interval is one-sided at level
`1 - alpha`: the unbounded end prints as `Inf` or `-Inf`. `p` is
[`pnorm()`](https://rdrr.io/r/stats/Normal.html) of the signed z
statistic for `"wald"` and, as in the two-sided case where `se`, `z` and
`p` stay Wald-based, for `"profile"` too - the profile changes the
BOUND, which is the matching endpoint of the two-sided `1 - 2 * alpha`
profile interval, and nothing else. For `"boot"` and the draws method
`p` is the tail proportion of the draws with the `(1 + k) / (1 + n)`
correction. Where brms reports the posterior probability of the
direction, this reports its frequentist complement: small `p` is
evidence for the stated direction. `">="` and `"<="` read as `">"` and
`"<"`.

## brms class and group shorthand

`class` and `group` prefix the bare names in the hypothesis, so
`hypothesis(fit, "Intercept - age > 0", class = "sd", group = "patient")`
tests `sd_patient__Intercept - sd_patient__age`. `class = "b"` (brms's
default) and `class = NULL` leave the names alone; a `class` without a
`group` prefixes `<class>_`, which is how a distributional coefficient
such as `sigma_Intercept` is named. A name already written in full keeps
its spelling, so the two can be mixed; but a name that exists only
WITHOUT the prefix is an error rather than a test of the unprefixed
parameter, so a wrong `class` or `group` cannot quietly answer a
different question.

Available names: the fixed-effect coefficients under their
[`vcov()`](https://rdrr.io/r/stats/vcov.html) row names with parentheses
stripped (`Intercept`, `x`, `sigma_Intercept`, ...), natural-scale
random-effect summaries `sd_<group>__<term>` and
`cor_<group>__<t1>__<t2>` (brms naming), and `sigma` when the residual
SD is a scalar. So an ICC is
`"sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)"`.
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
lists every usable name for a fit. The internal spelling that
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) print is accepted as well,
backquoted because it carries parentheses: `` "`(Intercept)` - x" ``.
The traffic runs the other way too: `confint(parm = )` and
`profile(parm = )` take these names, whenever one of them stands for a
single internal parameter.

## Which random-effect blocks contribute names

Every block whose covariance parameters ARE standard deviations and
correlations: the plain structures (`us`, `diag`, `homdiag`, `cs`,
`ar1`, `toep`, the spatial and reduced-rank ones) and the
known-structure blocks `gr(cov = )`, `gr(prec = )` and `equalto()`,
whose `sd_`/`cor_` names describe the WITHIN-level covariance that
multiplies the fixed relationship matrix. That is what makes
heritability-as-ICC writable directly:
`"sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2)"` on an animal
model fitted with `(1 | gr(id, cov = A))`. An `equalto()` block
estimates nothing, so its names are constants with zero variance.

An `|ID|`-merged block is ONE block, so it contributes one name per
merged coefficient, and the names carry the linear predictor they came
from just as correlated slopes do. A two-trait animal model written
`(1 | q | gr(id, cov = A))` in both formulas of an
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) gives
`sd_id__y1.muIntercept`, `sd_id__y2.muIntercept` and
`cor_id__y1.muIntercept__y2.muIntercept` - the last being the genetic
correlation between the traits.
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
prints them.

Two blocks on the same grouping factor with the same term name - an
animal model's `(1 | gr(id, cov = A)) + (1 | id)`, where the genetic and
permanent-environment terms both name the group `id` - collide on one
`sd_id__Intercept`, and the first block in formula order claims it. Give
the second term its own grouping column (a copy of the factor under
another name) when both are wanted by name.

Excluded: `s()`/`t2()` smooths, `gp()`/`hsgp()`, `car()` and `spde()`.
Their theta segments are not standard deviations - an inverse smoothing
parameter, lengthscales, a mixing proportion, a precision and an inverse
range - so there is no `sd_<group>__<term>` to name. Read those off
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md),
which reports each under its own label (`sd(gp)`, `range(gp)`,
`sd(car)`, ...).

## When a coefficient shadows a natural-scale name

Model terms and natural-scale summaries share one namespace here, and
the coefficient wins: a covariate literally named `sigma` makes
`"sigma = 0"` a test on ITS coefficient, not on the residual standard
deviation. The same holds for a coefficient that spells out
`sd_<group>__<term>`, `cor_...` or an autocorrelation name such as
`ar1`. The shadowed quantity keeps a name: prefix it with a dot,
`.sigma`, `.sd_g__Intercept`, `.ar1`. The dot spelling exists only where
a collision does, and `hypothesis()` says so once per call when one is
in play.
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
lists both names in that case.

## See also

[`vcov.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/vcov.frmtmb_fit.md)
with `full = TRUE` for the same joint covariance (fixed effects plus
covariance parameters, on their internal scale) as a matrix, which is
what the `"wald"` method uses here.

Methods:

- `"wald"` (default): delta-method z-test, finite-difference gradient
  against the joint parameter covariance (under REML, from the joint
  precision).

- `"profile"`: profile-likelihood interval via
  [`TMB::tmbroot()`](https://rdrr.io/pkg/TMB/man/tmbroot.html) with a
  `lincomb` direction. Only for hypotheses that are linear in the
  parameters, and only for ML fits; `se`, `z`, and `p` stay Wald-based -
  the method changes the interval.

- `"boot"`: parametric bootstrap through
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  (percentile interval; `p` is the two-sided percentile p-value, whose
  resolution is limited by `nsim`; `se` is the bootstrap SD). Handles
  any expression, including the variance-component names, whose sampling
  distributions Wald approximates poorly.

For a
[`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
result the Wald estimate and delta-method variance are computed per
imputation and pooled by Rubin's rules with Barnard-Rubin degrees of
freedom; the returned table carries `t` and `df` columns in place of `z`
(reference t distribution, not normal), and only Wald inference is
available.

## Examples

``` r
set.seed(4)
dd <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
                 g = factor(rep(1:10, 12)))
dd$y <- rnorm(120, 1 + 0.6 * dd$x1 + 0.4 * dd$x2 +
                rnorm(10, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x1 + x2 + (1 | g)) + gaussian(), data = dd)
hypothesis(fit, c("x1 - x2 = 0", "exp(Intercept)"))
#> Hypothesis tests (method = wald)
#>      hypothesis estimate     se   lwr    upr     z         p
#>     x1 - x2 = 0    0.370 0.1280 0.119 0.6209 2.889 3.859e-03
#>  exp(Intercept)    2.507 0.4829 1.561 3.4540 5.192 2.079e-07
# brms's directional form: one-sided p, one-sided interval
hypothesis(fit, "x1 > x2")
#> Hypothesis tests (method = wald)
#>  hypothesis estimate    se    lwr upr     z       p
#>     x1 > x2     0.37 0.128 0.1594 Inf 2.889 0.00193
#>   rows written with '<' or '>' are one-sided: p and the finite interval
#>   bound hold at level 0.95
# class/group name the natural-scale random-effect summaries
hypothesis(fit, "Intercept > 0", class = "sd", group = "g")
#> Hypothesis tests (method = wald)
#>     hypothesis estimate    se    lwr upr     z         p
#>  Intercept > 0   0.5494 0.152 0.2993 Inf 3.614 0.0001508
#>   rows written with '<' or '>' are one-sided: p and the finite interval
#>   bound hold at level 0.95
# variance-component expressions: an ICC with bootstrap intervals
hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)",
           method = "boot", nsim = 20, seed = 1)
#> Hypothesis tests (method = boot)
#>   bootstrap draws: 20 (0 failed or not converged)
#>                                         hypothesis estimate     se     lwr
#>  sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)   0.2694 0.1231 0.03343
#>     upr     z       p
#>  0.4187 2.189 0.09524
```
