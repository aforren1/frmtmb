# Hypothesis tests on parameter expressions

The frequentist analog of brms's `hypothesis()`: evaluates expressions
of the model parameters at the estimates and tests them against zero. A
hypothesis is `"expr"` (tested against 0) or `"expr = rhs"`, e.g.
`"x1 - x2 = 0"` or `"exp(Intercept) = 1"`.

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
  ...
)

# S3 method for class 'frmtmb_draws'
hypothesis(x, hypothesis, alpha = 0.05, ...)

# S3 method for class 'frmtmb_multiple'
hypothesis(x, hypothesis, alpha = 0.05, ...)
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

## Value

A `frmtmb_hypothesis` object: a data frame with one row per hypothesis
(`estimate`, `se`, `lwr`, `upr`, `z`, `p`) carrying the method payload
in attributes - the bootstrap draws matrix (`attr(., "draws")`) or the
profile curves (`attr(., "profiles")`).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
bootstrap distribution, the profile curve, or the implied Wald normal
density, one panel per hypothesis. Subsetting with `[` drops the
attributes; keep the full object for plotting.

## Details

Available names: the fixed-effect coefficients under their
[`vcov()`](https://rdrr.io/r/stats/vcov.html) row names with parentheses
stripped (`Intercept`, `x`, `sigma_Intercept`, ...), natural-scale
random-effect summaries `sd_<group>__<term>` and
`cor_<group>__<t1>__<t2>` (brms naming), and `sigma` when the residual
SD is a scalar. So an ICC is
`"sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)"`.
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
lists every usable name for a fit.

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
