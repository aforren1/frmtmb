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
```

## Arguments

- x:

  A `frmtmb_fit`.

- ...:

  Backend controls: passed to
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html) for
  `method = "profile"` (e.g. `ytol`, `ystep`, `maxit`, `parm.range`) and
  to [`frm_bootstrap()`](frm_bootstrap.md) for `method = "boot"` (e.g.
  `re.form = NULL` for a conditional bootstrap). Unused for `"wald"` (a
  warning).

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
  [`frm_bootstrap()`](frm_bootstrap.md) (percentile interval; `p` is the
  two-sided percentile p-value, whose resolution is limited by `nsim`;
  `se` is the bootstrap SD). Handles any expression, including the
  variance-component names, whose sampling distributions Wald
  approximates poorly.
