# Hypothesis tests on parameter expressions

The frequentist analog of brms's `hypothesis()`: evaluates expressions
of the model parameters at the estimates and tests them against zero. A
hypothesis is `"lhs = rhs"`, e.g. `"x1 - x2 = 0"` or
`"exp(Intercept) = 1"`, or brms's directional `"lhs > rhs"` /
`"lhs < rhs"`. Every hypothesis states a relation: a string with no `=`,
`<` or `>` is refused, as it is in brms, so write `"x1 = 0"` rather than
`"x1"`.

## Usage

``` r
hypothesis(x, ...)

# S3 method for class 'frmtmb_fit'
hypothesis(
  x,
  hypothesis,
  class = "b",
  group = "",
  scope = c("standard", "ranef", "coef"),
  alpha = 0.05,
  robust = FALSE,
  seed = NULL,
  method = c("wald", "profile", "boot"),
  nsim = 500,
  vcov = NULL,
  ...
)

# S3 method for class 'frmtmb_multiple'
hypothesis(
  x,
  hypothesis,
  class = "b",
  group = "",
  scope = c("standard", "ranef", "coef"),
  alpha = 0.05,
  robust = FALSE,
  seed = NULL,
  ...
)
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
  for `method = "boot"` (e.g. `re_formula = NULL` for a conditional
  bootstrap). Refused for `"wald"`.

- hypothesis:

  Character vector of hypotheses. Names on it become the `Hypothesis`
  labels, as in brms.

- class, group:

  brms's name prefix; see *brms class and group shorthand*. The default
  `class = "b"` is brms's.

- scope:

  brms's `"standard"` is the only scope a fit supports: `"ranef"` and
  `"coef"` evaluate the hypothesis on each group level's draws, which a
  maximum-likelihood fit does not have, and are refused by name.

- alpha:

  Test level; the reported interval covers `1 - alpha` for a two-sided
  row and `1 - 2 * alpha` for a directional one, as in brms.

- robust:

  brms's median-and-MAD switch. It summarizes draws and a fit has none,
  so `TRUE` is refused by name.

- seed:

  Optional seed for `method = "boot"`.

- method:

  `"wald"`, `"profile"`, or `"boot"`.

- nsim:

  Bootstrap draws for `method = "boot"`; all hypotheses share one
  bootstrap run.

- vcov:

  `method = "wald"` only: a covariance matrix over the whole outer
  parameter vector to use in place of the model-based one -
  [`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md)
  with `full = TRUE`, or a function of the fit returning such a matrix.
  The delta-method standard error is then the cluster-robust one, and a
  matrix carrying reference degrees of freedom switches the test to a
  `t` reference.

## Value

A `frmtmb_hypothesis` list in brms's shape; see *The returned object*.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) shows the
bootstrap distribution, the profile curve, or the implied Wald normal
density, one panel per hypothesis.

## The returned object

brms's SHAPE under frmtmb's own class: a list of class
`"frmtmb_hypothesis"` with the elements brms has, in brms's order. It
does not carry brms's `brmshypothesis` class. frmtmb owns
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for its own
class and a frmtmb fit is not a brms fit, so `is(x, "brmshypothesis")`
in a ported script is a rule-2 divergence like the other seventeen the
port ledger records.

- `hypothesis`: a data frame with brms's eight columns, one row per
  hypothesis. On a maximum-likelihood fit they mean:

  - `Hypothesis`: brms's label, `(lhs)-(rhs) > 0`, or the name given on
    the hypothesis vector.

  - `Estimate`: the expression at the estimates.

  - `Est.Error`: its delta-method standard error (`"wald"`,
    `"profile"`), the bootstrap standard deviation (`"boot"`), or the
    Rubin pooled standard error (a
    [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
    result).

  - `CI.Lower`, `CI.Upper`: brms's interval, which is central at
    `1 - alpha` for `"="` and central at `1 - 2 * alpha` for a
    directional row, so that its relevant end is the one-sided bound.
    Wald, profile-likelihood or bootstrap percentile, per `method`.

  - `Evid.Ratio`, `Post.Prob`: `NA`. They are posterior quantities, and
    a fit has no posterior.

  - `Star`: `"*"` when a two-sided row's interval excludes 0, or when a
    directional row's one-sided test rejects at level `alpha` (brms
    stars a posterior probability above `1 - alpha` there).

- `samples`: brms's frame of draws, columns `H1`, `H2`, ...: the
  bootstrap replicates for `"boot"`, and no rows otherwise.

- `prior_samples`: the same columns, all `NA`.

- `class`: the name prefix `class` and `group` produced, without its
  trailing underscores, as brms stores it.

- `alpha`.

What brms's frame has no column for rides on attributes, so the frame
keeps brms's shape: `attr(h, "test")` is the test statistic and its
p-value per row (`z`, or `t` with `df` for a pooled or
degrees-of-freedom-carrying covariance), `attr(h, "method")`, and the
method payload, `attr(h, "draws")` for the bootstrap matrix and
`attr(h, "profiles")` for the profile curves.

## Directional hypotheses

`"lhs > rhs"` and `"lhs < rhs"` test the same difference `(lhs) - (rhs)`
against zero with a one-sided alternative, so the reported `p` is the
one-sided tail probability. `p` is
[`pnorm()`](https://rdrr.io/r/stats/Normal.html) of the signed z
statistic for `"wald"` and, as in the two-sided case where the standard
error and the statistic stay Wald-based, for `"profile"` too: the
profile changes the INTERVAL and nothing else. For `"boot"` `p` is the
tail proportion of the replicates with the `(1 + k) / (1 + n)`
correction. Where brms reports the posterior probability of the
direction, this reports its frequentist complement: small `p` is
evidence for the stated direction. `">="` and `"<="` read as `">"` and
`"<"`.

## brms class and group shorthand

The hypothesis is read as `brms:::eval_hypothesis()` reads it. Every
variable in the text gets the prefix `class` and `group` make, and a
prefixed name the model does not have is refused with brms's message,
"Some parameters cannot be found in the model". brms's default
`class = "b"` reads `"x1 - x2 = 0"` as `b_x1 - b_x2`, so it refuses
`"b_x1 = 0"`, which it reads as `b_b_x1`.
`class = "sd", group = "patient"` reads `"Intercept - age > 0"` as
`sd_patient__Intercept - sd_patient__age`. `class = NULL` (or `""`)
takes every name as written, which is what `sigma` or `sd_g__Intercept`
needs, as in brms.

A name may carry what brms's names carry: `x:fe` is the interaction
coefficient `b_x:fe`, not R's `:` operator, and `ar[1]` is the
autocorrelation. brms's renaming (`:` to `___`, `[` and `]` to `.`, `,`
to `..`) is applied to the text before it is parsed, as in brms.

Available names, which
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
lists, are brms's: the coefficients (`b_Intercept`, `b_x`, `b_IxE2` for
`I(x^2)`, `b_sigma_Intercept` for a `sigma` formula, `bs_sx_1` for the
unpenalized part of `s(x)`); a distributional parameter nobody wrote a
formula for, on its natural scale (`sigma`, `shape`, `nu`, and
`sigma_ya` for response `y_a` of a multivariate model); the group-level
summaries `sd_<group>__<coef>` and `cor_<group>__<coef1>__<coef2>`,
where the coefficient of a distributional or nonlinear parameter carries
that parameter's name (`sd_g__sigma_Intercept`); the autocorrelation
parameters (`ar[1]`, `cosy`); and the residual correlations
`rescor__<resp1>__<resp2>`. So an ICC is
`hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0", class = NULL)`.
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
`"sd_id__Intercept^2 / (sd_id__Intercept^2 + sigma^2) = 0"` with
`class = NULL` on an animal model fitted with `(1 | gr(id, cov = A))`.
An `equalto()` block estimates nothing, so its names are constants with
zero variance.

An `|ID|`-merged block is ONE block, so it contributes one name per
merged coefficient, and the names carry the predictor they came from in
brms's spelling. A two-trait animal model written
`(1 | q | gr(id, cov = A))` in both formulas of an
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) gives
`sd_id__y1_Intercept`, `sd_id__y2_Intercept` and
`cor_id__y1_Intercept__y2_Intercept`, the last being the genetic
correlation between the traits.
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
prints them.

Two terms that give one grouping factor the same coefficient, an animal
model's `(1 | gr(id, cov = A)) + (1 | id)`, are refused when the model
is built, with brms's message "Duplicated group-level effects are not
allowed". Give the second term a copy of the factor under another name,
`(1 | gr(id, cov = A)) + (1 | id_pe)`, and the two are
`sd_id__Intercept` and `sd_id_pe__Intercept`.

Excluded: `s()`/`t2()` smooths, `gp()`/`hsgp()`, `car()` and `spde()`.
Their theta segments are not standard deviations: an inverse smoothing
parameter, lengthscales, a mixing proportion, a precision and an inverse
range. There is no `sd_<group>__<coef>` to name. Read those off
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md),
which reports each under its own label (`sd(gp)`, `range(gp)`,
`sd(car)`, ...).

## Names that would collide

brms's renaming can give two parameters one name, and frmtmb does what
brms does with each case:

- Within one predictor, a covariate whose renamed column repeats
  another's, `y ~ Intercept + x` (`(Intercept)` and `Intercept` are both
  `Intercept`), is refused with brms's "Internal renaming led to
  duplicated names".

- Across predictors, the later name takes brms's `__1` suffix: in
  `bf(y ~ sigma_z, sigma ~ z)` the covariate `sigma_z` of `mu` and the
  coefficient `z` of `sigma` are `b_sigma_z` and `b_sigma_z__1`.

- A group-level coefficient given twice on one group is refused, see
  above.

- Two responses brms spells alike, `y_a` and `ya`, are refused with
  brms's "Cannot use the same response variable twice".

- A group-level label given twice on draws, from levels `lvl 1` and
  `lvl.1`, takes brms's `__1` on the later level.

- An interaction group two of whose levels brms joins to one string
  (`1_2:3` and `1:2_3`) is refused: brms pools the two levels.

A coefficient and a natural-scale name cannot meet: every coefficient
starts `b_`, `bs_` or `bsp_`, and `sigma` does not.

A mixture's weights with no theta formula are brms's simplex,
`theta1 ... thetaK`, the mixing probabilities, computed together from
the `K - 1` estimated log ratios against the last component.

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
  parameters, and only for ML fits; the standard error and the test stay
  Wald-based, because the method changes the interval.

- `"boot"`: parametric bootstrap through
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  (percentile interval; `p` is the two-sided percentile p-value, whose
  resolution is limited by `nsim`; `Est.Error` is the bootstrap SD).
  Handles any expression, including the variance-component names, whose
  sampling distributions Wald approximates poorly.

For a
[`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
result the Wald estimate and delta-method variance are computed per
imputation and pooled by Rubin's rules with Barnard-Rubin degrees of
freedom; the test attribute carries `t` and `df` in place of `z`, and
only Wald inference is available.

## Examples

``` r
set.seed(4)
dd <- data.frame(x1 = rnorm(120), x2 = rnorm(120),
                 g = factor(rep(1:10, 12)))
dd$y <- rnorm(120, 1 + 0.6 * dd$x1 + 0.4 * dd$x2 +
                rnorm(10, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x1 + x2 + (1 | g)) + gaussian(), data = dd)
h <- hypothesis(fit, c("x1 - x2 = 0", "exp(Intercept) = 1"))
h
#> Hypothesis Tests for class b:
#>                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
#> 1              (x1-x2) = 0     0.37      0.13     0.12     0.62         NA
#> 2 (exp(Intercept))-(1) = 0     1.51      0.48     0.56     2.45         NA
#>   Post.Prob Star
#> 1        NA    *
#> 2        NA    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> Method: wald. Est.Error is the standard error; Evid.Ratio and
#> Post.Prob are NA, because a maximum-likelihood fit has no posterior.
#> '*': For one-sided hypotheses, the one-sided test rejects at level 0.05;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#>                Hypothesis     z        p
#>               (x1-x2) = 0 2.889 0.003859
#>  (exp(Intercept))-(1) = 0 3.122 0.001799
h$hypothesis$Est.Error
#> [1] 0.1280463 0.4829380
attr(h, "test")$p
#> [1] 0.003859403 0.001799292
# brms's directional form: one-sided p, and a 90% interval
hypothesis(fit, "x1 > x2")
#> Hypothesis Tests for class b:
#>      Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#> 1 (x1)-(x2) > 0     0.37      0.13     0.16     0.58         NA        NA    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> Method: wald. Est.Error is the standard error; Evid.Ratio and
#> Post.Prob are NA, because a maximum-likelihood fit has no posterior.
#> '*': For one-sided hypotheses, the one-sided test rejects at level 0.05;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#>     Hypothesis     z       p
#>  (x1)-(x2) > 0 2.889 0.00193
# class/group name the natural-scale random-effect summaries
hypothesis(fit, "Intercept > 0", class = "sd", group = "g")
#> Hypothesis Tests for class sd_g:
#>        Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
#> 1 (Intercept) > 0     0.55      0.15      0.3      0.8         NA        NA
#>   Star
#> 1    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> Method: wald. Est.Error is the standard error; Evid.Ratio and
#> Post.Prob are NA, because a maximum-likelihood fit has no posterior.
#> '*': For one-sided hypotheses, the one-sided test rejects at level 0.05;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#>       Hypothesis     z         p
#>  (Intercept) > 0 3.614 0.0001508
# variance-component expressions need class = NULL, as in brms: an
# ICC with bootstrap intervals
hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2) = 0",
           class = NULL, method = "boot", nsim = 20, seed = 1)
#> Hypothesis Tests for class :
#>                                            Hypothesis Estimate Est.Error
#> 1 (sd_g__Intercept^2/(sd_g__Intercept^2+sigma^2)) = 0     0.27      0.12
#>   CI.Lower CI.Upper Evid.Ratio Post.Prob Star
#> 1     0.03     0.42         NA        NA    *
#> ---
#> 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
#> Method: boot (20 bootstrap draws, 0 failed or not converged). Est.Error is the bootstrap SD; Evid.Ratio and
#> Post.Prob are NA, because a maximum-likelihood fit has no posterior.
#> '*': For one-sided hypotheses, the one-sided test rejects at level 0.05;
#> for two-sided hypotheses, the value tested against lies outside the 95%-CI.
#>                                           Hypothesis     z       p
#>  (sd_g__Intercept^2/(sd_g__Intercept^2+sigma^2)) = 0 2.189 0.09524
```
