# Fit a model

Fits a model specified with
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) by maximum
likelihood, using the Laplace approximation for random effects through
RTMB.

## Usage

``` r
frm(
  formula,
  data,
  family = NULL,
  REML = FALSE,
  start = NULL,
  control = frmtmb_control(),
  se = FALSE,
  na.action = stats::na.omit,
  lower = NULL,
  upper = NULL,
  priors = NULL,
  quadrature = FALSE,
  data2 = list(),
  dry_run = NULL,
  verbose = FALSE
)
```

## Arguments

- formula:

  A `frmtmb_formula` from
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) (with a
  family attached via `+`), or a plain formula combined with the
  `family` argument. With neither, the family is
  [`gaussian()`](https://rdrr.io/r/stats/family.html).

- data:

  A data frame. A `tibble`, a `data.table`, or a plain named list of
  equal-length columns is accepted as well, since each reaches
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  unchanged. A matrix column is a supported model variable and enters
  the design as its own block of columns (this is how a functional
  predictor or a [`cbind()`](https://rdrr.io/r/base/cbind.html) term is
  written). A list column is not a model variable and is refused by
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) if the
  formula names one, though it may sit unused in `data`.

- family:

  A family: a `frmtmb_family`, a
  [stats::family](https://rdrr.io/r/stats/family.html) object, a family
  constructor with or without its parentheses (for example
  [`gaussian()`](https://rdrr.io/r/stats/family.html), `poisson`,
  `cumulative`), or a family name as a string. It overrides a family
  already attached to `formula`; in a multivariate model it fills only
  the responses that have none. The default, `NULL`, means
  [`gaussian()`](https://rdrr.io/r/stats/family.html) - the brms, `lme4`
  and `glmmTMB` convention - so `frm(y ~ x, data = d)` is a linear
  model.

- REML:

  If `TRUE`, integrate the `mu` fixed effects out of the likelihood
  along with the random effects (restricted maximum likelihood). For a
  gaussian model the Laplace approximation is exact for that integral,
  so this IS classical Patterson-Thompson REML. For any other family no
  error-contrast derivation exists; what is computed is the
  Laplace-approximated integrated likelihood for the variance parameters
  under a flat prior on the fixed effects, which agrees to first order
  with the Cox-Reid adjusted profile likelihood (the general-model REML
  analogue). glmmTMB computes the same quantity. It usually reduces the
  finite-sample downward bias of variance components, which is REML's
  purpose, but it is an approximation stacked on an approximation, not
  an exactness result. The usual restriction carries over:
  [`anova()`](https://rdrr.io/r/stats/anova.html) compares REML fits
  only when their fixed-effect designs agree.

  Distributional-parameter coefficients (a `sigma ~ z` model's, for
  example) are deliberately NOT integrated: they belong to the
  variance-parameter set, exactly as a `varExp()` coefficient does in
  `nlme::gls(method = "REML")`, and the two agree to optimizer precision
  for gaussian models. Random effects appearing in a distributional
  parameter's own formula are integrated like any other latent.

- start:

  Optional named list of starting values; components must match the
  parameter template (`beta`, `betad`, `theta`).

- control:

  A list from
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md).

- se:

  If `TRUE`, run `RTMB::sdreport()` at fit time. The default (`FALSE`)
  defers it until standard errors are first needed (`summary`, `vcov`,
  `confint`, `predict(se.fit = TRUE)`), which cuts roughly a quarter off
  fit time in fit-and-predict or bootstrap loops. The deferred report is
  cached, so nothing is computed twice.

- na.action:

  How to handle missing values, as in
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) (default
  [stats::na.omit](https://rdrr.io/r/stats/na.fail.html)). Rows dropped
  for missingness are reported in a message; wrap the call in
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) to silence
  it.

- lower, upper:

  Optional named numeric vectors of hard box constraints on outer
  parameters (brms `lb`/`ub`), on the internal scale, e.g.
  `lower = c(b = 0)` for a nonlinear rate parameter. Names as in
  [`confint()`](https://rdrr.io/r/stats/confint.html) rows, with
  parentheses optional; a nonlinear parameter declared intercept-only
  (`b ~ 1`) may be named bare, as in that example, which resolves to
  `b_(Intercept)`. One that carries several coefficients is refused
  rather than resolved to one of them.

- priors:

  Optional
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  specification. This makes the fit MAP / regularized ML (glmmTMB's
  `priors=` in spirit): useful for stabilizing singular variance
  components or separating binomials. The reported logLik/AIC then
  include the prior terms and are penalized quantities, and
  [`anova()`](https://rdrr.io/r/stats/anova.html) comparisons across
  different priors are meaningless.

- quadrature:

  If `TRUE`, marginalize each scalar random effect by adaptive
  Gauss-Kronrod quadrature instead of the Laplace approximation (the
  `glmer(nAGQ = k)` analogue; matches it in tests). Worth it for
  Bernoulli responses with small clusters, where Laplace biases variance
  components. Scalar random-intercept models only, and not with `mi()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html), `REML = TRUE`, or
  `frmtmb_control(profile = TRUE)`. A plain Laplace fit runs first and
  the quadrature tape is built at its optimum: the Gauss-Kronrod
  rescaling is fixed when the tape is built, so the starting point
  decides whether the marginalized objective is finite at all. That fit
  also supplies the conditional modes, which the marginalized objective
  no longer carries, so
  [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) work as usual.

- data2:

  A named list of objects that are not columns of `data`: the adjacency
  matrix of `car()`, the mesh triple of `spde()`, and the matrices of
  `gr(prec = )`, `gr(cov = )` and `equalto()`. This is brms's `data2`
  argument, with one deliberate extension: brms accepts a bare name
  only, while frmtmb also evaluates compound expressions with `data2` in
  front of the data mask, so `gr(g, cov = solve(Q))` finds `Q` there.
  Each structural expression resolves from `data2` first, then `data`,
  then the formula environment; that last step is what a model written
  before `data2` relied on, so old code keeps working. Prefer `data2`:
  its objects are stored on the fit, so
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and a later
  [`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md),
  [`influence()`](https://rdrr.io/r/stats/lm.influence.html) or
  [`update()`](https://rdrr.io/r/stats/update.html) in a fresh session
  do not need the calling environment to still exist.

- dry_run:

  `"spec"` returns the parsed intermediate representation without
  touching `data`; `"frame"` returns the assembled design matrices and
  parameter template without fitting; `"objective"` additionally tapes
  the objective and returns an UNFITTED object carrying it, which is
  what
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  samples when it is given a formula rather than a fit. Methods that
  report a maximum-likelihood quantity refuse on that object.

- verbose:

  Report fit progress; a shortcut for
  `control = frmtmb_control(verbose =)`, whose value wins when both are
  given. See
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)
  for the levels and the output.

## Value

An object of class `frmtmb_fit`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(sleepstudy, package = "lme4")
fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
              data = sleepstudy)
summary(fit)
} # }
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
summary(fit)
#> Family: gaussian 
#> Formula: y ~ x + (1 | g) 
#> Method: ML   nobs: 100 
#> Groups: g, 10 
#> logLik: -145.746  AIC: 299.492  BIC: 309.913 
#> 
#> Random effects:
#>   1 | g 
#>         Name Std.Dev.
#>  (Intercept)  0.55725
#> 
#> Coefficients (mu):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)  1.13568    0.20133  5.6410 1.691e-08
#> x            0.66736    0.11269  5.9222 3.176e-09
#> 
#> Coefficients (sigma):
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept) -0.034743   0.074545 -0.4661   0.6412
fixef(fit)
#> $mu
#> (Intercept)           x 
#>   1.1356757   0.6673571 
#> 
#> $sigma
#> (Intercept) 
#> -0.03474316 
#> 
VarCorr(fit)
#>   1 | g 
#>         Name Std.Dev.
#>  (Intercept)  0.55725

# distributional regression: model sigma too
fit2 <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
anova(fit, fit2)
#> Likelihood-ratio tests
#> 
#>                            Df  logLik    AIC  Chisq Chi Df Pr(>Chisq)   
#> y ~ x + (1 | g)             4 -145.75 299.49                            
#> y ~ x + (1 | g), sigma ~ x  5 -141.62 293.23 8.2614      1    0.00405 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
