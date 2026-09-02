# Model diagnostics

Diagnostics for a frmtmb fit answer three separate questions: did the
optimizer find a real optimum, does the fitted family describe the data,
and are the reported uncertainties trustworthy. Each has its own tool.

``` r

library(frmtmb)
set.seed(1)
n <- 400
dd <- data.frame(x = rnorm(n), g = factor(rep(1:20, n / 20)))
dd$y <- rpois(n, exp(0.5 + 0.4 * dd$x + rnorm(20, 0, 0.5)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
```

## Did the optimizer converge?

[`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
reports the convergence code, the largest gradient component at the
optimum, whether the Hessian is positive definite, and covariance
parameters near a boundary (a variance component collapsing to zero is
the common case):

``` r

diagnose(fit)
#> Optimizer convergence code: 0 (relative convergence (4)) 
#> Max |gradient|: 6.943e-05 at (Intercept) 
#> Hessian positive definite: TRUE 
#> No convergence problems detected
```

## Convergence problems

The warning “Large maximum absolute gradient at the optimum” means the
optimizer stopped at a point where the objective still has slope above
`frmtmb_control(grad_tol = 1e-3)`. Work through these in order; each
step tells you something even when it does not fix the warning.

1.  **Identify the parameter.** `diagnose(fit)` prints the largest
    gradient component and its name. A coefficient of a continuous
    predictor points to scaling; a `theta_*` component points to a
    variance-parameter problem, and
    [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
    flags extreme values (a log-SD far below zero is a variance
    collapsing to its boundary, which is a model issue, not an optimizer
    issue).
2.  **Fix predictor scaling.** For a flagged coefficient,
    `frmtmb_control(autoscale = TRUE)` standardizes badly scaled columns
    internally and judges convergence in the rescaled units; rescaling
    the variable in the data works equally well.
3.  **Restart harder.** `frmtmb_control(restarts = 3)` re-launches the
    optimizer from the current optimum until the gradient passes;
    `optCtrl = list(iter.max = 5000, eval.max = 5000)` lifts the
    iteration caps for slow problems.
4.  **Compare optimizers.** `frm_allfit(fit)` refits with nlminb, optim,
    bobyqa, and nloptr and tabulates log-likelihoods. All optimizers
    agreeing on the log-likelihood to several decimals while one warns
    usually means a flat ridge near the optimum, not a wrong answer;
    disagreement means the warning is real.
5.  **Look at the profile.** `plot(profile(fit, "theta_1"))` shows
    whether the flagged parameter sits in a flat valley (benign) or on a
    slope (not converged).
6.  **Boundary fits.** When a random-effect SD is collapsing to zero,
    the likelihood surface is one-sided and the gradient check can never
    fully settle. Either simplify the random-effect structure, or keep
    the term and regularize with
    `priors = set_prior("exponential(1)", class = "sd")` (a MAP fit; the
    same move as `blme`’s covariance priors).
7.  **Judge marginal cases.** A gradient just above the default `1e-3`
    (say 0.001-0.01) with stable estimates across restarts and
    optimizers is typically flatness at machine precision for that
    problem’s scale, not non-convergence. The threshold is a heuristic;
    `frmtmb_control(grad_tol = )` tunes it, and steps 4-5 are the
    evidence for deciding.

## Does the family fit the data?

Three residual types, in increasing order of statistical care:

- `residuals(fit)` is response-scale (`y - E[y]`), in the response’s
  units. Good for effect-size intuition, useless as a distributional
  check for counts.
- `residuals(fit, type = "pearson")` standardizes by the family
  variance.
- `residuals(fit, type = "osa")` gives one-step-ahead quantile residuals
  via \[TMB::oneStepPredict()\]: each observation’s CDF position given
  the previous ones, mapped through the normal quantile function. They
  are standard normal under a correctly specified model regardless of
  family, and remain valid under correlated observations, where pearson
  residuals mislead.

``` r

r <- residuals(fit, type = "osa")
qqnorm(r); abline(0, 1)
```

![Normal quantile-quantile plot of the one-step-ahead residuals. The
points follow the identity line over most of the range, with small
departures in the tails.](diagnostics_files/figure-html/osa-qq-1.png)

Simulation-based residuals through DHARMa cover the same ground with a
richer test suite (dispersion, zero-inflation, outliers) and better
plots:

``` r

dh <- dharma_residuals(fit, nsim = 200, seed = 1)
plot(dh)
```

![The standard DHARMa panel pair. The left panel is a uniform
quantile-quantile plot of the scaled residuals with its three
goodness-of-fit tests. The right panel plots the scaled residuals
against the rank-transformed model
predictions.](diagnostics_files/figure-html/dharma-plot-1.png)

``` r

DHARMa::testDispersion(dh, plot = FALSE)
#> 
#>  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
#>  simulated
#> 
#> data:  simulationOutput
#> dispersion = 1.0486, p-value = 0.59
#> alternative hypothesis: two.sided
```

A misspecified model shows up immediately; here the same data fit
without the group effect is overdispersed:

``` r

fit_bad <- frm(bf(y ~ x) + poisson(), data = dd)
DHARMa::testDispersion(dharma_residuals(fit_bad, nsim = 200, seed = 1),
                       plot = FALSE)
#> 
#>  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
#>  simulated
#> 
#> data:  simulationOutput
#> dispersion = 1.661, p-value < 2.2e-16
#> alternative hypothesis: two.sided
```

Two quicker looks. `plot(fit)` draws Pearson residuals against fitted
values and their normal QQ plot. And `pp_check(fit)` overlays the
observed response with responses simulated from the fitted model,
through bayesplot’s `ppc_*` functions (brms users know this display):

``` r

bayesplot::pp_check(fit, ndraws = 20)
```

![Posterior predictive check. A thick dark density curve is the observed
count response. Twenty thin light curves are responses simulated from
the fitted model. The simulated curves surround the observed
one.](diagnostics_files/figure-html/pp-check-1.png)

## Are the standard errors trustworthy?

Wald standard errors and the Laplace approximation share failure modes:
variance components with few groups, binary data with tiny clusters.
[`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md)
samples the fitted objective with NUTS (initialized at the ML mode) and
compares posterior means and SDs against the ML estimates and Wald SEs;
large `z_shift` or `sd_ratio` far from 1 flags parameters where Wald
intervals should not be trusted - use `confint(method = "profile")` or
the posterior itself for those.

``` r

check_laplace(fit, chains = 2, iter = 1000)
```

The drawn posterior is also a diagnostic object in its own right;
bayesplot works on the draws:

``` r

ds <- frm_sample(fit, chains = 4,
                 priors = set_prior("normal(0, 5)", class = "b") +
                   set_prior("exponential(1)", class = "sd"))
bayesplot::mcmc_intervals(as.matrix(ds),
                          pars = c("(Intercept)", "x"))
```

Without priors,
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
samples with flat improper priors on the internal parameters; fine as a
diagnostic, fragile as inference. Supply priors (brms-style
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md))
for anything more.

## Variance components on their natural scale

[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
reports random-effect SDs and correlations with delta-method intervals
on interpretable scales:

``` r

confint_varcorr(fit)
#>   block        term type  estimate       lwr       upr
#> 1 1 | g (Intercept)   sd 0.4783034 0.3383263 0.6761937
```
