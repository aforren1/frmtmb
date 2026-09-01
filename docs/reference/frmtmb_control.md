# Control parameters for frmtmb fits

Control parameters for frmtmb fits

## Usage

``` r
frmtmb_control(
  optimizer = "nlminb",
  optCtrl = list(iter.max = 1000, eval.max = 1000),
  restarts = 1,
  grad_tol = 0.001,
  profile = FALSE,
  sparse_x = FALSE,
  autoscale = FALSE,
  check_nlev_1 = c("warning", "ignore", "stop"),
  check_olre = c("warning", "ignore", "stop"),
  verbose = NULL
)
```

## Arguments

- optimizer:

  `"nlminb"` (default), `"optim"` (L-BFGS-B), or a function with
  signature `(par, fn, gr, lower, upper, control)` returning a list with
  elements `par`, `objective`, `convergence` (0 = success), and
  optionally `message` - the hook for optimx, nloptr, and friends
  without frmtmb depending on them. Example:

      nlopt <- function(par, fn, gr, lower, upper, control) {
        r <- nloptr::nloptr(par, fn, gr, lb = lower, ub = upper,
                            opts = list(algorithm = "NLOPT_LD_LBFGS",
                                        xtol_rel = 1e-10, maxeval = 2000))
        list(par = r$solution, objective = r$objective,
             convergence = as.integer(r$status < 0), message = r$message)
      }
      frm(..., control = frmtmb_control(optimizer = nlopt))

- optCtrl:

  Control list passed to the built-in optimizers
  ([`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) /
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html)).

- restarts:

  Number of times to restart the optimizer from the current optimum
  while the gradient remains above `grad_tol`.

- grad_tol:

  Warn (and restart) if the maximum absolute gradient at the optimum
  exceeds this value.

- profile:

  Experimental: move the primary (`beta`) coefficients into the inner
  (Laplace) problem, TMB's `profile` argument - the analog of
  `glmmTMBControl(profile = TRUE)` and `glmer(nAGQ = 0)`. Speeds up
  models with many fixed effects, and like those it is an approximation:
  estimates differ slightly from the exact fit. Coefficient covariance
  comes from the joint precision. Not compatible with `REML = TRUE` or
  `quadrature = TRUE`; profile/uniroot
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  `hypothesis(method = "profile")` need a non-profiled fit.

- sparse_x:

  Build the parametric fixed-effect design matrices as sparse
  [`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
  objects, the analog of `glmmTMB(sparseX =)`. Worth it when a
  many-level fixed factor makes the dense design dominate memory;
  estimates are identical either way.
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) on the
  fit then returns a sparse matrix.

- autoscale:

  Standardize badly scaled continuous predictors internally (the lme4
  \>= 1.1.37 feature): fit a copy of the model with each qualifying
  fixed-effect column centered and scaled (scaled only, in a linear
  predictor without an intercept), map that optimum back to the original
  parameterization exactly, and warm-start the ordinary fit there.
  Reported results are always on the original scale, so every downstream
  method works unchanged; the cost is a second (cheap) optimization.
  Columns qualify when they are parametric and take more than two
  distinct values; intercepts, factor contrasts, smooth bases, and
  mo()/mi() columns are never touched, and the whole step is a silent
  no-op when nothing qualifies. Compatible with `profile = TRUE`. Under
  `priors` or bounds, the first stage applies them to the scaled
  coefficients; the second stage is the fit that is reported.

- check_nlev_1:

  What to do about a scalar random-effect term whose grouping factor has
  a single level: `"warning"` (default), `"ignore"`, or `"stop"`,
  following lme4's `lmerControl()` check vocabulary. Such a term has no
  variance to estimate - the single level is absorbed by the intercept -
  and its standard deviation collapses to zero. Structured blocks over
  several terms per level (`ar1()`, `us()`, the spatial covariance
  structures) are never flagged: one grouping level there is a single
  realization of a field, which is the normal way to write them.

- check_olre:

  What to do about an observation-level random effect - one level per
  row - on a gaussian, student or lognormal response: `"warning"`
  (default), `"ignore"`, or `"stop"`. Its variance is confounded with
  the residual standard deviation, so only their sum is identified and
  the split between them is arbitrary. The check is skipped when `sigma`
  is not free to absorb it - a `se()` response or a constant `sigma` -
  which is the random-effects meta-analysis, and for discrete families,
  where an observation-level term is the usual overdispersion model.

- verbose:

  Report fit progress through
  [`message()`](https://rdrr.io/r/base/message.html), one terse line per
  stage with its elapsed seconds, so a slow fit shows where the time
  went. `FALSE` (default) is silent and costs nothing. `TRUE` (or `1`)
  reports parsing, frame assembly, the autoscale pre-fit, tape
  construction, each optimizer run and restart, `sdreport()` when
  `se = TRUE`, and a closing line with the objective, the maximum
  absolute gradient, and the number of convergence warnings. The fit
  itself opens with a line naming the family, the mode (ML or REML, plus
  profile, quadrature, autoscale, priors), and the optimizer. `2` adds
  the optimizer's own trace, by setting `optCtrl$trace` unless you set
  it yourself. That trace is printed by
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) /
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) to standard
  output, not through
  [`message()`](https://rdrr.io/r/base/message.html), and a custom
  optimizer function receives `optCtrl` unchanged, so `verbose` does not
  reach it.

  Use [`suppressMessages()`](https://rdrr.io/r/base/message.html) to
  silence a verbose fit. The refit loops in
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md),
  [`influence()`](https://rdrr.io/r/stats/lm.influence.html), and
  [`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)
  force `verbose` off, so a verbose fit does not make them print
  hundreds of lines;
  [`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md) and
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  report normally.

## Value

A list of control settings.

## Examples

``` r
set.seed(1)
n <- 200
dd <- data.frame(x = rnorm(n), g = factor(rep(1:10, 20)))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# another optimizer, with its own control list
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           control = frmtmb_control(optimizer = "optim",
                                    optCtrl = list(maxit = 500)))
fit$opt$convergence
#> [1] 0

# a tighter gradient criterion, with restarts from the current
# optimum until it is met
frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
    control = frmtmb_control(grad_tol = 1e-4, restarts = 3))
#> frmtmb fit: y ~ x + (1 | g) 
#> Family: gaussian   Method: ML 
#> logLik: -296.63  AIC: 601.26  nobs: 200 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      1.4416      0.5196 
#>  sigma:
#> (Intercept) 
#>    -0.01237 
#> 
#> Random effects:
#>   1 | g 
#>         Name Std.Dev.
#>  (Intercept)  0.99749

# badly scaled predictors: fit an internally standardized copy first,
# then warm-start the reported fit from it
dd$xbig <- dd$x * 1e5
frm(bf(y ~ xbig + (1 | g)) + gaussian(), data = dd,
    control = frmtmb_control(autoscale = TRUE))
#> frmtmb fit: y ~ xbig + (1 | g) 
#> Family: gaussian   Method: ML 
#> logLik: -296.63  AIC: 601.26  nobs: 200 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)        xbig 
#>   1.442e+00   5.196e-06 
#>  sigma:
#> (Intercept) 
#>    -0.01237 
#> 
#> Random effects:
#>   1 | g 
#>         Name Std.Dev.
#>  (Intercept)  0.99749

# the object is a plain list, so it can be built once and reused
ctrl <- frmtmb_control(check_nlev_1 = "ignore")
ctrl$optimizer
#> [1] "nlminb"
```
