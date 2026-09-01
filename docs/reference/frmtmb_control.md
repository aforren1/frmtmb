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
  autoscale = FALSE
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

## Value

A list of control settings.
