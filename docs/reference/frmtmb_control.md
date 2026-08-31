# Control parameters for frmtmb fits

Control parameters for frmtmb fits

## Usage

``` r
frmtmb_control(
  optimizer = "nlminb",
  optCtrl = list(iter.max = 1000, eval.max = 1000),
  restarts = 1,
  grad_tol = 0.001
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

## Value

A list of control settings.
