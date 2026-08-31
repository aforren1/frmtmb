# Control parameters for frmtmb fits

Control parameters for frmtmb fits

## Usage

``` r
frmtmb_control(
  optCtrl = list(iter.max = 1000, eval.max = 1000),
  restarts = 1,
  grad_tol = 0.001
)
```

## Arguments

- optCtrl:

  Control list passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- restarts:

  Number of times to restart the optimizer from the current optimum
  while the gradient remains above `grad_tol`.

- grad_tol:

  Warn (and restart) if the maximum absolute gradient at the optimum
  exceeds this value.

## Value

A list of control settings.
