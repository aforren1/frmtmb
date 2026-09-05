# Numerical controls for [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

The generalized drift-diffusion likelihood has no closed form, so every
evaluation solves a partial differential equation on a grid. These are
the grid and what is done with the answer.

## Usage

``` r
gddm_control(
  dt = 0.01,
  ny = 201L,
  t_max = NULL,
  renormalize = TRUE,
  max_ndt = NULL,
  tridiagonal = c("recorded", "atomic")
)
```

## Arguments

- dt:

  Time step, in the units of the response. The default matches PyDDM's.

- ny:

  Number of interior spatial nodes between the boundaries. Odd is the
  natural choice, putting an unbiased start on a node.

- t_max:

  End of the modeled window, in the units of the response. `NULL`, the
  default, takes the smallest multiple of `dt` strictly above the
  largest response time. Set it explicitly to the experiment's response
  deadline when there is one, and whenever you will
  [`predict()`](https://rdrr.io/r/stats/predict.html) on new data whose
  largest response time differs.

- renormalize:

  Divide each condition's pair of densities by their own total mass.
  Leave this on. The discretized solve loses mass in a way that depends
  on the parameters, so a likelihood that does not renormalize rewards
  parameter values that absorb faster and biases the leak and the
  boundary height. Turning it off is provided so that the size of that
  bias can be measured, not because it is ever the better model.

- max_ndt:

  Upper bound for the non-decision time, in the units of the response.
  `NULL` takes the smallest response time. See
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md),
  whose `ndt` link this shares.

- tridiagonal:

  How the tridiagonal solve inside each step reaches the tape.
  `"recorded"`, the default, lets 'RTMB' record the sweep: the tape is
  slow to build and fast to run. `"atomic"` collapses the sweep into one
  node with a hand-written adjoint: the tape is quick to build and
  slower to run, because each evaluation calls back into R. Measured at
  the shipped grid over six conditions, the atomic builds about twelve
  times faster and evaluates about twelve times slower, so it pays only
  when the tape build dominates: a very fine grid, many conditions, or a
  fit that stops after a few tens of iterations. Both compute the same
  derivative.

## Value

A `gddm_control`.

## Cost

One evaluation costs one solve per **condition**, not per trial, so the
number of distinct parameter settings in the design is what this scales
with. One solve is `t_max / dt` time steps of work proportional to `ny`.
Halving `dt` doubles the cost and, because the scheme is second order in
time, divides the time-discretization error by about four. Doubling `ny`
doubles the cost of a step.

The tape is built once per fit and then evaluated once per iteration. At
the default `tridiagonal = "recorded"` the tape grows with `t_max / dt`,
with `ny` and with the number of conditions, which makes the build the
larger cost for a short fit and the evaluation the larger cost for a
long one. `tridiagonal = "atomic"` moves the balance the other way.

## Accuracy

With a constant drift and boundaries that do not move the model is the
Wiener model, so the solver can be checked against the closed form. At
the default grid, across a range of drifts, separations and starting
points and over decision times from 0.2 s on, the two agree to better
than 0.01 in the log density; a coarser grid is worse and a finer one
better, and the package's tests pin all three.

Short decision times are the limit. An implicit scheme spreads a little
probability everywhere immediately, where the true first-passage density
is exponentially small, so within a few steps of zero the density is far
larger than the truth. It is small in absolute terms and a lapse
component floors it, but a fit should not be asked to read it.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

## Examples

``` r
gddm_control(dt = 0.005, ny = 301)
#> <gddm_control: dt = 0.005, ny = 301, t_max = from data, renormalize = TRUE>
```
