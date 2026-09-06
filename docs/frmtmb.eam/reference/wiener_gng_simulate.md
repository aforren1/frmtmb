# Simulate from a go/no-go diffusion model

Draws from the generative process: run the two-boundary diffusion, and
record a response only if it reached the upper boundary before the
deadline. Nothing here evaluates the go density or the no-go
probability, so the simulator is an independent statement of the same
model and can be used to check the likelihood rather than merely to
agree with it.

## Usage

``` r
wiener_gng_simulate(
  n,
  mu = 1,
  bs = 1.4,
  ndt = 0.25,
  bias = 0.5,
  deadline = 1.5,
  sv = 0,
  sz = 0,
  st = 0
)
```

## Arguments

- n:

  Number of trials.

- mu:

  Drift rate.

- bs:

  Boundary separation.

- ndt:

  Non-decision time.

- bias:

  Relative start point in `(0, 1)`.

- deadline:

  The response deadline. One number, or one per trial.

- sv, sz, st:

  Across-trial variability in the drift rate, the relative start point
  and the non-decision time, as
  [`ddm_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ddm_simulate.md)
  takes them. Zero, the default, is the plain model.

## Value

A data frame with `responded` (1 for a go trial), `rt` (the response
time on a go trial, and the deadline on a no-go trial, which is the
placeholder
[`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
documents) and `deadline`.

## Examples

``` r
set.seed(1)
dat <- wiener_gng_simulate(1000, mu = 1.0, bs = 1.4, ndt = 0.25,
                           deadline = 1.5)
mean(dat$responded)
#> [1] 0.759
summary(dat$rt[dat$responded == 1])
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2950  0.4431  0.5709  0.6473  0.7767  1.4981 
```
