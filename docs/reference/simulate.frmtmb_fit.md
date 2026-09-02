# Simulate responses from a frmtmb fit

A [`trunc()`](https://rdrr.io/r/base/Round.html)ed response simulates by
rejection within its bounds, so every draw lies in `[lb, ub]` and
posterior-predictive checks
([`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md),
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md))
see the same support the likelihood was normalized on.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
simulate(object, nsim = 1, seed = NULL, re.form = NULL, censored = FALSE, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- nsim:

  Number of simulated response vectors.

- seed:

  Optional RNG seed. Follows the
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) contract:
  the global RNG state is restored afterwards, and the seed used is
  attached as the `"seed"` attribute.

- re.form:

  `NULL` (default) conditions on the estimated random effects; `NA`
  redraws them from their estimated distribution (marginal simulation).

- censored:

  Apply the fitted `cens()` mechanism to the draws (see Censored
  responses). Ignored without `cens()`.

- ...:

  Unused.

## Value

A data frame with `nsim` columns and a `"seed"` attribute.

## Structured draws

Most families draw each row on its own. Some cannot, and those go
through one implementation that
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
and
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
all reach (see `sim_ctx` in
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)):

- an [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) draw
  walks the hidden Markov chain forward per sequence and then emits from
  each row's state;

- a `mixture(groups = ~g)` draw takes one class per GROUP and then
  simulates each row from its group's component;

- a
  [`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
  draw takes a class per row and then a multivariate normal with that
  class's own covariance;

- a residual correlation term
  ([`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `cosy()`, ...) is
  one multivariate residual draw per group added to the mean predictor,
  so the draws carry the fitted autocorrelation.

A structured draw covers whole sequences or groups, so
[`trunc()`](https://rdrr.io/r/base/Round.html) rejection cannot resample
single rows within it (every structured model refuses
[`trunc()`](https://rdrr.io/r/base/Round.html) when the frame is
assembled) and `posterior_predict(newdata =)` is refused: the structure
indexes the rows the model was fitted on.

## Censored responses

On a `cens()` fit the default draws the LATENT, uncensored response: the
model describes the latent distribution, and censoring is a property of
the observation process, not of the response. This matches brms, whose
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
also ignores `cens()` (and whose
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
therefore drops the censored rows). The draws are then not comparable
with the observed values on censored rows, which is why
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
and `residuals(type = "osa")` refuse or skip them.

`censored = TRUE` applies the fitted censoring mechanism to each draw
instead, so the draws are directly comparable with the observed data:
every draw is recorded at the edge of the observation window it falls
outside, capped above by the right-censoring point and below by the
left-censoring point. Those points are the response values of the
censored rows, and they must be the same for every censored row on a
side (type-I censoring): with row-varying censoring times an uncensored
row's censoring point is unknown, so the mechanism cannot be applied to
its draws and the call is refused. Interval censoring has no
single-value representation and is refused too.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

# one column per draw; the seed used is attached
sims <- simulate(fit, nsim = 5, seed = 42)
str(sims)
#> 'data.frame':    100 obs. of  5 variables:
#>  $ sim_1: int  2 3 0 4 1 2 2 1 2 5 ...
#>  $ sim_2: int  1 0 0 2 3 6 2 3 2 0 ...
#>  $ sim_3: int  1 1 2 2 0 2 4 2 1 2 ...
#>  $ sim_4: int  0 1 0 2 2 5 1 1 1 3 ...
#>  $ sim_5: int  0 1 1 2 2 1 4 1 0 3 ...
#>  - attr(*, "seed")= num 42
#>   ..- attr(*, "kind")=List of 3
#>   .. ..$ : chr "Mersenne-Twister"
#>   .. ..$ : chr "Inversion"
#>   .. ..$ : chr "Rejection"
attr(sims, "seed")
#> [1] 42
#> attr(,"kind")
#> attr(,"kind")[[1]]
#> [1] "Mersenne-Twister"
#> 
#> attr(,"kind")[[2]]
#> [1] "Inversion"
#> 
#> attr(,"kind")[[3]]
#> [1] "Rejection"
#> 

# re.form = NA redraws the group effects, which is the right choice
# for a parametric bootstrap over new groups
sims_m <- simulate(fit, nsim = 5, re.form = NA, seed = 42)
apply(sims_m, 2, var) > apply(sims, 2, var)
#> sim_1 sim_2 sim_3 sim_4 sim_5 
#>  TRUE FALSE FALSE FALSE  TRUE 

# a posterior-predictive check by hand: does the fit reproduce the
# share of zeros in the data?
mean(dd$y == 0)
#> [1] 0.29
colMeans(simulate(fit, nsim = 20, seed = 1) == 0)
#>  sim_1  sim_2  sim_3  sim_4  sim_5  sim_6  sim_7  sim_8  sim_9 sim_10 sim_11 
#>   0.30   0.24   0.35   0.30   0.25   0.25   0.23   0.24   0.33   0.32   0.27 
#> sim_12 sim_13 sim_14 sim_15 sim_16 sim_17 sim_18 sim_19 sim_20 
#>   0.27   0.30   0.37   0.27   0.29   0.24   0.29   0.28   0.34 
```
