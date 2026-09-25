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
simulate(
  object,
  nsim = 1,
  seed = NULL,
  re_formula = NULL,
  censored = FALSE,
  newdata = NULL,
  allow_new_levels = FALSE,
  ...
)
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

- re_formula:

  Which group-level terms the draws condition on: `NULL` (default) all
  of them, `NA` (or `~0`, `~1`) none, and a one-sided formula the terms
  it names. A term not kept is redrawn from its estimated distribution
  in each replicate (see Group-level terms).

- censored:

  Apply the fitted `cens()` mechanism to the draws (see Censored
  responses). Ignored without `cens()`.

- newdata:

  Optional data frame to simulate the responses of, instead of the
  fitted rows (see New data).

- allow_new_levels:

  With `newdata`: draw the effect of a grouping level the fit never saw
  from its term's estimated distribution, rather than refuse it. brms's
  spelling.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

A data frame with `nsim` columns and a `"seed"` attribute, with one row
per fitted row, or per row of `newdata`.

## Structured draws

Most families draw each row on its own. Some cannot, and those go
through one implementation that
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
`posterior_predict()` and
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
all reach (see `sim_ctx` in
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)):

- a `mixture(groups = ~g)` draw takes one class per GROUP and then
  simulates each row from its group's component;

- a
  [`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
  draw takes a class per row and then a multivariate normal with that
  class's own covariance;

- a residual correlation term
  ([`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `cosy()`, ...) is
  one multivariate residual draw per group added to the mean predictor,
  so the draws carry the fitted autocorrelation;

- a family from an extension package draws through whatever its
  structure declares, by the same route.

A structured draw covers whole sequences or groups, so
[`trunc()`](https://rdrr.io/r/base/Round.html) rejection cannot resample
single rows within it (every structured model refuses
[`trunc()`](https://rdrr.io/r/base/Round.html) when the frame is
assembled). At `newdata` a structured FAMILY is refused, because its
structure indexes the rows the model was fitted on:
`mixture(groups = )`, a hidden Markov family, a learning family, and
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md),
which `predict(newdata = )` refuses as well. A residual correlation term
is rebuilt on the new rows instead: rows that share a group are drawn
jointly, and a time the fit never saw is refused. A lag is counted in
the FITTED time levels, as the likelihood counts it, so rows at times 2
and 4 are two levels apart even with nothing at time 3 in newdata. brms
counts a lag by a row's position among its group's newdata rows instead,
which makes the correlation of two rows depend on which other rows
newdata holds; frmtmb departs from it on purpose.

## Group-level terms

`re_formula` chooses which group-level terms the draws condition on,
read as
[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
reads it: `NULL` keeps every term, `NA`, `~0` and `~1` keep none, and a
one-sided formula keeps the terms it names (a term the fit does not have
is an error). A term that is kept enters at its estimated effects. A
term that is not kept is REDRAWN from its estimated distribution in
every replicate, which is lme4's unconditional simulation. That is where
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) differ: a prediction
is for an average group, so a dropped term contributes nothing there,
and a simulated response needs a group, so here it gets a new one. When
a formula keeps some columns of a term and drops others, as `~ (1 | g)`
does on a `(1 + x | g)` fit, the dropped columns are drawn given the
kept ones at their estimates.

A population smooth, `gp()` or `hsgp()` curve is not a group-level term
and is never redrawn. A factor-smooth term is, with the other
group-level terms.

## New data

With `newdata` the draws are for its rows. The response column is not
needed. At a grouping level the fit saw, a kept term enters at that
level's estimate, and a redrawn term shares one draw across the rows of
the level, so newdata must carry the grouping column. A level the fit
never saw is an error unless `allow_new_levels = TRUE`, which draws its
effect from the term's estimated distribution, as
[`predict()`](https://rdrr.io/r/stats/predict.html) does. Under
`re_formula = NA` (or `~0`, `~1`) every term is redrawn anyway, so an
unseen level is one more fresh level and needs nothing, except on a fit
with a factor-smooth term: there an unseen level takes the population
curve under `allow_new_levels = TRUE`, as in
[`predict()`](https://rdrr.io/r/stats/predict.html), and is not redrawn.

## Censored responses

On a `cens()` fit the default draws the LATENT, uncensored response: the
model describes the latent distribution, and censoring is a property of
the observation process, not of the response. This matches brms, whose
`posterior_predict()` also ignores `cens()` (and whose
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
single-value representation and is refused too. At `newdata` the same
fitted window applies.

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

# re_formula = NA redraws the group effects, which is the right choice
# for a parametric bootstrap over new groups
sims_m <- simulate(fit, nsim = 5, re_formula = NA, seed = 42)
apply(sims_m, 2, var) > apply(sims, 2, var)
#> sim_1 sim_2 sim_3 sim_4 sim_5 
#>  TRUE FALSE FALSE FALSE  TRUE 

# draws for rows the fit never saw, at a known group and a new one
nd <- data.frame(x = c(-1, 1), g = factor(c("3", "new")))
simulate(fit, nsim = 3, seed = 1, newdata = nd, allow_new_levels = TRUE)
#>   sim_1 sim_2 sim_3
#> 1     1     2     0
#> 2     3     1     1

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
