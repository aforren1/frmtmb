# Check the Laplace/Wald approximation against NUTS

Samples the fitted objective (see
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md))
and compares the ML estimates and sdreport standard errors against
posterior means and SDs. Close agreement supports the Laplace
approximation and Wald intervals; a posterior SD much larger than the
Wald SE, or a shifted mean, flags parameters where they are unreliable
(typically variance components with few groups).

## Usage

``` r
check_laplace(fit, chains = 2, iter = 1000, ...)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- chains, iter:

  Passed to
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Passed to
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

## Value

A data frame (one row per outer parameter) with columns `ml`,
`post_mean`, `wald_se`, `post_sd`, `z_shift` ((post_mean - ml)/post_sd)
and `sd_ratio` (post_sd/wald_se).

## Details

**It samples the density the fit maximized, and no other.** It hands
Stan `fit$obj` as it stands. For an ordinary fit that objective is the
LIKELIHOOD, so the run is flat; for a MAP fit (`frm(prior = )`) the
penalty is taped INTO that objective, so the run carries it, which is
right: the mode and the Wald standard errors this function compares
against come from the same penalized Hessian. "Flat" here means no prior
ADDED, not the fit's own penalty stripped - for the bare likelihood of a
MAP fit, write `frm_sample(fit, prior = "flat")` instead, and read it as
a different question.

What this function never adds is
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)'s
DEFAULT priors. Those defaults apply on both of
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)'s
routes, and this function opts out of them explicitly
(`frm_sample(.diagnostic = TRUE)`) rather than inheriting whatever the
default is: a prior nothing in the fit ever saw would change the very
thing being measured.

That is also why it samples CENTERED on an ordinary fit.
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)'s
non-centered parameterization (see Reparameterization there) is offered
only for blocks whose variance parameters carry a prior, and with the
defaults off none do: the flat prior that makes the comparison
meaningful is exactly the one that leaves a flat tail at `sd = 0` for a
non-centered chain to walk into. So the reparameterization default costs
this function nothing and changes nothing about it. Give the variance
parameters a prior through `prior =` and the run non-centers; but then
it is measuring the Laplace approximation of a different posterior,
which is usually not the question.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
# a binary GLMM with small clusters: the regime where the Laplace
# approximation and Wald intervals are least reliable
set.seed(4)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:30, 4)))
dd$y <- rbinom(120, 1,
               plogis(0.3 + 0.5 * dd$x + rnorm(30, 0, 1)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)

cl <- check_laplace(fit, chains = 1, iter = 500, refresh = 0)
cl
# |z_shift| well above 0 or sd_ratio far from 1 marks the parameters
# whose Wald interval to replace with a profile or bootstrap one
cl[abs(cl$z_shift) > 0.3 | cl$sd_ratio > 1.3, ]
}
#> frm_sample(): sampling stays centered: no random-effect block of this model has a non-centered form:
#>   1 | g [us]: its variance parameter has a flat prior here, and a non-centered chain walks the flat tail that opens at sd = 0. Give it a prior, set_prior(class = "sd"), which frm_sample() supplies by default unless prior = "flat" turned it off
#> Warning: The largest R-hat is 1.08, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> check_laplace(): the chain mixed too poorly to judge the approximation (bulk ESS under 100 for Intercept, theta_1). Rerun with more iterations before reading z_shift or sd_ratio
#> [1] parameter ml        post_mean wald_se   post_sd   z_shift   sd_ratio 
#> [8] ess_bulk 
#> <0 rows> (or 0-length row.names)
# }
```
