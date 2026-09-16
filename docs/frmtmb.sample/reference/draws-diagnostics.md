# Sampler diagnostics and MCMC plots

`rhat()` and `neff_ratio()` are the convergence diagnostics brms
reports, computed by the posterior package on these draws: `rhat()` is
the rank-normalized split-R-hat and `neff_ratio()` is
`min(ess_bulk, ess_tail) / ndraws`. `nuts_params()` and
`log_posterior()` delegate to bayesplot's `stanfit` methods on the
`stanfit` inside the draws object, which is what brms does too, so every
`bayesplot::mcmc_nuts_*()` display works. `mcmc_plot()` is brms's
spelling for "call a bayesplot `mcmc_*` function on these draws";
[`pairs()`](https://rdrr.io/r/graphics/pairs.html) is
[`bayesplot::mcmc_pairs()`](https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html).

## Usage

``` r
mcmc_plot(object, ...)

# S3 method for class 'frmtmb_draws'
mcmc_plot(
  object,
  pars = NA,
  type = "intervals",
  variable = NULL,
  regex = FALSE,
  fixed = FALSE,
  ...
)

# S3 method for class 'frmtmb_draws'
pairs(x, variable = NULL, ...)

nuts_params(object, ...)

# S3 method for class 'frmtmb_draws'
nuts_params(object, ...)

log_posterior(object, ...)

# S3 method for class 'frmtmb_draws'
log_posterior(object, ...)

rhat(x, ...)

# S3 method for class 'frmtmb_draws'
rhat(x, pars = NULL, regex = FALSE, ...)

neff_ratio(object, ...)

# S3 method for class 'frmtmb_draws'
neff_ratio(object, pars = NULL, regex = FALSE, ...)
```

## Arguments

- object, x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- ...:

  For `mcmc_plot()` and
  [`pairs()`](https://rdrr.io/r/graphics/pairs.html), passed to the
  bayesplot function; for `nuts_params()` and `log_posterior()`, passed
  to bayesplot's own `stanfit` method, so `nuts_params(x, "stepsize__")`
  reaches its `pars`. **`rhat()` and `neff_ratio()` read nothing from
  it**: they take `pars` and `regex` and no more, so an argument they do
  not have, such as `variable = "x"`, is accepted and IGNORED rather
  than refused, and the whole set of variables comes back. brms errors
  on that call. A refusal is coming from `frm_check_dots()` (plan item
  2.5e); until it lands, this is the accurate statement of what happens.

- pars:

  Which variables to report, in brms's spelling. The rule differs by
  method; see *Two different `pars` rules, both brms's*.

- type:

  The bayesplot function to call, without the `mcmc_` prefix (default
  `"intervals"`).

- variable:

  For `mcmc_plot()` and
  [`pairs()`](https://rdrr.io/r/graphics/pairs.html), the variables to
  use, by name; it defaults to everything except the group-level modes
  and `lp__`. `rhat()` and `neff_ratio()` do not take it, because brms's
  do not: their selector is `pars`. Naming it on either of those two is
  silently ignored today; see `...`.

- regex:

  For `rhat()` and `neff_ratio()`, `TRUE` makes `pars` a regular
  expression; for `mcmc_plot()`, it makes `variable` one.

- fixed:

  For `mcmc_plot()`, `TRUE` matches `pars` by exact name rather than as
  a regular expression.

## Value

A ggplot object, or the diagnostic data frame / vector bayesplot
returns.

## Details

All of these report the frmtmb draws-side parameter names (no
parentheses), not Stan's `par[1]`, except `nuts_params()`, whose rows
are the sampler's own quantities and not model parameters.

## Which R-hat this is

`rhat()` follows brms, whose `rhat.brmsfit()` is
`posterior::summarise_draws(rhat = posterior::rhat)`. That is the
rank-normalized split-R-hat of Vehtari et al. (2021), the maximum of the
bulk and tail quantities, and it is NOT the classic split-R-hat that
`rstan::summary()` reports. The two disagree by about the size of the
excess over 1 that either of them reports, so `rhat(ds)` and
`ds$stanfit` do not agree and are not meant to.

`summary(ds)` agrees with `rhat(ds)`, because it reports the same three
posterior quantities `brms:::summary.brmsfit()` reports, under the same
column names: `Rhat`, `Bulk_ESS` and `Tail_ESS`. The sampler's own
classic split-R-hat and `n_eff` are in
`rstan::summary(ds$stanfit)$summary` for anyone who wants them.

## Two different `pars` rules, both brms's

brms spells two different selectors `pars`, and this page carries both
because it documents methods on either side of the line.

`mcmc_plot()` takes brms's deprecated alias of `variable`: `NA` is every
variable, a string is a regular expression unless `fixed = TRUE`, and
anything that is neither `NA` nor character is refused.
[`as.mcmc()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-as_draws.md)
and
[`posterior_interval()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/sample-posterior_summary.md)
take the same one.

`rhat()` and `neff_ratio()` do not. brms's `rhat.brmsfit()` passes
`variable = pars` straight to
[`as_draws_array()`](https://mc-stan.org/posterior/reference/draws_array.html),
so `NULL` is every variable, a string is an EXACT variable name,
`regex = TRUE` makes it a regular expression, and a name that is not
there is an error. These two follow that rule, which is why their
default is `NULL` and not `NA`.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    requireNamespace("bayesplot", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  mcmc_plot(ds)
  mcmc_plot(ds, type = "trace", variable = "x")
  head(rhat(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>       Intercept               x sigma_Intercept            b[1]            b[2] 
#>        1.014249        1.010188        1.003701        1.012647        1.006744 
#>            b[3] 
#>        1.004775 
# }
```
