# Convert draws to a posterior draws object

brms's converters, with brms's arguments and brms's output. The
`as_draws_*()` family takes `variable` (exact names unless `regex`) and
returns a posterior draws object;
[`as_draws()`](https://mc-stan.org/posterior/reference/draws.html) is
brms's
[`as_draws_list()`](https://mc-stan.org/posterior/reference/draws_list.html).
[`as.matrix()`](https://rdrr.io/r/base/matrix.html),
[`as.array()`](https://rdrr.io/r/base/array.html) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) return
brms's unclassed objects, with brms's deprecated `pars` and `subset`
accepted under brms's own warning.

## Usage

``` r
# S3 method for class 'frmtmb_draws'
as_draws(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'frmtmb_draws'
as.data.frame(
  x,
  row.names = NULL,
  optional = TRUE,
  pars = NA,
  variable = NULL,
  draw = NULL,
  subset = NULL,
  ...
)

# S3 method for class 'frmtmb_draws'
as.array(x, pars = NA, variable = NULL, draw = NULL, subset = NULL, ...)

# S3 method for class 'frmtmb_draws'
as_draws_matrix(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'frmtmb_draws'
as_draws_array(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'frmtmb_draws'
as_draws_df(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'frmtmb_draws'
as_draws_list(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'frmtmb_draws'
as_draws_rvars(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

as.mcmc(x, ...)

# S3 method for class 'frmtmb_draws'
as.mcmc(
  x,
  pars = NA,
  fixed = FALSE,
  combine_chains = FALSE,
  inc_warmup = FALSE,
  ...
)

# S3 method for class 'frmtmb_draws'
as.matrix(x, pars = NA, variable = NULL, draw = NULL, subset = NULL, ...)
```

## Arguments

- x:

  A `frmtmb_draws` object.

- variable:

  Variables to keep, by exact name unless `regex`.

- regex:

  If `TRUE`, `variable` is a regular expression.

- inc_warmup:

  Accepted for brms's signature and only `FALSE` is supported: a
  `frmtmb_draws` keeps the post-warmup draws alone.

- ...:

  For [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
  [`as.array()`](https://rdrr.io/r/base/array.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html), the
  `regex`, `fixed` and `inc_warmup` brms passes on; anything else is
  refused by name, rather than silently changing nothing.

- row.names, optional:

  Accepted for the generic and unused, as in brms.

- pars:

  Variables to keep, in brms's spelling: `NA` (the default) for all of
  them, otherwise a character vector matched as a regular expression
  unless `fixed = TRUE`. The argument sits in brms's own second
  position, so `as.mcmc(x, TRUE)` is refused here exactly as brms
  refuses it.

- draw:

  Draws to keep, by index.

- subset:

  brms's deprecated alias of `draw`; it warns.

- fixed:

  If `TRUE`, `pars` is matched by exact name.

- combine_chains:

  If `TRUE`, one `mcmc` object over the pooled draws; otherwise an
  `mcmc.list` with one component per chain, which is what coda's
  diagnostics (`gelman.diag()`) need.

## Value

A
[`posterior::draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html):
one column per sampled variable and one row per draw.

## Examples

``` r
# \donttest{
if (requireNamespace("posterior", quietly = TRUE) &&
    requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(9)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

  # hands the draws to the posterior package, keeping the frmtmb
  # parameter names
  dm <- as_draws(ds)
  posterior::summarise_draws(dm)
  # which is what variables() lists
  head(variables(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] "b_Intercept"      "b_x"              "sigma"            "r_g[1,Intercept]"
#> [5] "r_g[2,Intercept]" "r_g[3,Intercept]"
# }
```
