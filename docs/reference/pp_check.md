# Predictive check against simulated responses

The frequentist analog of brms's `pp_check()`: responses are simulated
from the fitted model (marginally over the random effects) and handed to
the corresponding bayesplot `ppc_*` function (bayesplot must be
installed, but not necessarily attached).

## Usage

``` r
pp_check(object, ...)

# S3 method for class 'frmtmb_fit'
pp_check(
  object,
  type = "dens_overlay",
  ndraws = 10,
  prefix = c("ppc", "ppd"),
  group = NULL,
  x = NULL,
  newdata = NULL,
  resp = NULL,
  ...,
  re_formula = NA
)
```

## Arguments

- object:

  A `frmtmb_fit` for a univariate model.

- ...:

  Passed to the `ppc_*` function.

- type:

  The bayesplot check, i.e. the part after `ppc_` (`"dens_overlay"`,
  `"hist"`, `"stat"`, `"stat_grouped"`, `"scatter_avg"`, ...). With
  `prefix = "ppc"` a name that
  [`bayesplot::available_ppc()`](https://mc-stan.org/bayesplot/reference/available_ppc.html)
  does not list is refused, as brms refuses it.

- ndraws:

  Number of simulated response vectors.

- prefix:

  `"ppc"` (the default) plots the observed response against the
  simulated ones; `"ppd"` plots the simulated ones alone, through
  bayesplot's `ppd_*` function of the same name.

- group:

  The name of a model variable to stratify by, for the `*_grouped`
  types, which need it. The name is looked up in the model frame, as
  brms looks it up in the model's data, so a column the model does not
  use is refused.

- x:

  The name of a model variable for the types that take an `x`
  (`"intervals"`, `"ribbon"`, `"error_scatter_avg_vs_x"`, ...), looked
  up as `group` is.

- newdata, resp:

  brms's arguments. A fit simulates only its own rows, so a `newdata` is
  refused. `resp` is accepted and ignored, which is what brms does with
  it on a model that has one response; a multivariate fit is refused
  before `resp` could select one.

- re_formula:

  The random-effect switch, in brms's spelling (`pp_check()` is a brms
  function). On a fit it is passed to
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and defaults to
  `NA`, which simulates new random effects; on draws it is passed to
  `posterior_predict()` and defaults to `NULL`, because a draw already
  carries its own. lme4's `re.form` is refused. brms honors it on
  `pp_check()` and warns that it ignored it, which is a leak through its
  dots rather than a decision to copy.

## Value

A ggplot object, as returned by the bayesplot `ppc_*` function that
`type` selects.

## Types a fit cannot draw

The `loo_*` types weight posterior draws by Pareto-smoothed importance
sampling, and a maximum-likelihood fit has no posterior draws, so they
are refused on a fit. `pp_check()` on the draws of
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
is the route to them. `"error_binned"` is refused on an ordinal,
categorical or multinomial fit, as brms refuses it.

## Examples

``` r
if (requireNamespace("bayesplot", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  # the observed density against draws from the fit
  pp_check(fit, ndraws = 20)

  # any bayesplot ppc_* check, named by its suffix. A statistic the
  # model was not fitted to is the informative one: here, the share
  # of zeros, which is how zero inflation shows up.
  pp_check(fit, type = "stat", stat = function(y) mean(y == 0),
           ndraws = 50)
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
