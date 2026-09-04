# Expected-value and predictive draws from sampled parameters

`posterior_epred()` evaluates the response-scale expectation per draw;
`posterior_predict()` additionally simulates responses from the family,
giving the posterior predictive distribution. Both condition on each
draw's own random effects (`re_formula = NA` drops them; `re.form` is
the accepted alias, see *Argument spellings*).

## Usage

``` r
posterior_epred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_epred(
  object,
  newdata = NULL,
  resp = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  ndraws = NULL,
  ...
)

posterior_linpred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  resp = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  dpar = NULL,
  ndraws = NULL,
  ...
)

posterior_predict(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_predict(
  object,
  newdata = NULL,
  resp = NULL,
  re_formula = arg_unset(),
  re.form = arg_unset(),
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- ...:

  Unused.

- newdata, resp:

  As in
  [`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html).

- re_formula:

  The random-effect switch, in brms's spelling: `NULL` (the default)
  conditions on each draw's own random effects, `NA` or `~0` gives the
  population-level quantity. Its meaning is
  [`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html)'s
  `re.form`; see *Argument spellings*.

- re.form:

  lme4's spelling of `re_formula`, accepted as an alias. Pass one or the
  other, not both.

- ndraws:

  Number of draws to use (default: all).

- transform:

  For `posterior_linpred()`: if `TRUE`, apply the inverse link (the
  value of the `mu` dpar on its natural scale, brms's convention; unlike
  `posterior_epred()` this is not the response mean for zero-inflated
  and similar families).

- dpar:

  For `posterior_linpred()`: which distributional parameter's linear
  predictor to evaluate.

## Value

A draws-by-observations matrix; for a categorical outcome
`posterior_epred()` returns a draws-by-observations-by-categories array
(see the section below).

## Categorical outcomes

An ordinal family predicts a DISTRIBUTION per observation, not one
number: each draw's `predict(type = "response")` is an `n x K` matrix of
category probabilities. Those stack into a 3-D
`draws x observations x categories` array. `dimnames` are
`list(NULL, <observation names or NULL>, <category levels>)`, so
`ep[, , "high"]` is the draws-by-observations matrix for one category
and `ep[k, , ]` is draw `k`'s own `n x K` prediction, the matrix
`predict(type = "response")` returns. Every `ep[k, i, ]` sums to 1 for
an ordinal family.

This is brms's convention:
[`?brms::posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.html)
documents "an S x N x C array" for categorical and ordinal models and an
S x N matrix otherwise, and frmtmb follows brms spelling for brms-origin
functions. Any family whose per-draw response-scale prediction is a
matrix takes the array shape; every family that predicts one number per
observation keeps the plain `draws x observations` matrix.

`posterior_predict()` is unaffected for an ordinal or categorical family
(it draws one category per observation), and so is
`posterior_linpred()`, which is a statement about one distributional
parameter and stays an `n`-column matrix of the latent predictor. What
does take the array shape in `posterior_predict()` is a matrix-valued
RESPONSE:
[`frmtmb::multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
counts,
[`frmtmb::mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html)
draws and
[`frmtmb.latent::lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.html)
item codes give one row per observation, so the draws stack into
`draws x observations x columns`.

## Structured draws

`posterior_predict()` uses the same simulator
[`simulate()`](https://rdrr.io/r/stats/simulate.html) does, including
the structured families
([`frmtmb.latent::hmm()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/hmm.html),
`mixture(groups = )`,
[`frmtmb::mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html))
and residual correlation terms; see the Structured draws section of
[`frmtmb::simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.html).
Those draws index the rows the model was fitted on, so `newdata` is
refused for them.

## Argument spellings

frmtmb answers to two dialects, and this family sits on the seam. The
rule is that a brms-NAMED function speaks brms's argument names, while
frmtmb's own fit surface
([`frmtmb::predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.html),
[`frmtmb::simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.html),
[`frmtmb::frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.html))
keeps lme4's, because that is the heritage each name comes from and a
reader should be able to tell which library a call was written against.

`posterior_epred()` and its relatives are brms functions, so the
random-effect switch is `re_formula`. They also SHIPPED taking lme4's
`re.form`, so that spelling keeps working and means exactly the same
thing: both names feed one internal setting, and whichever one is given
wins. brms does the same on `posterior_epred.brmsfit()`, which carries
`re_formula` and `re.form` side by side.

Giving both at once is refused rather than resolved. Two names for one
setting supplied together is a question about what was meant, and
guessing at it would silently ignore one of them.

The literal default of both formals is an internal "not supplied" marker
rather than a value, because `NULL` (keep the random effects) and `NA`
(drop them) are both real settings here and neither can double as
"unset". The behavior when neither is given is unchanged: `NULL` on
every draws method, `NA` on
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.html)
for a fit.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rpois(80, exp(0.3 + 0.4 * dd$x + rnorm(8, 0, 0.5)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

nd <- data.frame(x = c(-1, 0, 1),
                 g = factor(1, levels = levels(dd$g)))

# the expected response per draw: uncertainty in the mean
ep <- posterior_epred(ds, newdata = nd)
apply(ep, 2, quantile, c(0.025, 0.5, 0.975))

# the predictive distribution adds the family's own noise, so its
# intervals are wider
pp <- posterior_predict(ds, newdata = nd)
apply(pp, 2, quantile, c(0.025, 0.5, 0.975))

# the linear predictor itself, on the link scale by default
head(posterior_linpred(ds, newdata = nd, ndraws = 5))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0, 2.5)
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#>             [,1]       [,2]      [,3]
#> [1,] -0.34673055 -0.1204018 0.1059269
#> [2,]  0.17182952  0.3451118 0.5183942
#> [3,] -0.12957334  0.2246253 0.5788240
#> [4,]  0.09809158  0.3436258 0.5891600
#> [5,]  0.07641280  0.3309746 0.5855364
# }
```
