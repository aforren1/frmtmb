# Define a model family

Constructs a family object for use with
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md). The
log-density function must be vectorized and AD-compatible: it is
evaluated on RTMB 'advector' objects during taping, so it must use
RTMB-overloaded operations (RTMB and RTMBdist `d*` functions, plain
arithmetic) and must not branch on parameter values.

## Usage

``` r
frmtmb_family(
  family,
  dpars,
  links,
  lpdf,
  valid_y = NULL,
  init_dpars = list(),
  type = "continuous",
  post = list(),
  sim = NULL,
  sim_ctx = NULL,
  sim_refusal = NULL,
  primary_dpars = "mu",
  lcdf = NULL,
  extra_pars = NULL,
  drop_intercept = FALSE,
  structure = NULL
)

custom_family(
  family,
  dpars,
  links,
  lpdf,
  valid_y = NULL,
  init_dpars = list(),
  type = "continuous",
  post = list(),
  sim = NULL,
  sim_ctx = NULL,
  sim_refusal = NULL,
  primary_dpars = "mu",
  lcdf = NULL,
  extra_pars = NULL,
  drop_intercept = FALSE,
  structure = NULL
)
```

## Arguments

- family:

  Character name of the family.

- dpars:

  Character vector of distributional parameter names. The first entry
  must be `"mu"`.

- links:

  Named list mapping each dpar to a link name (see
  `frmtmb:::frmtmb_links`) or a link object.

- lpdf:

  Function `(y, dpars, aterms)` returning the vectorized log-density.
  `dpars` is a named list of advector vectors; `aterms` is a named list
  of numeric addition-term values (for example `trials`).

- valid_y:

  Optional function `(y, aterms)` that signals an error for invalid
  responses. Called once at assembly time.

- init_dpars:

  Optional named list of functions `(y, aterms)` giving a response-scale
  starting value per dpar (applied to the intercept through the link).

- type:

  One of `"continuous"`, `"discrete"`, `"ordinal"`, `"categorical"`. The
  last two say the modelled response is a distribution over `1..K`
  categories rather than a number, so
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns an
  `n x K` probability matrix; `"ordinal"` shares one latent predictor
  across the categories and `"categorical"` gives each non-reference
  category its own.

- post:

  Named list of numeric helper functions used by
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `mean_fn(dpars, aterms)` (the response mean), `var_fn(dpars, aterms)`
  (for pearson residuals) and `dev_fn(y, dpars, aterms)` (the unit
  deviance
  `2 * (loglik of the saturated fit - loglik at the fitted value)`, for
  `residuals(type = "deviance")`). A family that omits one is refused by
  the method that needs it.

- sim:

  Optional numeric simulator `(dpars, aterms, n)` returning `n` response
  draws; used by [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  and
  [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md).
  It is stateless and rowwise: it sees the distributional parameters and
  nothing else. A family whose extra parameters (`extra_pars`) enter the
  draw declares a fourth argument `extra` instead.

- sim_ctx:

  Optional structured simulator `(ctx)` for a family whose draws are not
  rowwise (see Structured simulators). It takes precedence over `sim`.

- sim_refusal:

  Optional one-sentence reason why the family has no simulator, appended
  to the refusal each entry point raises. Use it when the omission is a
  decision rather than a gap.

- primary_dpars:

  Which dpars receive the main model formula (default `"mu"`). Families
  with several location predictors (for example multinomial's
  per-category `mu2`, `mu3`, ...) list them all; these live in the
  `beta` parameter vector and are integrated out under REML.

- lcdf:

  Optional vectorized AD log-safe CDF `(q, dpars, aterms)` returning
  probabilities; enables `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) addition terms.

- extra_pars:

  Optional function `(y, aterms)` returning a named list of numeric
  starting vectors for family-level parameters outside the dpar system
  (for example ordinal thresholds). They join the parameter template
  under their own names and reach `lpdf` as its fourth argument.

- drop_intercept:

  If `TRUE`, the intercept column is removed from the main formula's
  design matrix (ordinal families: thresholds take its place).

- structure:

  Optional
  [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
  for a family whose likelihood does not factorize over rows (a
  group-level
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md),
  a hidden Markov chain). It carries the non-rowwise log-likelihood, the
  frame block that likelihood reads, and the capability flags that say
  which post-fit methods the family can answer. `lpdf` stays required
  even then, for the rowwise contract, and may be a stub that refuses.

## Value

An object of class `frmtmb_family`.

## Structured simulators

Some families cannot draw a response one row at a time: a group-level
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
draws one class per group, an
[`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) walks a
Markov chain per sequence, and a
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
draw needs the class covariances, which are family-level extras rather
than dpars. Those families supply `sim_ctx(ctx)` instead of
`sim(dpars, aterms, n)`.

`ctx` is a list with `fit` (any object carrying `spec`, `frame` and
`estimates` - a fitted model, one posterior draw, or the de novo shim),
`family`, `rspec`, `resp`, `dpars` (the evaluated numeric distributional
parameters), `aterms`, `n`, `extra` (the family-level extra parameters)
and the frame structures `autocor` and `block` (the structured family's
own data; see
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)).
Read its fields with `[[ ]]`.

The same `sim_ctx()` serves \[simulate()\], \[posterior_predict()\] and
\[frm_simulate()\]. Because a structured draw covers whole sequences or
groups, [`trunc()`](https://rdrr.io/r/base/Round.html) rejection and
`newdata` cannot apply to it and are refused.

\[ \]: R:%20 \[simulate()\]: R:simulate() \[posterior_predict()\]:
R:posterior_predict() \[frm_simulate()\]: R:frm_simulate()

## Tape-safe scope

`lpdf` and `lcdf` run with RTMB's tape-safe
[`c()`](https://rdrr.io/r/base/c.html), `[<-` and `diag<-` in scope
automatically (the `"c" <- RTMB::ADoverload("c")` boilerplate is spliced
in unless the function already binds it), so base spellings keep the
automatic-differentiation class. A helper the density CALLS still needs
its own bindings: lexical scope does not travel into other functions.

## Examples

``` r
# a custom family is a plain R log-density over taped parameters
dd <- data.frame(y = rbinom(100, 5, 0.4),
                 size = 5, x = rnorm(100))
fam <- custom_family(
  "vbinom", dpars = "mu", links = list(mu = "logit"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
  },
  type = "discrete"
)
fit <- frm(bf(y | vint(size) ~ x) + fam, data = dd)
fixef(fit)
#> $mu
#> (Intercept)           x 
#> -0.58378094 -0.02406052 
#> 
```
