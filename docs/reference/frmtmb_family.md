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
  required_aterms = character(0),
  family_finalize = NULL,
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
  required_aterms = character(0),
  family_finalize = NULL,
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
  `posterior_predict()` and
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

- required_aterms:

  Character vector of addition-term values the density cannot do
  without, named as they reach `aterms`: `"vint1"`, `"vreal2"`,
  `"trials"`. Frame assembly refuses a model that omits one, naming the
  family, the term and the spelling that supplies it. Without the
  declaration an absent term reaches the density as `NULL`, arithmetic
  on it gives `numeric(0)`, and the log-likelihood becomes a sum over
  nothing: a fit that returns, with a log-likelihood of zero. Declare
  every per-row datum the density indexes.

- family_finalize:

  Optional function `(fam, y, aterms)` returning a family. It runs once
  at frame assembly, after the response is coerced and validated and
  before any link is used, and whatever it returns is the family the
  rest of the fit sees. It is how a family derives a link bound, a
  default, or an extra slot FROM the data instead of asking the user for
  a quantity the framework already holds (see Deriving a family from the
  data).

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
draws one class per group, a
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
draw needs the class covariances, which are family-level extras rather
than dpars, and a hidden Markov family walks a chain per sequence. Those
families supply `sim_ctx(ctx)` instead of `sim(dpars, aterms, n)`.

`ctx` is a list with `fit` (any object carrying `spec`, `frame` and
`estimates` - a fitted model, one posterior draw, or the de novo shim),
`family`, `rspec`, `resp`, `dpars` (the evaluated numeric distributional
parameters), `aterms`, `n`, `extra` (the family-level extra parameters)
and the frame structures `autocor` and `block` (the structured family's
own data; see
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)).
Read its fields with `[[ ]]`.

The same `sim_ctx()` serves \[simulate()\], `posterior_predict()` and
\[frm_simulate()\]. Because a structured draw covers whole sequences or
groups, [`trunc()`](https://rdrr.io/r/base/Round.html) rejection and
`newdata` cannot apply to it and are refused.

\[ \]: R:%20 \[simulate()\]: R:simulate() \[frm_simulate()\]:
R:frm_simulate()

## Slot call order

The order the slots run in is part of the contract, because a family
that derives anything from the data depends on it. Measured on an
instrumented family, one fit of one response:

1.  `valid_y(y, aterms)`, once, at frame assembly, with the response
    coerced and the addition terms evaluated.

2.  `family_finalize(fam, y, aterms)`, once, immediately after, and
    still before any link function is called.

3.  `aterm_data()` and `extra_pars(y, aterms)`, once each, at assembly.

4.  `init_dpars[[dpar]](y, aterms)`, once per dpar, when the starting
    values are built; each value goes straight through that dpar's
    `linkfun`.

5.  `linkinv` and `lpdf`, on the tape, from then on.

`post$mean_fn`, `post$var_fn`, `post$dev_fn` and `sim` are never called
by [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md). They
run only when a post-fit method asks for them.

So a link may read anything `valid_y` or `family_finalize` computed, and
neither of those may read anything a link produced.

## Deriving a family from the data

Some families are not fully determined until the response is in hand. A
shifted family whose density is zero below a non-decision time wants a
link bounded above by `min(y)`, so that the constraint is structural
rather than left to the optimizer. The bound is a property of the data,
and the family object is built before
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) sees any.

`family_finalize` closes that gap. It receives the family, the coerced
response and the evaluated addition terms, and returns the family the
fit will use:

    shifted <- frmtmb_family(
      "shifted", dpars = c("mu", "ndt"),
      links = list(mu = "log", ndt = "log"),
      lpdf = function(y, dpars, aterms) { ... },
      family_finalize = function(fam, y, aterms) {
        # a logit onto (0, min(y)): ndt cannot reach the smallest
        # observation, whatever the optimizer tries
        ub <- min(y)
        fam$links$ndt <- list(
          name    = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
          linkfun = function(mu) log(mu / (ub - mu)),
          linkinv = function(eta) ub / (1 + exp(-eta)),
          mu_eta  = function(eta) {
            p <- 1 / (1 + exp(-eta))
            ub * p * (1 - p)
          }
        )
        fam
      }
    )

Replacing a link this way replaces it everywhere: the starting values,
the tape, and every post-fit method read the finalized family. The
alternative an extension reached for before this slot existed was to
have `valid_y` write the bound into an environment the link closures
read at call time, which works only for as long as the call order
happens to hold and leaves the family object lying about what it is.

## Tape-safe scope

`lpdf` and `lcdf` run with RTMB's tape-safe
[`c()`](https://rdrr.io/r/base/c.html), `[<-` and `diag<-` in scope
automatically (the `"c" <- RTMB::ADoverload("c")` boilerplate is spliced
in unless the function already binds it), so base spellings keep the
automatic-differentiation class. A helper the density CALLS still needs
its own bindings: lexical scope does not travel into other functions.

## See also

[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for a likelihood that does not factorize over rows,
[`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.md)
for giving the family's per-row data a name of its own instead of
`vint()`,
[`frmtmb_register_compat()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.md)
for telling
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
what the family does and does not combine with, and
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the accessors a family outside frmtmb may use

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
  # the density indexes vint1, so a model without it is refused
  # rather than fitted against a zero-length log-likelihood
  required_aterms = "vint1",
  type = "discrete"
)
fit <- frm(bf(y | vint(size) ~ x) + fam, data = dd)
fixef(fit)
#> $mu
#> (Intercept)           x 
#> -0.58378094 -0.02406052 
#> 
```
