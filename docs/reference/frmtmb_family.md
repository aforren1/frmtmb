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
  lccdf = NULL,
  required_aterms = character(0),
  accepts_aterms = NULL,
  se_dpar = NULL,
  exclusive_aterms = list(),
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
  lccdf = NULL,
  required_aterms = character(0),
  accepts_aterms = NULL,
  se_dpar = NULL,
  exclusive_aterms = list(),
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

  Named list mapping each dpar to a link name or a link object.
  [frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md)
  lists the names and says what a link object has to carry.

- lpdf:

  Function `(y, dpars, aterms)` returning the vectorized log-density.
  `dpars` is a named list of advector vectors; `aterms` is a named list
  of numeric addition-term values (for example `trials`). Read `dpars`
  by name, never by position: it also carries reserved entries that are
  not distributional parameters. Take `log(mu)`, `log(1 - mu)` and
  `1 - mu` from
  [frmtmb-robust-dpars](https://aforren1.github.io/frmtmb/reference/frmtmb-robust-dpars.md)
  rather than writing them out, because an inverse link saturates and
  the plain arithmetic returns `NaN` for the value and the gradient
  alike in the tail.

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

  `post$fit_check(fit, resp)` is different in kind: it is run once, when
  a fit FINISHES, and its return value is discarded. It is where a
  family says something about where the optimizer landed, which nothing
  else can: [`logLik()`](https://rdrr.io/r/stats/logLik.html) reads the
  optimizer's own value, so a family whose likelihood is floored or
  degenerate in some region had no way to report it. Warn from it rather
  than stopping; a hook that throws is caught, reported as a warning
  naming the family, and the fit is returned regardless.

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
  [`trunc()`](https://rdrr.io/r/base/Round.html) addition terms. It is
  the only thing either term asks of a family, whatever the family's
  `type`: a discrete family that supplies one is censored under the
  inclusive convention (see Censoring a discrete response), which reads
  a lower bound as `F(q - 1)` and so calls this at one below a recorded
  value.

- lccdf:

  Optional vectorized AD LOG SURVIVOR function `(q, dpars, aterms)`
  returning `log(1 - F(q))` directly. A family that declares it scores a
  RIGHT-censored row from it instead of from `log(1 - F)`, which cannot
  be accurate once `F` rounds to one (see Right censoring and the
  representable tail). It is optional and independent of `lcdf`: a
  family that supplies only `lccdf` accepts right censoring and refuses
  left censoring, interval censoring and
  [`trunc()`](https://rdrr.io/r/base/Round.html), each by name.

- required_aterms:

  The addition-term values the density cannot do without, named as they
  reach `aterms`: `"vint1"`, `"vreal2"`, `"trials"`. A character vector
  names the terms it needs ALL of. A LIST adds the alternative: an
  element of length one is a term that must be there, and an element of
  length more than one is a set of spellings, ANY one of which will do,
  so `required_aterms = list(c("dec", "vint1"), "vreal1")` reads as one
  of `dec` or `vint1`, and `vreal1`. A family that reads the same datum
  from either of two terms declares it that way instead of hand-rolling
  the refusal.

  Frame assembly refuses a model that leaves a requirement unmet, naming
  the family, the terms and the spelling that supplies them; for a set
  of alternatives it names all of them and writes the first into the
  example formula. Without the declaration an absent term reaches the
  density as `NULL`, arithmetic on it gives `numeric(0)`, and the
  log-likelihood becomes a sum over nothing: a fit that returns, with a
  log-likelihood of zero. Declare every per-row datum the density
  indexes.

- accepts_aterms:

  The addition terms this family reads or lets the core act on, named as
  a formula writes them and without parentheses:
  `c("weights", "trials", "cens")`. Frame assembly refuses any other
  term on the response, by name, and lists the ones the family takes.
  `character(0)` declares a family that takes none. `NULL`, the default,
  accepts every registered term, which is what a family written before
  this argument existed keeps.

  `required_aterms` is a conjunction of what the density cannot do
  without; this is the complementary allow-list, and the two are read
  together, so a required term need not be repeated here. Without a
  declaration an unread term is parsed, stored on the fit and silently
  ignored: `wiener()` accepted a `vint()` it cannot use, and `lba()`
  accepted a `dec()`, both giving a fit bit-identical to the one without
  the term.

  Naming `"se"` here is more than an allow-list entry: it is how a
  family OPTS IN to `se()`. That term is the one whose entire effect is
  inside the density - the core hands over `aterms[["se"]]` and maps out
  the residual scale it replaces, and does nothing else with it - so a
  family that does not declare it is refused the term rather than given
  one it would ignore. `NULL`, which accepts every other term, declares
  nothing and so does not open this one. Read `aterms[["se"]]` as the
  known standard deviation, and `aterms[["se_sigma"]]` to honor
  `se(x, sigma = TRUE)`, which asks for the known and estimated scales
  in quadrature.

- se_dpar:

  The dpar that a known standard error replaces, named so that the core
  can map it out: `se_dpar = "tau"` maps out `tau` exactly as the
  convention maps out `sigma`. `NA` declares that `se()` replaces NO
  dpar, which is the shape of a family whose whole scale IS the known
  standard error. `NULL`, the default, reads the convention: the dpar
  named `sigma`, if the family has one.

  Only a family that declares `se()` reaches this. `se()` without
  `sigma = TRUE` says the residual scale is known, so the dpar it
  replaces has to stop being estimated. A dpar the density never reads
  is a flat direction and a NaN standard error, so a declaring family
  with no `sigma` and another free dpar is refused: the core cannot tell
  a second SCALE from a genuine SHAPE. This argument is how the family
  says which one it has. `se_dpar = NA` is a promise that every
  remaining dpar is read alongside the known standard error, the way a
  skew or a tail index is.

- exclusive_aterms:

  Sets of addition-term values that say the SAME thing to the density,
  so that at most one of each set may be supplied. Named as
  `required_aterms` names them, in values rather than in terms:
  `exclusive_aterms = list(c("dec", "vint1"))`. A character vector is
  one such set; a list is several. The FIRST name of a set is the
  spelling the density reads, and the refusal tells the user to keep it.

  This is what an allow-list cannot express, because both spellings are
  legitimately on it. `wiener()` reads its boundary indicator from
  `dec()` and falls back to `vint1`, so `rt | dec(u) + vint(1 - u)`
  passed every guard and fitted with a log-likelihood bit-identical to
  the `dec()`-only model: the user said one thing twice, and
  inconsistently, and nothing complained. Declare the alternatives of an
  any-of `required_aterms` group here when the density reads only one of
  them; leave them undeclared when it reads both, as `gddm()` does,
  where `dec()` and `vint1` carry different data.

  Read WITH `required_aterms`, not instead of it: an any-of group in one
  and the same set in the other together mean "exactly one". Declaring a
  set that holds two values `required_aterms` demands TOGETHER is
  refused here, because no model could then be fitted.

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

## Right censoring and the representable tail

Core forms a right-censored row's contribution as `log(Fub - F(y))`,
which without truncation is `log(1 - F)`. A double cannot represent the
complement of a probability that has rounded to one, so past that point
the contribution is not merely inaccurate: it is CONSTANT, and its
gradient is exactly zero. An optimizer then prices such a row the same
however far it moves, and fits every other row as if the survivor were
free. That failure is silent - the fit converges, and
[`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) report the floored number.

`lccdf` removes the class for right censoring by giving core the
quantity it actually needs. Measured on a standard normal tail, one
process:

|       |                       |                                                |
|-------|-----------------------|------------------------------------------------|
| **z** | **log(1 - pnorm(z))** | **pnorm(z, lower.tail = FALSE, log.p = TRUE)** |
| 8.0   | -35.013               | -35.013                                        |
| 8.3   | -Inf                  | -37.494                                        |
| 37    | -Inf                  | -689.031                                       |
| 500   | -Inf                  | -125007.13                                     |

Both the value and the derivative are exact on the right, over the whole
range.

The built-in families that declare `lccdf` are
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`lognormal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`exponential()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`weibull()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
and
[`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).
The last three have `log S` in closed form (`-q/mu`, `-(q/scale)^shape`,
`-H0(t) * mu`); the first two use
[`pnorm()`](https://rdrr.io/r/stats/Normal.html)'s own log upper tail.

Two families with a CDF do NOT declare one, for measured reasons.
[`inverse.gaussian()`](https://rdrr.io/r/stats/family.html) gains
nothing: `RTMBdist::pinvgauss(lower.tail = FALSE, log.p = TRUE)` is
computed on the probability scale and reaches `-Inf` at the same
`log S = -34` that `log(1 - F)` does.
[`poisson()`](https://rdrr.io/r/stats/family.html) is censored (see
Censoring a discrete response) and would use the slot, but cannot
declare one: `RTMB::ppois(lower.tail = FALSE, log.p = TRUE)` is exact in
R (`-2773.28` at `q = 700`, `lambda = 5`, where `log(1 - F)` is `-Inf`)
and does not TAPE - inside `MakeADFun` it reaches
[`stats::ppois`](https://rdrr.io/r/stats/Poisson.html) and errors with
"Non-numeric argument to mathematical function". An exact discrete log
survivor has to be written out before poisson can have one.

`lccdf` fixes RIGHT censoring and nothing else. Left censoring is still
`log(F(y) - Flb)`, interval censoring is still a difference of CDFs, and
the truncation normalizer is still `log(Fub - Flb)`, so a LEFT-TRUNCATED
survival model - delayed entry, which is routine - meets the identical
representability problem from the other side. Closing that needs a
windowed log-difference slot, and this is the first step rather than the
last one.

## Censoring a discrete response

A censoring bound on a discrete response NAMES a value the response can
take, and the value is INCLUDED in the event:

|                |                |                    |
|----------------|----------------|--------------------|
| **code**       | **means**      | **scored as**      |
| `0` "none"     | `Y == y`       | `f(y)`             |
| `-1` "left"    | `Y <= y`       | `F(y)`             |
| `1` "right"    | `Y >= y`       | `1 - F(y - 1)`     |
| `2` "interval" | `y <= Y <= y2` | `F(y2) - F(y - 1)` |

Equivalently: every LOWER edge enters the CDF as `F(edge - 1)`, and
upper edges are unchanged because `F` already includes its argument. It
is the rule `trunc(lb = )` has always followed on a discrete response,
so one number means one thing however a response is bounded, and it is
what a count recorded as "5 or more" means.

It DIFFERS from brms for RIGHT and INTERVAL censoring only, where brms
emits `poisson_lccdf(y | mu)`, that is `P(Y > y)`, and reads an interval
as `(y, y2]`. LEFT censoring is `P(Y <= y)` in both packages and agrees
exactly. Migrating a right- or interval-censored count model changes its
log-likelihood; on 200 poisson draws at `lambda = 4` right censored at
6, the two readings differ by 20.8 log units and 3.7 percent of the
estimate. Subtract one from every right-censored and interval lower
bound to reproduce a brms fit.

The divergence is deliberate, because **brms is internally inconsistent
here and frmtmb cannot be both.** brms's own discrete truncation emits
`poisson_lccdf(lb - 1 | mu)`, an INCLUSIVE lower bound, `P(Y >= lb)`. So
in brms `trunc(lb = 6)` means `Y >= 6` while a right-censored row
recorded at 6 means `Y > 6`: one number, two meanings, on one response.
frmtmb's [`trunc()`](https://rdrr.io/r/base/Round.html) reproduces brms
bit for bit (poisson on `y = 3,4,5,6` at `b0 = log 4`, `trunc(lb = 2)`
gives 6.9990718955 under both, against 6.2954805379 for the exclusive
reading), so its `cens()` had to choose between matching brms's
censoring and matching its own truncation. It matches its own.

It is also the only reading consistent with the censored SIMULATOR,
which predates all of this: `simulate(censored = TRUE)` caps a draw with
`pmin(pmax(y, lo), hi)`, recording the value `k` exactly when the latent
draw is `>= k`. Measured on 4000 draws from a fit censored at 7, the
simulated mass at that point is 0.13250, against `P(Y >= 7) = 0.12190`
inclusive and `P(Y > 7) = 0.05763` exclusive. The other reading would
silently decouple the likelihood from the simulator that
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
rests on.

Two consequences worth knowing. A one-point interval (`y2 == y`) is
legal and is the exact observation `P(Y = y)`, where on a continuous
response it is refused as an event of probability zero. And the shift
assumes the support is the unit integer lattice, so a non-integer
censoring bound is refused rather than moved onto a point the family has
no mass at.

`residuals(type = "osa")` is refused on a censored discrete fit:
inclusive bounds make an uncensored row's support `[lo + 1, hi - 1]`
rather than the `[lo, hi]` the one-step window is built on.

## Tape-safe scope

`lpdf` and `lcdf` run with RTMB's tape-safe
[`c()`](https://rdrr.io/r/base/c.html), `[<-` and `diag<-` in scope
automatically (the `"c" <- RTMB::ADoverload("c")` boilerplate is spliced
in unless the function already binds it), so base spellings keep the
automatic-differentiation class. A helper the density CALLS still needs
its own bindings: lexical scope does not travel into other functions.

## See also

[frmtmb-robust-dpars](https://aforren1.github.io/frmtmb/reference/frmtmb-robust-dpars.md)
for the accessors a density uses to stay exact where an inverse link
saturates,
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
for the accessors a family outside frmtmb may use after a fit

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
