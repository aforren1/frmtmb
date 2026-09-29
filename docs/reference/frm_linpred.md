# Predictions from a frmtmb fit

Predictions from a frmtmb fit

## Usage

``` r
frm_linpred(
  object,
  newdata = NULL,
  type = c("link", "response", "conditional", "zprob", "zlink", "disp"),
  dpar = NULL,
  resp = NULL,
  re_formula = NULL,
  se.fit = FALSE,
  allow_new_levels = FALSE,
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- newdata:

  Optional data frame to predict on. Defaults to the training data.

- type:

  `"link"` for the linear predictor, `"response"` for the expected
  response (which equals
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on the
  training data; for zero-inflated, hurdle, and similar families this is
  the response mean, not the `mu` dpar). When `dpar` is given,
  `"response"` is that dpar on its natural scale. The glmmTMB spellings
  `"conditional"` (the `mu` dpar on its natural scale),
  `"zprob"`/`"zlink"` (the zero-inflation/hurdle probability on the
  response/link scale), and `"disp"` (the dispersion dpar) are accepted
  as aliases.

- dpar:

  Which distributional parameter to predict; defaults to the family's
  first location parameter (`"mu"` for most families).

- resp:

  For multivariate fits: which response to predict (defaults to the
  first).

- re_formula:

  Which group-level terms enter the prediction. `NULL` (default) keeps
  all of them; `NA` keeps none, which is the population-level
  prediction. A one-sided formula keeps the terms it names; see *A
  one-sided `re_formula`*. See *What `re_formula = NA` drops* for what
  that means when the model has smooths.

- se.fit:

  If `TRUE`, return a list with elements `fit` and `se.fit`
  (delta-method standard errors accounting for fixed-effect and
  random-effect uncertainty). Exact `gp()` terms predict unseen
  positions by kriging: the conditional mean at the fitted kernel, with
  the GP conditional variance added to the standard errors.

- allow_new_levels:

  Predict unseen grouping-factor levels at the population level instead
  of erroring. A factor-smooth term (`bs = "fs"`) follows the same rule:
  a level it never saw contributes nothing, which leaves the population
  curve. A smooth with an `re` basis or margin cannot, because its
  design has one column per fitted level and no zero row, so an unseen
  level there is refused by name.

  For a `(x | g)` term it also makes the grouping COLUMN optional, as it
  does in brms (`validate_newdata()`: "grouping factors do not need to
  be specified by the user if new levels are allowed"). A column
  `newdata` does not carry is filled with `NA`, so every row is an
  unseen level. Without `allow_new_levels` a missing grouping column is
  refused by name, and the refusal offers this argument and
  `re_formula = NA`, which drops the random effects instead. A SMOOTH
  term's grouping column is never optional: `re_formula` does not remove
  the term, so the basis has to be rebuilt at a named level.

- ...:

  Refused. An argument this method does not have is an error naming it,
  and the two lme4 spellings that were live in 0.57.0 (`re.form`,
  `allow.new.levels`) are refused by name with the brms spelling that
  replaced them.

## Value

A numeric vector, or a list when `se.fit = TRUE`. For an ordinal family
with `type = "response"`, an `n x K` matrix of category probabilities.

## Details

When the fixed-effect design was rank deficient, the aliased columns
were dropped at fit time and some coefficient combinations are not
estimable. Rows of `newdata` that load on a dropped direction get `NA`
(and `NA` standard errors), with one warning naming the dropped columns;
every other row is unaffected. The test is the one
[`stats::predict.lm()`](https://rdrr.io/r/stats/predict.lm.html) uses: a
row is non-estimable when it is not orthogonal to the null space of the
fitted design, up to a relative tolerance of `1e-8`. Two limits follow.
It is a numerical test, so near-aliased designs sit on a threshold
rather than a clean yes/no. And it covers the parametric fixed-effect
block only: smooth null-space, `gp()`, `mo()` and `mi()` columns are
appended after the rank check and are never dropped.

## Truncated responses

For a response with [`trunc()`](https://rdrr.io/r/base/Round.html)
bounds, `type = "response"` (and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)) report the
truncated mean `E[Y | lb <= Y <= ub]`, matching the likelihood the model
was fitted with. Predictions of a distributional parameter
(`type = "link"`, `dpar = `, or `type = "conditional"`) stay
**untruncated**: they are statements about the latent parameter, not
about the observed, truncated response. Bounds are re-evaluated on
`newdata` the same way `trials()` and `se()` are: a literal bound
carries over unchanged, and a bound given as a variable must be a column
of `newdata` of the right length.

## Ordinal responses

[`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
and
[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
have no mean on the response scale, so `type = "response"` (and its
alias `type = "conditional"`) returns an `n x K` matrix of category
probabilities instead of a vector - the brms
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) convention -
with the response's own factor levels as column names. The rows sum to
one. `cs()` category-specific terms are honored: they enter each
threshold separately and are re-evaluated on `newdata`.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
same matrix, so the usual `frm_linpred(type = "response") == fitted()`
identity holds here too.

`type = "link"` (the default) and `dpar = "mu"` still give the latent
linear predictor, which is where the fixed-effect coefficients live and
where `se.fit` is available. `se.fit` on the response scale is refused:
the prediction is a K-vector per row, not one number. `emmeans` and
[`insight::get_predicted()`](https://easystats.github.io/insight/reference/get_predicted.html)
stay on that latent scale, which is the `mode = "latent"` convention for
`clm`-like models.

`type = "conditional"` is glmmTMB's name for the conditional MEAN, so it
gives the category probabilities here too rather than the linear
predictor: an ordinal response has no mean, and answering a question
about a mean with a latent predictor is the confusion this section
exists to remove. Ask for the predictor by name (`type = "link"`, or
`dpar = "mu"`) when that is what you want.

## A one-sided `re_formula`

A formula keeps the group-level terms it names and drops the others,
with brms's rule (`update_re_terms()`). A term in the formula is matched
to a term of the fit that has the same grouping factor and whose columns
include the formula term's columns. So on a fit with
`(1 + x | g) + (1 | h)`:

- `~ (1 | g)` keeps the intercept of `g` and drops its slope and `h`.

- `~ (0 + x | g)` keeps the slope of `g` alone.

- `~ (1 + x | g) + (1 | h)` keeps everything, and the answer is
  identical to `re_formula = NULL`.

- `~0` and `~1` name no group-level term, so they are `NA`.

The same term is kept in every distributional parameter that has it.
brms's spellings are read too: an id (`(1 | p | g)`), a `gr()` wrapper,
`||` and a nested group (`a/b`). A term outside the bars is ignored.

A dropped term is dropped everywhere: from the estimate, from the
standard error, and from
[`predict()`](https://rdrr.io/r/stats/predict.html)'s draws. Its
grouping column is not needed in `newdata`, and a level of it the fit
never saw is not an error.

One case is refused: a formula term that matches no term of the fit is
an error that names it, because a misspelled grouping factor would
otherwise change the answer with nothing said; brms drops such a term
silently. A `car()` or `spde()` field is a group-level term named by its
grouping factor, `(1 | loc)`, so a formula keeps or drops it like any
other. A smooth is not a group-level term: it stays at every
`re_formula`, a partial one included, as in brms.

## What `re_formula = NA` drops

`re_formula = NA` (equivalently `~0` or `~1`) asks for the
POPULATION-level prediction. It removes the `(x | g)` group-level terms
and nothing else. Every SMOOTH stays in.

Dropped:

- `(1 | g)`, `(x | g)`, and the structured spellings of them (`gr()`,
  `cs()`, [`ar()`](https://rdrr.io/r/stats/ar.html), `mm()`, `car()`,
  `spde()`, ...).

Kept:

- `s(t)`, `s(t, by = x)`, `t2()` and every other smooth. A smooth's
  wiggly part is stored as a random-effect block because that is how a
  penalty is written as a mixed model, but the term is part of the
  formula and the population prediction is the fitted curve, not the
  null-space line through it.

- `s(t, g, bs = "fs")`, `s(g, bs = "re")`, `s(x, g, bs = "re")` and
  `t2(t, g, bs = c("cr", "re"))`: a smooth whose basis gives each level
  of a grouping factor its own curve is still a smooth. It is kept, at
  each row's own level.

- `gp()` and `hsgp()` terms.

This is brms's rule.
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
at `re_formula = NA` is BITWISE the same as at `re_formula = NULL` on a
fit whose only group-indexed content is `s(g, bs = "re")`,
`s(x, g, bs = "fs")` or a `t2()` with an `re` margin, and differs on
`s(x) + (1 | g)` (measured: `dev/resmooth-brms.txt`). Through 0.64.0
frmtmb dropped those three, which answered a different model: on
`y ~ s(x, g, bs = "fs")` the whole fitted structure went, so
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
drew a FLAT line and the default
[`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
compared the data against draws whose per-row spread was 2.93 times the
fitted sigma.

The test is what the basis MEANS, not the `bs` string, because the
grouping factor still has to be named for `newdata` and still may not be
a level the fit did not see.

In `mgcv`'s spelling the answer is `predict.gam()` with the group-level
intercept excluded and every smooth left in, which
`tests/testthat/test-smooth-population.R` asserts (the two packages'
independent fits agree to about 8e-7 relative).

A kept factor-smooth term reads its grouping column from `newdata` at
every `re_formula`, and `re_formula = NA` is not a way around a missing
column or an unseen level. brms refuses both too, and not through its
group-level machinery: the grouping factor of a smooth is an ordinary
predictor there, so a missing column is a missing variable and a new
level is "New factor levels are not allowed" whatever `allow_new_levels`
says. frmtmb refuses by name, and keeps one opt-in brms does not have:
under `allow_new_levels = TRUE` an unseen level of an `fs` term takes
the population curve, because mgcv's `fs` basis returns a zero row for a
level it does not know. An `re` basis or margin has one column per
fitted level and no such row, so an unseen level there is refused either
way.

Because the basis is now rebuilt at every `re_formula`, one fittable
term is left with no `newdata` route at all: a factor smooth that also
carries a `by =` factor, `s(x, g, bs = "fs", by = f)`. mgcv's
random-effect split of that basis is one this version cannot invert, so
`predict(newdata = )`, `fitted(newdata = )`, `frm_linpred(newdata = )`
and
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
all stop with that named error at every `re_formula`, where
`re_formula = NA` used to answer by dropping the term. What it answered
was the intercept at every row, so the refusal replaces a flat line and
not a curve. In-sample
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`frm_linpred()` on such a fit work as they always did.

[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
passes `re_formula = NA`, so it follows this rule too: on a model with a
factor-smooth term the displayed curve is the curve of the grouping
factor's REFERENCE level, which is where the display holds every
predictor it is not varying, and the grouping factor is not offered as
an effect to plot.

## Standard errors of the expected response

For a family whose mean is the `mu` dpar, `se.fit` on A dpar whose
response scale is not its own link inverse, such as a
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
mixing weight reporting the softmax, takes the delta method through that
transform with respect to its OWN predictor. That is exact for a
two-component mixture and conservative for three or more; see
[`?mixture`](https://aforren1.github.io/frmtmb/reference/mixture.md).

`type = "response"` is the usual one-predictor delta method:
`|dmu/deta| * se(eta)`.

When the mean is a function of several dpars (zero-inflated and hurdle
families, `lognormal`, a `trials()` binomial, or any
[`trunc()`](https://rdrr.io/r/base/Round.html)ed response), the delta
method runs jointly over every dpar's linear predictor: `se^2 = g' V g`,
where row `i` of `g` stacks `dm_i/deta_k` times the design row of
predictor `k`, and `V` is the joint covariance of all the coefficients
([`vcov()`](https://rdrr.io/r/stats/vcov.html)'s `jointPrecision` block,
so the cross-predictor covariances and the shared random-effect block
are included). The gradients `dm/deta_k` are central differences of the
family mean, taken one predictor at a time with a relative step.

The estimate is at the random-effect modes, and the standard error
carries their uncertainty, through the joint covariance, as `se.fit`
does for the linear predictor. Unseen grouping levels
(`allow_new_levels = TRUE`) add their block's marginal variance,
propagated through the same gradients.

## See also

[`fitted.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/fitted.frmtmb_fit.md)
for the same expected response in brms's four-column shape,
[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
for the predictive summary, and
[frmtmb-scales](https://aforren1.github.io/frmtmb/reference/frmtmb-scales.md),
which states which scale every method reports. The default here is the
LINK scale.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

# the link scale by default; "response" is what fitted() estimates
head(frm_linpred(fit))
#>            1            2            3            4            5            6 
#> -0.705162793 -0.075119938 -0.231249581  0.978192337 -0.003131653  0.940684161 
max(abs(frm_linpred(fit, type = "response") -
          fitted(fit)[, "Estimate"]))
#> [1] 0

# re_formula = NA drops the random effects: the population prediction
nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
frm_linpred(fit, newdata = nd, re_formula = NA, type = "response")
#>         1         2         3 
#> 0.9340309 1.3398092 1.9218728 

# delta-method standard errors, on whichever scale was asked for
p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
cbind(fit = p$fit, se = p$se.fit)
#>          fit        se
#> 1 -0.8399281 0.3626996
#> 2 -0.4791552 0.3357936
#> 3 -0.1183823 0.3299798

# a level the fit never saw errors unless it is allowed explicitly,
# in which case it is predicted at the population level
nd_new <- data.frame(x = 0, g = factor("new"))
try(frm_linpred(fit, newdata = nd_new))
#> Error : New levels in grouping factor `g`: new. Use allow_new_levels = TRUE to predict them at the population level
frm_linpred(fit, newdata = nd_new, allow_new_levels = TRUE)
#>         1 
#> 0.2925272 

# a distributional parameter instead of the mean
fit2 <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)
head(frm_linpred(fit2, dpar = "sigma", type = "response"))
#>        1        2        3        4        5        6 
#> 1.511247 1.889479 1.426552 2.788566 1.967020 1.432528 
```
