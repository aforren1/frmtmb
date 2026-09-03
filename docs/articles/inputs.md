# Inputs and preprocessing

This is a reference page. It says what
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) accepts,
what it does to what you give it, what reaches the output, and what it
costs to run. For a tour of the modeling grammar, see
[`vignette("frmtmb")`](https://aforren1.github.io/frmtmb/articles/frmtmb.md).

## Terminology

The package uses a small vocabulary consistently. These are the terms.

**Response.** The left side of a model formula, plus its addition terms
([`weights()`](https://rdrr.io/r/stats/weights.html), `trials()`,
`cens()`, [`trunc()`](https://rdrr.io/r/base/Round.html), `se()`,
`rate()`). A multivariate model has more than one response.

**Distributional parameter (dpar).** Every parameter of the response
distribution, not only the mean. A Gaussian response has `mu` and
`sigma`; a negative binomial has `mu` and `shape`. Auxiliary parameters
are dpars like any other: each has a link, each takes a formula with the
full predictor grammar, and each can be fixed to a constant with
`bf(..., sigma = 1)`.

**Linear predictor.** One response paired with one dpar gives one linear
predictor, keyed `<response>.<dpar>`. Its value is
`eta = X %*% beta + Z %*% b + offset`, and the dpar is
`link_inverse(eta)`. A model has as many linear predictors as it has
(response, dpar) pairs.

**Link.** The monotone function that maps a dpar onto the unbounded
scale the linear predictor lives on. Links come from an AD-safe
registry, not from
[`stats::make.link()`](https://rdrr.io/r/stats/make.link.html), because
the C level clamping in the base links is invisible to automatic
differentiation.

**Outer and inner parameters.** The optimizer moves the *outer*
parameters: the fixed-effect coefficients (`beta`, `betad`), the
covariance and smoothing parameters (`theta`), and any family extras
such as ordinal thresholds. The *inner* parameters are the random-effect
coefficients `b`, which include the wiggly parts of smooths and the
values of Gaussian process fields. They are not optimized in the outer
loop; they are integrated out.

**Laplace approximation.** The method that integrates the inner
parameters out of the joint likelihood to leave a marginal likelihood in
the outer parameters alone. It comes from TMB and RTMB.
[`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md)
audits it against NUTS draws on the same objective.

**ML and REML.** Under maximum likelihood, only `b` is inner.
`REML = TRUE` moves the fixed effects into the inner set as well, so
they are integrated out too. This is why bounds and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) blocks that name fixed
effects are refused under REML: those coefficients are no longer outer
parameters.

**Covariance structure (covstruct).** The parameterization of one
random-effect block: `us`, `diag`, `cs`, `ar1`, `ou`, `exp`, `gau`,
`mat`, `rr`, `car`, `spde`, and the rest. It decides how many `theta`
entries the block owns and how they map to a covariance matrix.

## How a formula becomes design matrices

Model building has four stages. Each one has an inspectable output, and
`dry_run =` stops the pipeline at the first two.

    bf(y ~ x + (1 | g), sigma ~ s(z)) + gaussian()    user grammar
          |  bf(), mvbf(), +.frmtmb_formula
          v
    frmtmb_formula                                    unevaluated spec
          |  parse_spec()                             [R/parse.R]
          v
    frmtmb_spec                                       data-free description
          |  assemble_frame(spec, data)               [R/frame.R]
          v
    frmtmb_frame                                      X, Z, blocks, parameters
          |  build_objective(frame)                   [R/objective.R]
          v
    nll(pars) closure                                 data baked in as constants
          |  RTMB::MakeADFun(nll, template, random = "b")
          v
    obj  ->  nlminb(obj$par, obj$fn, obj$gr)  ->  sdreport(obj)

**Stage 1, parse.** `parse_spec()` reads the formula alone. It never
touches the data. It splits each linear predictor into a parametric part
(bars and special terms removed), a list of random-effect terms, a list
of smooth terms, and the special terms `mo()`, `mi()`, `cs()`, `gp()`,
`car()` and `spde()`. Because the output is data free, a spec is a pure
description of the model.

**Stage 2, assemble.** `assemble_frame()` binds the spec to a data frame
and produces the matrices.

- The parametric part goes through
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) and
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  The resulting `terms`, `xlevels` and `contrasts` are stored per linear
  predictor, and prediction on new data reapplies exactly those.
  `frmtmb_control(sparse_x = TRUE)` swaps in
  [`Matrix::sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html);
  the two paths are checked against each other so the set of dropped
  aliased columns cannot differ.
- Random-effect terms go through
  [`reformulas::mkReTrms()`](https://rdrr.io/pkg/reformulas/man/mkReTrms.html),
  which returns lme4’s transposed `Zt`. frmtmb transposes each term’s
  slice back out and stores `Z` with observations in rows.
- Smooths go through `mgcv::smoothCon(absorb.cons = TRUE, modCon = 3)`
  and `mgcv::smooth2random(type = 2)`. The fixed part is appended to
  `X`; the wiggly part becomes a one-variance block in `Z`, so the
  smoothing parameter is estimated as a variance component. Prediction
  uses
  [`mgcv::PredictMat()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html)
  with the stored transforms. `modCon = 3` makes that basis the fitted
  basis: a `t2()` smooth otherwise carries a second, sum-to-zero
  prediction constraint that `PredictMat()` honors and the fit does not.
  It changes no fitted quantity. `s()` and `t2()` are supported, `te()`
  is not.
- Gaussian processes, `car()` and `spde()` build their own sparse
  columns and carry their auxiliary data (distances, adjacency,
  finite-element matrices) on the block.
- [`offset()`](https://rdrr.io/r/stats/offset.html) terms are summed per
  linear predictor and added to `eta` after `X %*% beta + Z %*% b`.
  [`weights()`](https://rdrr.io/r/stats/weights.html) is not a design
  matrix at all: it is an addition term that multiplies the
  per-observation log density.

Every linear predictor’s `Z` spans the whole random-effect coefficient
vector, so a `Z` column index *is* a coefficient index. Columns are
ordered level major and coefficient minor inside a block, at
`(level - 1) * dim + offset + coefficient`. Terms linked with an `|ID|`
key are merged into one block, which is only a matter of where their
columns land.

**Stage 3, generate.** `build_objective()` closes over the assembled
matrices and returns a plain R function of the parameter list. Family,
link and covariance dispatch all resolve before taping, so the closure
stays vectorized and free of branches on parameters.

**Stage 4, tape and optimize.** `RTMB::MakeADFun()` tapes the closure
with `random = "b"` (and `"beta"` under REML),
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html) moves the outer
parameters, and `sdreport()` supplies standard errors.

You can look at the first two stages without fitting:

``` r

set.seed(1)
d <- data.frame(y = rnorm(30), x = rnorm(30), g = factor(rep(1:6, 5)))
rownames(d) <- paste0("obs", seq_len(30))

spec <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d, dry_run = "spec")
names(spec$responses$y$dpars)
#> [1] "mu"    "sigma"
spec$responses$y$dpars$mu$fixed
#> ~x

frame <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d, dry_run = "frame")
names(frame$linpreds)
#> [1] "y.mu"    "y.sigma"
c(X = dim(frame$linpreds[["y.mu"]]$X),
  Z = dim(frame$linpreds[["y.mu"]]$Z))
#> X1 X2 Z1 Z2 
#> 30  2 30  6
lengths(frame$par_template)
#>  beta betad     b theta 
#>     2     1     6     1
```

`beta` holds the coefficients of the mean linear predictors, `betad`
those of the other dpars, `b` the random-effect coefficients, and
`theta` the covariance and smoothing parameters. Here the outer problem
has four parameters and the inner problem six.

## Argument types and lengths

[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) is the
only entry point. There is no matrix-only interface: the model is always
a formula and a data frame.

| Argument | Type and length | On mismatch |
|----|----|----|
| `formula` | One [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) object with a family attached by `+`, or one plain two-sided formula with `family =` given separately. | A formula with no left side, or no family from either source, is an error. |
| `data` | One data frame. Its [`nrow()`](https://rdrr.io/r/base/nrow.html) sets the number of observations. Every variable in every linear predictor must resolve in it. | A missing variable is an error naming the variable and the linear predictor. |
| `family` | A `frmtmb_family`, a [`stats::family`](https://rdrr.io/r/stats/family.html) object or its constructor function, or a length-one family name. Overrides any family carried by `formula`. | An unsupported name errors and lists every supported name. |
| `REML` | Length-one logical. | Moves the fixed effects into the inner set, which makes [`drop1()`](https://rdrr.io/r/stats/add1.html) and fixed-effect bounds refuse. |
| `start` | Named list matching the parameter template. Use `frm(..., dry_run = "frame")$par_template` to see the names and lengths. | An unknown component is an error that points at the template. |
| `control` | The list from [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md). | See below. |
| `se` | Length-one logical, default `FALSE`. Defers `sdreport()`, which is roughly a quarter of the fit time. |  |
| `na.action` | A function, as in [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Default [`stats::na.omit`](https://rdrr.io/r/stats/na.fail.html). | `na.fail` errors on any `NA`. `na.exclude` drops rows for the fit and pads [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`predict()`](https://rdrr.io/r/stats/predict.html) and [`simulate()`](https://rdrr.io/r/stats/simulate.html) back to the input length with `NA` in the original positions. |
| `lower`, `upper` | Named numeric vectors over the outer parameters, on the internal scale, using the names [`confint()`](https://rdrr.io/r/stats/confint.html) prints. | A name that is not an outer parameter is an error. Under `REML = TRUE` or `frmtmb_control(profile = TRUE)` the fixed effects are not outer parameters, so bounds naming them are refused rather than silently applied elsewhere. |
| `prior` | A [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md) specification, or a `brmsprior` object brms built, which is translated row by row. Makes the fit penalized (MAP), so the reported log likelihood is penalized too. Spelled as brms spells it; the `priors` of releases before 0.43 is gone, and a call still using it fails as an unused argument. | [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md) checks the argument count of each distribution. A brms prior class with no faithful frmtmb spelling is refused by name rather than translated into a different density. |
| `quadrature` | Length-one logical. Scalar random-intercept models only. | Any other structure, and `mi()`, [`trunc()`](https://rdrr.io/r/base/Round.html), `REML = TRUE` or `profile = TRUE`, is an error naming the conflict. |
| `data2` | Named list of objects that the formula names and `data` cannot hold: the adjacency matrix of `car()`, the mesh triple of `spde()`, and the matrices of `gr(prec = )`, `gr(cov = )` and `equalto()`. | Anything that is not a named list is an error. A name a formula asks for and does not find is an error naming it. |
| `dry_run` | `NULL`, `"spec"`, `"frame"` or `"objective"`. | Any other value is ignored and the model fits normally. |
| `verbose` | `FALSE`, `TRUE`, or an integer level: 0 silent, 1 stage progress, 2 optimizer trace. An explicit `control$verbose` wins. | Never errors. `NA` and negative values become 0. |

[`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)
takes `optimizer` (`"nlminb"`, `"optim"`, or a function), `optCtrl`,
`restarts`, `grad_tol`, `profile`, `sparse_x`, `autoscale`,
`check_nlev_1`, `check_olre` and `verbose`.

[`predict()`](https://rdrr.io/r/stats/predict.html) takes `newdata` (a
data frame), `type`, `dpar`, `resp`, `re.form`, `se.fit` and
`allow_new_levels`. Note the spelling `re.form`, which follows lme4
rather than brms `re_formula`. Unknown arguments in `...` produce a
warning, not an error, because `allow.new.levels` is accepted as an lme4
and glmmTMB spelling.

## Predictor classes

### Accepted

| Class | Treatment |
|----|----|
| `numeric` | Used directly. |
| `factor` | Contrasts are read from the column and stored, then reapplied to new data. The package never creates a factor for a predictor; that would change the model silently. |
| `ordered` factor | Kept ordered, with its polynomial contrasts. `mo()` requires an ordered factor with at least three levels. |
| `logical` | Handled as a two-level factor by [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html), giving one `TRUE` column. Also accepted directly as a `mo()` or `mi()` interaction multiplier. |
| `character`, as a grouping variable | Accepted. An integer, a factor and a character grouping variable give the same fit. |
| `character`, in a fixed-effect term | Works, but by delegation: [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) coerces it to a factor with the default contrasts. Level order then follows the locale collation, not your data. Convert to a factor yourself if the level order matters. |
| Matrix-valued terms | [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()`, [`scale()`](https://rdrr.io/r/base/scale.html) and matrix covariates for functional regression are supported. Their basis is frozen at fit time through `predvars`, so new data is projected onto the fitted basis rather than refitted. |
| `Date`, `POSIXct`, `difftime` | Accepted as a number. [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) reduces the column to its underlying value, so the coefficient is per day (`Date`), per second (`POSIXct`), or per unit of the `units` attribute (`difftime`). Frame assembly reports the coercion, naming the column, its unit and its origin. See “Dates and times” below. |

### Rejected

These raise an error rather than fitting something surprising.

| Input | Behavior |
|----|----|
| List column in any term | Error from [`model.frame()`](https://rdrr.io/r/stats/model.frame.html): `invalid type (list) for variable`. There is no frmtmb-level check, so the message names the variable but not the model. |
| `character` or `factor` multiplier in a `mo()` or `mi()` interaction | Error. [`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a character vector gives all `NA`, and the only symptom would be an `NA/NaN gradient evaluation` from the optimizer. The same type gate runs again on `newdata`. |
| Unordered factor in `mo()`, or fewer than three categories | Error. |
| Factor response with a family that is neither ordinal nor binomial | Error. |
| `Inf` or `NaN` in the response | Error. This is deliberately separate from `NA`, which `na.action` owns. |
| Matrix response under `mi()` | Error. |
| `cs()` grouping that is not a single factor | Error. |
| `gr(cov =)` or `gr(prec =)` that is not square, or whose dimnames miss a grouping level | Error. |
| `car(type =)` outside `"escar"`, `"icar"`, `"esicar"`, `"bym2"` | Error listing the four types. |

One case warns instead of failing: `ar1()` and `hetar1()` correlate
levels by position, so whole-number factor levels with a gap in them
(`1, 2, 4`) get a warning that the gap counts as one step, with a
pointer to `ou(num_factor(...))` for irregular spacing.

``` r

# a list column stops at model.frame(), before any frmtmb check
d$lst <- I(lapply(seq_len(30), function(i) 1:2))
frm(bf(y ~ lst) + gaussian(), data = d)
#> Error: invalid type (list) for variable 'lst'
```

### Dates and times

A `Date`, `POSIXct` or `difftime` predictor is fitted as the number
underneath it. The model is the one you asked for, but it is expressed
in an origin you did not choose:
[`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a `Date` counts
days from 1970-01-01, and on a `POSIXct` it counts seconds. So the slope
is per day or per second, and the intercept is the fitted value on
1970-01-01.

That intercept is the problem. A modern date is about 18000 as a number,
and a modern timestamp about 1.7e9, which puts the intercept far outside
the data and makes the objective badly conditioned. The same model on
days-since-the-first-day converges where the raw column reports false
convergence, with an identical slope:

``` r

d <- data.frame(day = as.Date("2020-01-01") + 0:59)
d$y <- rnorm(60, 1 + 0.02 * as.numeric(d$day - min(d$day)), 0.5)

frm(bf(y ~ day) + gaussian(), data = d)
#> Date/time column used as a number: day (Date, days since 1970-01-01).
#> The coefficient is per unit of that origin and the intercept is the
#> value at it, which is far outside the data and can stop the optimizer
#> converging. Center the column, for example as.numeric(x - min(x)), to
#> put the intercept back in range.
#> Warning: Optimizer did not report convergence: false convergence (8)
#> (Intercept) = -366.16, day = 0.0201

# the same fit, with the intercept back inside the data
d$days <- as.numeric(d$day - min(d$day))
frm(bf(y ~ days) + gaussian(), data = d)
#> (Intercept) = 1.051, days = 0.0201
```

Frame assembly emits that message once per fit. It is a message and not
a warning because the coercion is deliberate and correct in meaning;
[`suppressMessages()`](https://rdrr.io/r/base/message.html) silences it
once you have centered the column or if you want the epoch origin.

## Character options and case

All character options are matched with
[`match.arg()`](https://rdrr.io/r/base/match.arg.html), so matching is
case sensitive and accepts unambiguous prefixes. `"response"` and
`"resp"` both work for `predict(type =)`; `"Response"` does not.

| Function | Argument | Values |
|----|----|----|
| [`predict()`](https://rdrr.io/r/stats/predict.html) | `type` | `"link"`, `"response"`, `"conditional"`, `"zprob"`, `"zlink"`, `"disp"` |
| [`residuals()`](https://rdrr.io/r/stats/residuals.html) | `type` | `"response"`, `"pearson"`, `"deviance"`, `"osa"` |
| [`confint()`](https://rdrr.io/r/stats/confint.html) | `method` | `"wald"`, `"Wald"`, `"profile"`, `"uniroot"`, `"boot"` |
| [`drop1()`](https://rdrr.io/r/stats/add1.html) | `test` | `"none"`, `"Chisq"` |
| [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md) | `method` | `"wald"`, `"profile"`, `"boot"` |
| [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md) | `method` | `"epred"`, `"predict"` |
| [`anova()`](https://rdrr.io/r/stats/anova.html) on [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md) | `method`, `use` | `"D3"`, `"D1"`, `"D2"`; `"likelihood"`, `"wald"` |
| [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md) | `check_nlev_1`, `check_olre` | `"warning"`, `"ignore"`, `"stop"` |
| [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md) | `class` | `"b"`, `"Intercept"`, `"sd"`, `"cor"`, `"theta"` |
| [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md) | `band` | `"wald"`, `"profile"`, `"boot"` |
| [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) | `init` | `"stationary"`, `"estimated"`, `"uniform"` |
| [`loo_compare()`](https://aforren1.github.io/frmtmb/reference/loo.md) | `criterion` | `"loo"`, `"waic"` |
| [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md) | `on_error` | `"penalize"`, `"error"` |
| [`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md) | `type` | `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"` |

`confint(method =)` is the one place that accepts two spellings of the
same choice: `"Wald"` is listed alongside `"wald"` and normalized
immediately, because both spellings are common in other packages.

Options that do not use
[`match.arg()`](https://rdrr.io/r/base/match.arg.html) are matched
exactly, with no prefix matching: `car(type =)`, a family given as a
string, `getME(name =)`, `frm(dry_run =)` and
`frmtmb_control(optimizer =)`.

The single case-insensitive option in the package is the code vocabulary
of `cens()`. There, `"right"`, `"Right"` and `"R"` are the same code,
because a capitalized word reads as a spelling variant rather than as
garbage. The numeric codes `0`, `-1`, `1`, `2` mean the same thing.

## Naming collisions

Four vocabularies meet in a model: the columns of `data`, the
distributional parameters of the family, the nonlinear parameters a
formula declares, and the display names that
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
and
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
use. They namespace apart by construction, so a column may carry any
name at all. The rules below say what happens where two of them do meet.

**Columns and coefficients never collide.** A coefficient is stored
under its dpar prefix, so a covariate named `sigma` in the mean is
`sigma` in `fixef(fit)$mu` while the residual standard deviation is
`(Intercept)` in `fixef(fit)$sigma`. Both fit, and neither reads the
other.

**A nonlinear parameter may not also be a column.**
`bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)` on data that has a
column `a` is refused: the body would use the parameter and ignore the
column silently. Rename one of them. A declared nonlinear parameter that
no body uses is refused too, and so is a dpar formula for a parameter
the family does not have (`bf(y ~ x, shape ~ 1)` on
[`gaussian()`](https://rdrr.io/r/stats/family.html), which names `sigma`
as the one available).

**In a nonlinear body, a column wins over a distributional parameter.**
A name in an
[`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) body
resolves in this order: a nonlinear parameter of the model, then a
column of `data`, then another distributional parameter of the same
response (read per row, after its link inverse). The dpar reference is
an extension of brms, which has only the first two; putting the column
ahead of it is what keeps every body brms accepts meaning here what it
means there. So `nlf(sigma ~ ls + th * log(abs(mu)))` is a variance
function of the model’s own mean, unless `data` carries a column `mu`,
in which case it is a variance function of that column.

**The internal `.eta_<dpar>` names are unreachable from `data`.** The
objective keeps each linear predictor beside its parameter under a
reserved `.eta_` name, for the families that need the untransformed
scale. Those live in an internal list rather than in the model frame, so
a column named `.eta_mu` is an ordinary covariate and gets an ordinary
coefficient.

**In
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
a coefficient shadows a natural-scale name.**
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
evaluates expressions in one flat namespace that holds the coefficients
and the natural-scale summaries together, and the coefficient wins:

``` r

# a covariate literally named sigma
fit <- frm(bf(y ~ sigma + (1 | g)), family = gaussian(), data = dd)
hypothesis(fit, "sigma = 0")
#> hypothesis() reads 'sigma' as the coefficient of the model term of
#> that name, not the residual standard deviation; that quantity is
#> available as '.sigma'.
```

The message names both meanings and fires once per call. The shadowed
quantity keeps a name of its own: prefix it with a dot. The same holds
for a coefficient that spells out `sd_<group>__<term>`, `cor_...`, or an
autocorrelation name such as `ar1`; `.sd_g__Intercept` and `.ar1` then
reach the natural-scale value. The dot spelling exists only where a
collision does, and
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
lists both names in that case.

## What reaches the output

### Preserved

- **Row names.** The row names of `data` reach
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) and the
  `Z` matrix from `getME(fit, "Z")`. They are carried by the design
  matrix, not reassigned, so a row name survives a subset and a reorder.
- **Coefficient names.** [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md)
  share one naming scheme.
- **Factor levels and contrasts.** The level set and the contrast matrix
  of every factor are stored per linear predictor and reapplied to
  `newdata`. A level missing from `newdata` is still a known level; a
  level that is new is an error unless `allow_new_levels = TRUE`.
- **Data-dependent bases.**
  [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()` and
  [`scale()`](https://rdrr.io/r/base/scale.html) on the predictor side
  are frozen through `predvars`, so `newdata` uses the fitted centering,
  scaling and knots.
- **Row positions under `na.exclude`.** Output is padded back to the
  input length with `NA` in the dropped positions.

``` r

fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
head(names(fitted(fit)), 3)
#> [1] "obs1" "obs2" "obs3"
head(names(predict(fit)), 3)
#> [1] "obs1" "obs2" "obs3"
```

### Not preserved

- **Attributes on the response.** The response reaches the objective as
  `as.numeric(as.vector(y))`, or `storage.mode(y) <- "double"` for a
  matrix response. A [`scale()`](https://rdrr.io/r/base/scale.html)d
  response therefore arrives without its `"scaled:center"` and
  `"scaled:scale"` attributes, and a one-column matrix response arrives
  without its `dim`. If you need to invert the transformation later,
  keep the attributes yourself.
- **Column attributes such as `units` or `haven` labels.** Nothing in
  the package inspects, preserves or warns about them.
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) and
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) drop
  them, and the coefficient is reported in whatever numeric units the
  column held.
- **`dimnames` on a matrix covariate.** Design column names come from
  the term labels, not from the covariate’s own dimnames.
- **Rows dropped by `na.omit`.** They are absent from
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html). Use
  `na.action = stats::na.exclude` to get them back as `NA`.

## Runtime and how it scales

A fit has four costs, and only one of them grows quickly.

| Stage | Behavior with data size |
|----|----|
| parse | Independent of the data. Effectively free. |
| frame | Linear in `n` and in the number of random-effect columns. |
| tape | One pass over the objective. Independent of the number of optimizer iterations. |
| optimize | (number of outer iterations) x (cost of one objective evaluation). This dominates. |
| sdreport | One Hessian and one sparse solve at the optimum. Skipped unless `se = TRUE` or a method needs it. |

One objective evaluation costs an inner Newton solve, whose work is the
sparse Cholesky factorization of the inner Hessian. That matrix is
(number of random-effect coefficients) square, and its sparsity comes
from the grouping structure, so the cost tracks the number and the
overlap of the random-effect levels rather than `n` alone. The number of
outer iterations tracks the number of outer parameters and the
conditioning of the problem, not the number of rows.

Measured on
[`lme4::InstEval`](https://rdrr.io/pkg/lme4/man/InstEval.html) (n =
73421, 2972 students, 1128 lecturers) for
`y ~ service + (1 | s) + (1 | d)`, which tapes to 5 outer and 4100 inner
parameters with a 4100 x 4100 inner Hessian holding 77521 nonzeros:

| Stage    | Seconds |
|----------|---------|
| parse    | 0.00    |
| frame    | 0.14    |
| tape     | 0.66    |
| optimize | 22.78   |
| restart  | 0.44    |
| sdreport | 3.84    |
| total    | 28.63   |

Everything outside `optimize` is roughly fixed: parse, frame and tape
together are about 0.8 s, and `sdreport()` about 3.8 s. `optimize` is 80
percent of the fit, and 100 percent of `optimize` is objective
evaluation: `nlminb` itself costs nothing measurable on a five parameter
problem. The first `obj$fn` call is about 1.6 s because the inner solve
starts cold; every later one is about 0.17 s because it warm starts from
the previous point, and `obj$gr` is about 0.18 s. This model takes 40
outer iterations, so halving the iteration count roughly halves the fit.

Fit times across sizes, all wall clock seconds on the same machine:

| Model                                           | frmtmb | glmmTMB | lme4 |
|-------------------------------------------------|--------|---------|------|
| LMM, n = 180 (`sleepstudy`)                     | 0.19   | 0.08    | 0.03 |
| LMM, n = 50000, 200 groups, `(x | g)`           | 4.61   | 2.86    | 0.35 |
| Poisson GLMM, n = 100000, 500 groups            | 6.00   | 5.18    | 9.89 |
| Negative binomial with `dispformula`, n = 20000 | 1.31   | 0.71    |      |
| Zero-inflated Poisson, n = 20000                | 0.47   | 0.58    |      |
| LMM, n = 73421 (`InstEval`)                     | 21.50  |         | 4.42 |

Read this as context, not as a race.
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) is about 4.9
times faster on `InstEval` because it profiles both the fixed effects
and the residual scale out of the deviance, leaving a two parameter
problem over a purpose-built sparse Cholesky, while frmtmb hands
`nlminb` five parameters over a generic Laplace tape. The two agree on
the objective to four decimal places (118865.3077). frmtmb is faster
than `lme4` on the large Poisson GLMM, where `glmer`’s adaptive
quadrature and penalized iteratively reweighted least squares scale
worse than Laplace with exact derivatives.

Three controls change the runtime measurably on `InstEval`:

- `frmtmb_control(profile = TRUE)` takes it from 21.50 s to 13.54 s, and
  should help more as the number of fixed effects grows.
- `frmtmb_control(optimizer = "optim")` takes it to 10.44 s, but it is
  slower than `nlminb` on the intercept-only version of the same model
  at matched tolerance, so it is a per-model option rather than a better
  default.
- `REML = TRUE` takes it to 12.06 s.

`frmtmb_control(sparse_x = TRUE)` changes nothing here, because `X` is
73421 by 2. It is for genuinely sparse and wide fixed-effect designs.

## Reproducing the timings

Every number above comes from a script in `dev/` of the source
repository at <https://github.com/aforren1/frmtmb>. `dev/` is not part
of the installed package, so clone the repository to run them.

| Script | Produces |
|----|----|
| `dev/bench-env.R` | machine, R and package versions |
| `dev/bench-setup.R` | one-time install into a scratch library |
| `dev/bench-op-baseline.R` | the per-stage table |
| `dev/bench-op-profile.R` | the split of `optimize` into `obj$fn` and `obj$gr` |
| `dev/bench-op-ceiling.R` | standalone per-call costs |
| `dev/bench-op-context.R` | the frmtmb, glmmTMB and lme4 comparison |
| `dev/bench-op-tight.R` | the same comparison at matched convergence tolerance |
| `dev/bench-rtmbp-core.R` | the parallel-backend appendix |

The comparison table is one command:

``` sh
Rscript dev/bench-op-context.R
```

Run `Rscript dev/bench-setup.R` once first, and edit the scratch library
path at the top of each script to a directory you can write to. Run the
timing scripts one at a time; the numbers move by about 15 percent with
another R process on the machine.

The headline frmtmb against lme4 comparison does not need the harness.
This one line reproduces it from an installed frmtmb:

``` sh
Rscript -e 'library(frmtmb); data(InstEval, package = "lme4"); print(system.time(frm(bf(y ~ service + (1 | s) + (1 | d)) + gaussian(), data = InstEval))); print(system.time(lme4::lmer(y ~ service + (1 | s) + (1 | d), data = InstEval, REML = FALSE)))'
```

The findings behind these numbers, including why parallelizing the outer
optimizer does not help and why an RTMB tape cannot be shipped to a
worker process, are written up in `dev/benchmarks.md`.
