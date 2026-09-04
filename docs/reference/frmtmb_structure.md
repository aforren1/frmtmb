# Declare a non-rowwise likelihood to the core

A `frmtmb_structure()` is what a family carries when its likelihood does
not factorize over the rows of the data: a group-level
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md), a
hidden Markov chain, a latent class measurement model. It is one object
with one contract, so the core needs no branch per family and a
structured family can live in another package. Attach it with
`frmtmb_family(structure = )`.

## Usage

``` r
frmtmb_structure(
  frame_vars = NULL,
  keep_na = FALSE,
  check_spec = NULL,
  frame_block = NULL,
  check_frame = NULL,
  check_fit = NULL,
  loglik = NULL,
  unit = NULL,
  fitted_mean = NULL,
  fitted_var = NULL,
  latent_probs = NULL,
  sim_ctx = NULL,
  supports = list(),
  refusals = list()
)
```

## Arguments

- frame_vars:

  `function(fam)` returning a list of language objects whose variables
  must be in the model frame but belong to no linear predictor: a
  grouping column, a time column, a sequence id. It runs before any data
  is seen, so it may read only the family.

- keep_na:

  `TRUE` means an `NA` in the response is data the family reads, so the
  row survives `na.action`. `NA`s in every other variable still drop the
  row. A family that keeps them must handle them in `frame_block` (a
  mask, a placeholder) or in `loglik`.

- check_spec:

  `function(resp, spec, av)` run before the generic addition-term
  guards, so a structured family refuses a term in its own words rather
  than through a missing CDF further down. It sees the response spec,
  the whole spec (for the univariate check) and the evaluated addition
  terms. It returns nothing and stops to refuse.

- frame_block:

  `function(resp, spec, av, mf, y, n)` run once at frame assembly, after
  `y` is coerced and before the random-effect blocks are built. It
  returns the block; see the Block section.

- check_frame:

  `function(spec, frame)` run after the predictors and random-effect
  blocks exist, for a refusal that depends on the design rather than on
  the response. It sees the assembled frame without the parameter
  template.

- check_fit:

  `function(resp, frame, template)` run once the starting values are
  resolved and before the optimizer, for a refusal or a warning that
  depends on WHERE the fit starts rather than on the model. `template`
  is the starting parameter list. The capability flags cover the fitting
  options a structure cannot answer at all; this covers a start the
  structure can answer only badly, such as a start sitting on a
  label-symmetry axis.

- loglik:

  `function(y, dpars, aterms, weights, block, extra)` returning the
  taped log-likelihood of the WHOLE response as one AD scalar:
  `sum(log-likelihood)`, not a negative. `y` is the response after any
  `block[["y"]]` replacement; `dpars` are the evaluated distributional
  parameters on the natural scale, one AD vector of length `n` or `1`
  each; `weights` are the effective row weights with cluster weights
  folded in, or `1`; `extra` are the family's extra parameters as an AD
  list, in the order `extra_pars` declared them. The family decides what
  a row weight means for a likelihood that is not rowwise, and may have
  refused weights in `check_spec`. It must not call
  [`RTMB::OBS()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html).

- unit:

  One noun phrase naming the smallest independent unit of `loglik`, as
  it should read in the middle of a sentence: "a hidden-Markov
  sequence", "a group-level mixture". The core quotes it where it must
  explain that a per-OBSERVATION quantity does not exist, which today is
  [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) and
  [`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md)
  refusing a pointwise log-likelihood matrix. Ignored when `loglik` is
  `NULL`.

- fitted_mean, fitted_var:

  `function(fit, block)` giving the conditional mean and variance of
  each row GIVEN the whole observed response, for
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  `predict(type = "response")` on the training data, and pearson
  residuals. `NULL` means "use the rowwise family mean", which is what a
  group-level mixture wants. A family with no mean supplies a function
  that stops. Both run at the estimates, outside the tape.

- latent_probs:

  `function(fit, block)` returning one matrix of posterior latent-state
  probabilities with column names, `n` rows or one row per group.

- sim_ctx:

  `function(ctx)` drawing the whole response, the structured-simulator
  contract of
  [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  with `ctx[["block"]]` carrying the block. One implementation serves
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  `posterior_predict()` and
  [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md).

- supports:

  Named logical vector or list of capability flags; see Capability
  flags. Unnamed entries and unknown names are refused.

- refusals:

  Named list of one message per `FALSE` flag.

## Value

An object of class `frmtmb_structure`.

## Details

Every slot defaults to the rowwise behavior, which is what a family that
changes the likelihood and nothing else needs: give `loglik` and leave
the rest alone.

## Extending frmtmb from another package

This is the seam for a family whose likelihood does not factorize over
rows. Attach a structure to a family and the core needs no branch that
names the family, which is what lets the family ship somewhere else. The
other seams are registries, each filled from a contributing package's
`.onLoad()`:
[`frmtmb_register_frame_check()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_frame_check.md)
for a feature that is NOT a family and still has to refuse a data
problem only the assembled frame shows;
[`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.md)
for an addition term of the family's own, so its per-row data keeps the
spelling its literature uses instead of `vint()`; and
[`frmtmb_register_compat()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.md)
for the rows
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
reports about it.
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
is the read-only accessor set every kind of extension uses to reach a
fit. The worked example of the frame-check seam is the `frmtmb.ode`
package, whose `frm_ode()` refuses a dynamics input that varies inside a
solve group; its registration is three lines of `.onLoad()`.

`loglik = NULL` keeps the family's own rowwise `lpdf` and makes the
structure a CAPABILITY DECLARATION instead. That is what a family whose
likelihood does factorize per row but which still has to refuse things
in its own words wants: `lca()` is rowwise per subject and still refuses
`REML`, `residuals(type = "osa")` and random effects in its gating
predictor, and each of those refusals is a property of the family rather
than of the method that meets it.

## The block

`frame_block()` returns the **block**: a plain list of DATA, stored at
`fit$frame$blocks[[response]]` and handed back to `loglik()`,
`fitted_mean()`, `fitted_var()`,
[`latent_probs()`](https://aforren1.github.io/frmtmb/reference/latent_probs.md)
and `sim_ctx()`. It is saved inside the fit and rebuilt by
[`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md), so it
must hold no AD values and no closure that captures the model frame.
Read and write it with `[[ ]]` only: `$` partial matching is how a `mix`
read once returned `mix_g`.

Three names in it are reserved, because the core reads them:

- `y`:

  If present, replaces the response vector for every later stage, which
  is how a family fills a placeholder in for an `NA` it keeps. Otherwise
  the response is unchanged.

- `miss`:

  A logical `n`-vector. Residuals are `NA` at these rows. Optional.

- `mask`:

  A numeric 0/1 `n`-vector the family multiplies into its own density.
  The core does not read it; the name is reserved so that every
  structured family spells it the same way.

Everything else in the block belongs to the family.

\[ \]: R:%20

## Capability flags

`supports` is a named logical vector or list. With a `loglik`, every
name defaults to `FALSE`: a non-rowwise likelihood starts fully refused
and opts in. Without one (a capability declaration for a family that
keeps its rowwise likelihood), every name defaults to `TRUE` and the
family names only what it refuses:

- `reml`:

  `REML = TRUE`.

- `quadrature`:

  `quadrature =` other than the Laplace default.

- `profile`:

  `frmtmb_control(profile = TRUE)`.

- `newdata_response`:

  `predict(newdata =, type = "response")`.

- `se_fit_response`:

  `predict(se.fit = TRUE, type = "response")`.

- `re_form`:

  `re.form =` in [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- `conditional_effects`:

  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md).

- `osa`:

  `residuals(type = "osa")`.

- `deviance`:

  `residuals(type = "deviance")`.

- `multivariate`:

  [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) and
  `rescor = TRUE`.

- `cens_trunc`:

  `cens()` and [`trunc()`](https://rdrr.io/r/base/Round.html).

- `mi`:

  `mi()` on the same response.

- `cluster_robust`:

  [`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md).

`FALSE` refuses the capability outright. `TRUE` means only that the
structure does not stand in the way: the core's ordinary rules still
apply, so a mixture that sets `deviance = TRUE` is still refused a
deviance residual by the family's missing `dev_fn`, with the message
that refusal has always carried.

Prediction on the link scale with `dpar =` is always available and is
not a flag: the linear predictors are rowwise and belong to the core.

`reml`, `quadrature` and `profile` are read at fit time; the rest by the
method each one names. `multivariate`, `cens_trunc` and `mi` describe
model shapes a family will usually want to refuse from `check_spec`
instead, where it can say which addition term was the problem.

`refusals[[flag]]` is the WHOLE message a user sees when that capability
is asked for, so a family explains its own refusal in its own words
instead of through a generic sentence. A flag with no string of its own
gets a generic one naming the family. A name may carry a context suffix,
`"re_form.simulate"`, when one flag is refused for two different reasons
in two places; the core falls back to the bare flag name.

## See also

[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the accessors a slot may use to read a fit,
[`latent_probs()`](https://aforren1.github.io/frmtmb/reference/latent_probs.md)
for the generic `latent_probs` answers,
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md),
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)

## Examples

``` r
# a structure whose likelihood is the rowwise one, written out: the
# smallest thing the protocol accepts
st <- frmtmb_structure(
  loglik = function(y, dpars, aterms, weights, block, extra) {
    sum(weights * RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE))
  },
  supports = list(conditional_effects = TRUE)
)
st
#> <frmtmb_structure>
#>   slots:    loglik
#>   keeps NA: FALSE
#>   supports: conditional_effects
```
