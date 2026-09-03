# Accessors for a structured family written outside frmtmb

The read-only surface a
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
slot may use. A structure's `fitted_mean`, `fitted_var`, `latent_probs`
and `sim_ctx` run at the estimates and outside the AD tape, and they are
handed a `frmtmb_fit`; these are how such a slot reads what it needs out
of one without reaching into the fit's internal layout. The rest of a
family's needs are already public:
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md),
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md),
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md),
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
among them.

## Usage

``` r
response_mean(fam, dpars, aterms)

as_frmtmb_family(x)

eval_dpars(fit, b = fit$estimates[["b"]])

single_response(fit, what)

fit_extras(fit)

dpar_linpred(frame, params, resp, dpar)
```

## Arguments

- fam:

  A `frmtmb_family`.

- dpars:

  Named list of numeric dpar values, as one element of `eval_dpars()`.

- aterms:

  Evaluated addition terms for the response,
  `fit$frame$aterm_values[[resp]]`.

- x:

  A family: a `frmtmb_family`, a
  [`stats::family()`](https://rdrr.io/r/stats/family.html) object, a
  family constructor, or a family name.

- fit:

  A `frmtmb_fit`.

- b:

  Random-effect vector to evaluate at, defaulting to the fit's own
  estimates. `NULL` drops the random-effect contribution.

- what:

  The calling function, as it should open a refusal: `"latent_probs()"`.

- frame:

  A `fit$frame`, or the frame a `check_fit` slot is given.

- params:

  A parameter list in the frame's layout: a fit's `estimates`, or the
  starting template a `check_fit` slot is given.

- resp, dpar:

  The response name and the distributional parameter name, as
  `eval_dpars()` returns them.

## Value

`single_response()` a response spec; `eval_dpars()` a nested list of
numeric vectors; `fit_extras()` a named list or `NULL`; `dpar_linpred()`
a numeric vector or `NULL`; `response_mean()` a numeric vector;
`as_frmtmb_family()` a `frmtmb_family`.

## Details

They are documented together because they are one contract with one
promise: what they return is stable, so a family in another package can
be written against them. Everything else in a fit, a frame or a family
object is internal and may be renamed without notice.

- `single_response()`:

  The one response spec of a univariate fit, or the refusal a method
  owes a multivariate one. `what` names the caller and opens the
  message. Almost every structured family is univariate, so this is the
  first line of most slots. It was named `uni_resp()` while it was
  internal.

- `eval_dpars()`:

  The distributional parameters of every response at the estimates, on
  their natural (response) scale, evaluated over the training rows:
  `out[[response]][[dpar]]` is a numeric vector of length `n` or `1`.
  `b` supplies a different random-effect vector, which is how a family
  evaluates at the modes of a resampled fit. This is the same evaluation
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) run.

- `fit_extras()`:

  The family's extra (non-dpar) parameters at the estimates, as a named
  list in the order `extra_pars` declared them, or `NULL` when the
  family declared none. Item profiles, ordinal thresholds and class
  covariances arrive here.

- `dpar_linpred()`:

  The FIXED-effect part of one dpar's linear predictor, on the link
  scale, at a parameter list in the frame's own layout, or `NULL` when
  that response and dpar have no linear predictor with columns. Random
  effects are excluded on purpose: this exists for a `check_fit` slot,
  which runs before the inner modes are solved and cannot see them.

- `response_mean()`:

  The expected response of a family at given dpar values, with
  [`trunc()`](https://rdrr.io/r/base/Round.html) bounds applied when the
  response carries them. A structured family composing per-state or
  per-class means calls this once per state rather than reaching into
  the wrapped family's `post` slot.

- `as_frmtmb_family()`:

  The internal representation of a family given as a `frmtmb_family`, a
  [`stats::family()`](https://rdrr.io/r/stats/family.html), a bare
  family constructor or a name. A family that WRAPS another one calls
  this on its argument, so that `hmm(2, gaussian)`, `hmm(2, gaussian())`
  and `hmm(2, "gaussian")` all mean the same thing.

## See also

[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for the protocol these serve, and
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
for the family object they read

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(50))
dd$y <- rnorm(50, 1 + 2 * dd$x, 0.5)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)

rspec <- single_response(fit, "my_family_probs()")
dp <- eval_dpars(fit)[[rspec$resp_name]]
str(dp)
#> List of 2
#>  $ mu   : Named num [1:50] -0.178 1.424 -0.591 4.215 1.712 ...
#>   ..- attr(*, "names")= chr [1:50] "1" "2" "3" "4" ...
#>  $ sigma: Named num [1:50] 0.479 0.479 0.479 0.479 0.479 ...
#>   ..- attr(*, "names")= chr [1:50] "1" "2" "3" "4" ...
head(response_mean(rspec$family, dp, list()))
#>          1          2          3          4          5          6 
#> -0.1776896  1.4240554 -0.5912754  4.2151815  1.7124624 -0.5613002 
fit_extras(fit)                       # NULL: no extra parameters
#> NULL
head(dpar_linpred(fit$frame, fit$estimates, rspec$resp_name, "mu"))
#>          1          2          3          4          5          6 
#> -0.1776896  1.4240554 -0.5912754  4.2151815  1.7124624 -0.5613002 
as_frmtmb_family(gaussian)
#> <frmtmb family> gaussian
#>   dpars: mu (identity), sigma (log)
```
