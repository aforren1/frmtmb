# The two fitting options every mixture-type family refuses, in its own name.

A mixture likelihood is invariant to permuting its components, so the
location coefficients enter a multimodal objective. Both `REML` and
`profile = TRUE` integrate those coefficients out with a Laplace
approximation about a single inner mode, which is not defined here: the
inner Newton solve walks between the component modes and the fit either
dies at "NA/NaN gradient evaluation" or reports an optimum with a
gradient near 1e9. Quadrature is unaffected, because it marginalizes the
random effects, not the coefficients.

The conservative all-`FALSE` default is right for a family whose
likelihood does not factorize: it starts fully refused and opts in. It
is exactly wrong for a family whose likelihood DOES factorize and which
carries a structure only to hold the two or three refusals that belong
to it, because there every other capability already works and a flag
would only take it away. This is the starting point for that second
kind, so such a family names what it refuses and nothing else.

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
`posterior_predict()`,
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
among them.

## Usage

``` r
response_mean(fam, dpars, aterms)

mixture_multimodal_refusals(what)

mixture_posterior(fit)

as_frmtmb_family(x)

eval_dpars(fit, b = fit$estimates[["b"]])

single_response(fit, what)

fit_extras(fit)

dpar_linpred(frame, params, resp, dpar)

structure_supports_all(...)

frame_block_of(frame, resp)
```

## Arguments

- fam:

  A `frmtmb_family`.

- dpars:

  Named list of numeric dpar values, as one element of `eval_dpars()`.

- aterms:

  Evaluated addition terms for the response,
  `fit$frame$aterm_values[[resp]]`.

- what:

  For `single_response()`, the calling function as it should open a
  refusal: `"latent_probs()"`. For `mixture_multimodal_refusals()`, the
  family as the refusal should name it: `"an lca() family"`.

- fit:

  A `frmtmb_fit`.

- x:

  A family: a `frmtmb_family`, a
  [`stats::family()`](https://rdrr.io/r/stats/family.html) object, a
  family constructor, or a family name.

- b:

  Random-effect vector to evaluate at, defaulting to the fit's own
  estimates. `NULL` drops the random-effect contribution.

- frame:

  A `fit$frame`, or the frame a `check_fit` slot is given.

- params:

  A parameter list in the frame's layout: a fit's `estimates`, or the
  starting template a `check_fit` slot is given.

- resp, dpar:

  The response name and the distributional parameter name, as
  `eval_dpars()` returns them.

- ...:

  Capability flags to override, named as in the `supports` table of
  [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md).

## Value

`single_response()` a response spec; `eval_dpars()` a nested list of
numeric vectors; `fit_extras()` a named list or `NULL`; `dpar_linpred()`
a numeric vector or `NULL`; `response_mean()` a numeric vector;
`as_frmtmb_family()` a `frmtmb_family`; `frame_block_of()` a list or
`NULL`; `structure_supports_all()` a named list of logicals;
`mixture_posterior()` a matrix of class probabilities;
`mixture_multimodal_refusals()` a list of two refusal strings.

## Details

`what` is the family as the sentence should name it. These used to be
one message listing every mixture-type family, raised from one gate in
fit.R that knew all three by name; each family states its own now.

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

- `frame_block_of()`:

  The block a `frame_block` slot built for one response, or `NULL` when
  it built none. This is the read half of that slot: the block is
  written once at frame assembly and read back by every slot that runs
  at the estimates, and where a frame keeps its blocks is not a layout
  an extension should spell for itself.

- `structure_supports_all()`:

  A `supports` list with every capability TRUE and the named ones
  overridden. The starting point for the second kind of structure: a
  family whose likelihood DOES factorize over rows and which declares a
  structure only to carry the two or three refusals that belong to it.
  Writing the TRUEs out by hand instead would silently take away each
  new capability frmtmb adds.

- `mixture_posterior()`:

  The posterior class probabilities of any family that implements the
  `fam$mix` component interface, which is what
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md),
  [`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
  and a latent-class family all do. A family implementing that interface
  uses this as its whole `latent_probs` slot rather than deriving the
  same quantity a second time.

- `mixture_multimodal_refusals()`:

  The `refusals` entries for `reml` and `profile` owed by a family whose
  likelihood is invariant to permuting its components. `what` is the
  family as the sentence should name it, for example
  `"an lca() family"`. The fact refused is the same one in every such
  family, so the sentence is issued from one place rather than copied.

## See also

[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for the protocol these serve,
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
for the family object they read, and the registries an extension fills
from its own `.onLoad()`:
[`frmtmb_register_frame_check()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_frame_check.md),
[`frmtmb_register_aterm()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_aterm.md)
and
[`frmtmb_register_compat()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_compat.md)

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
