# The structured-family protocol

Status: design, not implemented. Written 2026-09-03 against v0.43.0.

## Problem

Three families have a likelihood that does not factorize over rows:
`mixture(groups =)`, `hmm()` and `lca()`. Each one reaches the core
through its own named slot and its own branch:

| concern | mixture | hmm | lca |
|---|---|---|---|
| family marker | `fam$mix`, `fam$mix_groups` | `fam[["hmm"]]` | `fam$lca` + borrowed `fam$mix` |
| frame slot | `frame$mix_g[[r]]` | `frame$hmm_g[[r]]` | `frame$mix_g[[r]]` |
| frame variables | frame.R:794, 986 | frame.R:797, 991 | none |
| NA rows kept | no | frame.R:1029 | `fam$na_response` |
| pre-aterm checks | frame.R:1288 (inline) | `hmm_check_aterms()` | `valid_y()` |
| block builder | frame.R:1302 (inline) | `hmm_frame_block()` | none |
| post-predictor checks | none | none | `check_lca_structure()` |
| fit-time refusals | `has_mixture()` gate, fit.R:551 | `hmm_check_fit()` | via `has_mixture()` |
| objective branch | objective.R:325 | objective.R:313 | via mixture |
| fitted / predict | rowwise | `hmm_mean_response()` | `post$mean_fn` refuses |
| residuals | rowwise | predict.R:2213 | predict.R:2248 |
| simulate | `sim_ctx` reads `ctx$mix_g` | `sim_ctx` | rowwise `sim` |
| conditional effects | allowed | `ce_hmm_check()` | allowed |
| latent probabilities | `mixture_probs()` | `hmm_probs()` | `lca_probs()` |

Core files that name one of the three: fit.R, frame.R, objective.R,
parse.R, predict.R, conditional-effects.R, families.R. A fourth
structured family would add a fourth column and touch all of them.

## Goal

One named slot, `fam$structure`, built by one exported constructor,
consumed by core at a fixed set of call sites, with no family name in
core. `mixture()` moves onto the protocol first and keeps its test
suite unchanged; that is the acceptance test for the seam. `hmm()` and
`lca()` follow and can then live in another package.

Out of scope: R-side residual structures (`autocor`) and `frm_ode()`.
Autocor has the same shape (it replaces `lpdf` for a response) but hangs
off the response spec, not the family, and should be revisited once
this protocol is in. ODE is a nonlinear-body seam, not a family one.

## The constructor

```r
frmtmb_structure(
  frame_vars    = NULL,   # function(fam) -> list of language objects
  keep_na       = FALSE,  # logical: response NAs survive na.action
  check_spec    = NULL,   # function(resp, spec, av)
  frame_block   = NULL,   # function(resp, spec, av, mf, y, n) -> block
  check_frame   = NULL,   # function(spec, frame)
  loglik        = ,       # function(y, dpars, aterms, weights, block, extra) -> AD scalar
  fitted_mean   = NULL,   # function(fit, block) -> numeric(n) or NULL
  fitted_var    = NULL,   # function(fit, block) -> numeric(n) or NULL
  latent_probs  = NULL,   # function(fit, block) -> matrix
  sim_ctx       = NULL,   # function(ctx) -> response draw
  supports      = list(), # named logicals, see below
  refusals      = list()  # named strings, one per FALSE in supports
)
```

`loglik` is the only required slot. Everything else defaults to "the
rowwise behavior", which is what a structured family that only changes
the likelihood needs.

The constructor validates types and the `supports` names, stamps class
`frmtmb_structure`, and is the object `frmtmb_family(structure =)`
stores. `frmtmb_family()` gains one argument. The family's own `lpdf`
stays required for the rowwise contract and may be a refusing stub, as
`hmm()`'s is today.

## Slot contracts

### `frame_vars(fam)`

Returns language objects whose variables must be in the model frame
but belong to no linear predictor: a grouping column, a time column, a
sequence id. Core adds them to `nonpredictor_frame_vars()` and to the
frame formula's parts (today frame.R:790 and frame.R:984). Called
before any data is seen, so it may only read the family object.

### `keep_na`

`TRUE` means an NA in the response is data the family reads, so the
row survives `na.action`. Replaces both `fam$na_response` and the hmm
special case at frame.R:1029. NAs in every other variable still drop
the row. A family that keeps NAs must handle them in `frame_block`
(mask, placeholder) or in `loglik`.

### `check_spec(resp, spec, av)`

Runs before the generic aterm guards at frame.R:1131, so a structured
family refuses in its own words rather than through a missing CDF.
Sees the response spec, the whole spec (for the univariate check), and
the evaluated aterms. Returns nothing; stops to refuse.

### `frame_block(resp, spec, av, mf, y, n)`

Runs once at frame assembly, after `y` is coerced and before random
effects are built. Returns a plain list, the **block**, stored at
`frame$blocks[[resp_name]]`. Rules:

- Data only. No AD values, no closures that capture `mf`. The block is
  saved inside the fit and rebuilt by `refit()`.
- `[[` access only, on both sides; `$` partial matching is how
  `ctx$mix` once read `ctx$mix_g`.
- Reserved names core reads:
  - `y`: if present, replaces the response vector for every later
    stage (placeholder-filled, as hmm does). Otherwise `y` is unchanged.
  - `miss`: logical `n`-vector. Core sets residuals to `NA` at these
    rows and excludes them from `napred()` padding. Optional.
  - `mask`: numeric 0/1 `n`-vector the family multiplies into its
    emission density. Core does not read it; it is reserved so every
    structured family spells it the same way.
- Everything else is the family's own.

### `check_frame(spec, frame)`

Runs after the predictors and random-effect blocks exist (today
frame.R:2243), for refusals that depend on the design rather than the
response: `lca()` refusing random effects, a future family refusing a
smooth in a state predictor. Sees the assembled frame minus the
parameter template.

### `loglik(y, dpars, aterms, weights, block, extra)`

The taped log-likelihood of the whole response, as one AD scalar.
Called from the objective's response loop in place of the rowwise
`row_lpdf()` sum. Arguments:

- `y`: the response after any `block$y` replacement.
- `dpars`: the evaluated distributional parameters, one AD vector per
  dpar, length `n` or 1, on the natural scale.
- `aterms`: the evaluated aterms for this response.
- `weights`: the effective row weights, cluster weights folded in, or
  `1`. The family decides what a row weight means for a non-rowwise
  likelihood; it may have refused weights in `check_spec`.
- `block`: the frame block.
- `extra`: the family's extra parameters as an AD list, in the order
  `extra_pars()` declared them.

Must return `sum(log-likelihood)`, not a negative. Must not touch
`RTMB::OBS()`: one-step-ahead residuals are refused generically for any
structured family (see `supports`).

### `fitted_mean(fit, block)` and `fitted_var(fit, block)`

The conditional expectation and variance of each row given the whole
observed response, for `fitted()`, `predict(type = "response")` on the
training data, and pearson residuals. `NULL` means "use the rowwise
family mean", which is what a group-level mixture wants. A family with
no mean (lca) supplies a function that stops, as `post$mean_fn` does
today. Both run at the estimates, outside the tape.

### `latent_probs(fit, block)`

One `n`-by-`K` (or `n_groups`-by-`K`) matrix of posterior latent-state
probabilities, with column names. One exported generic,
`latent_probs(fit)`, dispatches here. `mixture_probs()`, `hmm_probs()`
and `lca_probs()` become thin aliases that check the family and call
it. Decoding passes that are not probabilities (`hmm_viterbi()`) stay
family exports.

### `sim_ctx(ctx)`

Unchanged from the current structured-simulator contract, with one
rename: `ctx[["block"]]` replaces `ctx[["mix_g"]]`. The context is
built identically by `simulate()`, `posterior_predict()` and
`frm_simulate()`, so a structured family has one simulator.

### `supports` and `refusals`

Named logicals with these names and defaults:

| name | default | what `FALSE` refuses |
|---|---|---|
| `reml` | `FALSE` | `REML = TRUE` |
| `quadrature` | `FALSE` | `quadrature =` other than Laplace |
| `profile` | `FALSE` | `frmtmb_control(profile = TRUE)` |
| `newdata_response` | `FALSE` | `predict(newdata =, type = "response")` |
| `se_fit_response` | `FALSE` | `predict(se.fit = TRUE, type = "response")` |
| `re_form` | `FALSE` | `re.form =` in predict and simulate |
| `conditional_effects` | `FALSE` | `conditional_effects()` |
| `osa` | `FALSE` | `residuals(type = "osa")` |
| `deviance` | `FALSE` | `residuals(type = "deviance")` |
| `multivariate` | `FALSE` | `mvbf()` and `rescor = TRUE` |
| `cens_trunc` | `FALSE` | `cens()`, `trunc()` |
| `mi` | `FALSE` | `mi()` on the same response |

Link-scale prediction with `dpar =` is always available and is not a
flag: the linear predictors are rowwise and belong to core.

`refusals[[name]]` is the sentence core appends after its generic
lead-in ("... is not available for a `hmm(3, gaussian)` family: "). The
hand-written explanations in predict.R and conditional-effects.R move
here verbatim, so the messages users see do not change. A `FALSE` with
no sentence gets a generic one. The defaults are conservative on
purpose: a new structured family starts fully refused and opts in.

The mixture reference implementation sets `conditional_effects`,
`newdata_response`, `se_fit_response` and `re_form` to `TRUE`, because
its per-row mean is rowwise.

## Core call sites after the refactor

Exactly these, each a null check on `fam$structure` followed by one
slot call. No family name appears in any of them.

| file | replaces |
|---|---|
| frame.R `nonpredictor_frame_vars()` | mix_groups and hmm expr lines |
| frame.R frame-formula parts | same |
| frame.R NA exemption | `na_response` and hmm filter |
| frame.R before aterm guards | `hmm_check_aterms()` |
| frame.R response loop | inline mixture checks, `hmm_frame_block()`, `mix_g` and `hmm_g` |
| frame.R after predictors | `check_lca_structure()` |
| frame.R return | `blocks =` replaces `mix_g =` and `hmm_g =` |
| fit.R | `hmm_check_fit()` and the `has_mixture()` REML/quadrature gate become `supports` checks |
| objective.R response loop | both branches become one `st$loglik()` call |
| predict.R `predict()` | hmm response branch becomes `supports` checks plus `fitted_mean` |
| predict.R `fitted()` | same |
| predict.R `residuals()` | hmm and lca branches become `supports` checks, `fitted_mean`, `fitted_var`, `block$miss` |
| predict.R `simulate()` | hmm refusals become `supports` checks |
| conditional-effects.R | `ce_hmm_check()` becomes a `supports` check |
| families.R `sim_context()` | `block =` replaces `mix_g =` |

Untouched, because they are already family-generic: `default_forms`
(parse.R), `extra_pars`, `primary_dpars`, `init_dpars`, `valid_y`,
`post`, `sim`, `sim_ctx`.

## What the extension package needs exported

Functions hmm.R and lca.R call today that are internal:

- `as_frmtmb_family()`, `response_mean()`, `eval_dpars()`,
  `fit_extras()`: read-only accessors a structured family needs to
  compute means and probabilities at the estimates.
- `uni_resp()`: rename to something public, `single_response()`.
- `linpred_key()`: check whether hmm needs it or a public accessor
  covers it.
- `assemble_frame()` from lca.R: find out why. A family should not
  assemble frames. Probably `lca_profiles()` rebuilding an item table;
  it should read the fit's frame instead.
- `has_mixture()` from lca.R: goes away with the `supports` flags.

Everything else hmm and lca call (`frm`, `bf`, `mvbf`, `fixef`,
`set_prior`, `frm_allfit`, `frm_simulate`, `frm_sample`,
`conditional_effects`, `hypothesis`, `posterior_predict`,
`check_laplace`) is already exported.

## Boundary test

One test in the core suite: grep every file in `R/` except the
structured families' own for the strings `hmm`, `lca`, `mix_g`,
`hmm_g`, `mixture`, and fail on a hit outside roxygen text. Run it
from the start of the refactor so the branches cannot creep back.

## Order of work

1. `frmtmb_structure()` constructor and validator. `frmtmb_family()`
   gains `structure =`. No behavior change.
2. Frame: `blocks` slot, `frame_vars`, `keep_na`, `check_spec`,
   `frame_block`, `check_frame` call sites. Old slots still populated.
3. Objective: the single `loglik` branch, with mixture and hmm
   forwarding to it. Old branches deleted.
4. Predict, fitted, residuals, simulate, conditional effects: the
   `supports` checks and `fitted_mean` / `fitted_var` / `block$miss`.
5. `mixture()` onto the protocol. `test-lca.R`, `test-mvn-mixture.R`
   and the mixture tests in `test-families.R` must pass unchanged.
   Delete `mix_g`.
6. `hmm()` onto the protocol. `test-hmm.R` unchanged. Delete `hmm_g`
   and every hmm string in core. Boundary test goes green.
7. `lca()` onto the protocol; drop its borrowed `fam$mix`.
8. `latent_probs()` generic; aliases.
9. Export the accessor list above, document `frmtmb_structure()` as
   the extension API with hmm as the worked example.
10. Move hmm.R, lca.R and their tests to the extension package.

Steps 1 through 4 are pure refactor with the old slots alive, so each
can land green on its own.

## Open questions

- Should `loglik` receive `weights` at all? hmm refuses them, mixture
  multiplies. Passing them and documenting the choice is simpler than
  a second flag, but a family that ignores them silently is a trap.
  Proposal: pass them, and have the validator warn when `loglik`'s
  formals omit `weights`.
- `fitted_mean` on newdata. hmm cannot, mixture can. The
  `newdata_response` flag covers the refusal, but a family that opts
  in must also accept a newdata block, which `frame_block` would need
  to build without a `y`. Defer until a family needs it.
- Whether autocor joins the protocol as a "residual structure" with
  the same `loglik` signature. It would remove the third
  response-loop branch and the `acs` frame slot. Worth a separate
  note once the family protocol is in.
