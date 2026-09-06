# The structured-family protocol

Status: approved 2026-09-03; COMPLETE. Steps 1 through 9 implemented
(1-5 at v0.44.0 and v0.45.0, 6-9 after); step 10, the move to another
package, landed after v0.47.0.
Written 2026-09-03 against v0.43.0; revised same date with the
packaging decision (monorepo), the boundary-test ratchet, and the
interop split notes.

Revised 2026-09-05 with the two factorization slots (`loglik_row`,
`loglik_group`), the block's reserved `group`, and the
addition-term allow-list. Those are additions to a protocol that was
otherwise complete, not a redesign of it.

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

This goal is met. Both families left for `extensions/frmtmb.latent` in
step 10, carrying their tests and their compatibility rules with them
and reaching the core only through its exported extension API; the
boundary test is no longer a ratchet with exempt files but a
zero-tolerance assertion over the whole of `R/`.

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
  loglik_row    = NULL,   # same signature -> one value per ROW
  loglik_group  = NULL,   # same signature -> one per block$group LEVEL
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

Two more slots landed on 2026-09-05, `loglik_row` and
`loglik_group`, which say how finely the whole-response likelihood
factorizes; see their section below and the addition-term allow-list
note at the end of this document.

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
  - `group`: factor or integer `n`-vector, the family's own
    independent unit per row. What `loglik_group` returns one value
    per, and what core checks the importance grouping against.
    Optional, and required with `loglik_group`.
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

### `loglik_row(...)` and `loglik_group(...)`

ADDED 2026-09-05, after the reinforcement-learning lane found the gap
and the review corrected the fix. `loglik` returns one number for the
whole response, which is everything the objective needs and less than
everything else needs. Three consumers want the pieces:
`frm(importance =)` resamples one grouping level at a time, `loo()` and
`waic()` leave one unit out at a time, and a deviance residual compares
one row against its saturated fit.

Both slots take `loglik`'s arguments and return the same quantity
factorized instead of summed. A family declares whichever factorization
it HAS, and the finest one it has:

- `loglik_row`, one value per row: the conditional log-density of that
  row given whatever the family's factorization conditions on. It must
  sum to `loglik`, and each value must depend only on its own group's
  random effects. The core sums it into groups with the sparse
  indicator `imp_group_map()` already builds.
- `loglik_group`, one value per level of the block's `group`: for a
  family whose finest factorization IS the group, which a forward
  recursion over a sequence is.

Neither is a substitute for `unit`, and this is the correction the
review made to the first proposal. `unit` declares what may honestly be
LEFT OUT; these declare how finely the likelihood factorizes. They are
different questions and the reference consumer is the counterexample
that proves it: `rw_delta`'s finest factorization is the TRIAL, and its
`unit` is one subject's whole sequence, because dropping a trial
changes every later trial's value store. A slot keyed on `unit` would
have the protocol promising a per-row quantity out of per-group data.

Declaring either without `loglik` is refused. A family whose likelihood
is already rowwise has these quantities through its own `lpdf`, and a
slot there would be a second definition of the same numbers with
nothing keeping the two equal. `lca()` is that family and declares
neither; `hmm()` declares `loglik_group` and no `loglik_row`, because a
row's emission density is not its contribution to the likelihood.

STACKING is part of the contract, not an accident of the caller. The
importance correction evaluates the whole design once per draw, stacked
at `ridx <- rep.int(seq_len(n), nd)`, so it calls these slots with `y`,
the dpars and the addition terms each repeated `nrep` times: entry `j`
is original row `((j - 1) %% n) + 1` of replicate `((j - 1) %/% n) + 1`,
with `n = length(block[["group"]])` and `nrep = NROW(y) / n`.
`loglik_row` returns `n * nrep` values in that order and `loglik_group`
returns `ng * nrep`, replicate-major. A sequential family runs its
recursion over unit-crossed-with-replicate rather than over unit; where
its loop is already vectorized across units, that is the same loop over
a longer vector.

A family that gets the stacking wrong is caught rather than believed.
`imp_verify()` already compares the correction's per-group values with
the plain objective, per group and in total, at the first freeze, and
it is what turns "the family says it factorizes" into a checked claim.

SATURATED VALUES. A deviance residual needs `2 * (saturated - fitted)`,
and no conditional log-density supplies the first half: the saturated
fit has one parameter per observation and only the family knows what
its density reaches there. Rather than a third slot, `loglik_row` may
attach `attr(x, "saturated")`, a numeric vector of the same length,
which the core reads outside the tape. Without it, deviance stays
refused with a sentence naming what is missing. The alternative
considered and rejected was to assume a saturated log-density of zero:
it is right for a Bernoulli trial, wrong for a Poisson count, and
nothing visible from the core tells the two apart, so the assumption
would have produced a plausible wrong number in silence.

### The block's `group`

A fourth reserved block name, beside `y`, `miss` and `mask`. A factor
or integer `n`-vector giving the family's own independent unit for each
row. It is what `loglik_group` returns one value per, in level order
for a factor and `sort(unique())` order otherwise; the order is fixed
by the core rather than by each family, because the core aligns two
groupings by these codes.

The GROUPING-ALIGNMENT check is required, because getting it wrong is
silent. The family's units and the correction's grouping levels must be
the same partition of the rows: `rw_delta(subject = id)` under
`(1 | id)` is the aligned case, and the same family under `(1 | item)`
gives a per-subject likelihood against a per-item proposal, which does
not add up in either direction. Equality, not refinement: a family unit
spanning two groups has no separable integrand, and a group holding two
family units would have the family concatenate two sequences that never
met. `check_importance_scope()` refuses by name, in groups and levels
rather than in rows.

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
structured families' own for the tokens `hmm`, `lca`, `mix_g`,
`hmm_g` (word boundaries, outside roxygen text), and fail on a hit.

It could not start at zero, because the branches it polices existed
until steps 6 and 7 deleted them. So it started as a RATCHET: the test
pinned the current file-by-token hit inventory exactly, failed on any
NEW hit, and each protocol step that deleted a branch shrank the
pinned inventory in the same commit. The exempt list outlived the
inventory by three steps, because `R/hmm.R` and `R/lca.R` were
themselves in `R/` until step 10 moved them out.

As built, the test is now the boundary rather than a countdown to it:
every token maps to an empty list of allowed homes, so a hit anywhere
in core `R/` fails. Its positive control moved with the families. A
scanner that matched nothing anywhere would pass the boundary
vacuously, and the control used to prove otherwise by reading the
exempt homes inside core; it now reads
`extensions/frmtmb.latent/R`, and it is the CONTROL that skips when
`extensions/` is absent, never the boundary assertion.

`mixture` tokens
are not policed: mixture is the in-core reference implementation and
stays.

## Packaging (decided 2026-09-03)

The split, when it goes through, is a MONOREPO: one git repository,
one directory per package. Precedent is kaskr's own layout (RTMB and
TMB each live as a subdirectory of a monorepo). The mechanics are not
onerous for anyone:

- CRAN users are unaffected: `install.packages()` knows nothing about
  repository layout.
- Development installs are `remotes::install_github("aforren1/frmtmb",
  subdir = "<pkg>")`, one documented line.
- r-universe builds every package in a monorepo natively, so the
  Additional_repositories story for off-CRAN pieces is unchanged.
- CI partitions by path filters; the test files already map one to
  one, so the suites split for free.

Package names as built: `frmtmb` the core, `frmtmb.ode`,
`frmtmb.sample`, and `frmtmb.latent` for the two structured families.
They took a package of their own rather than riding in
`frmtmb.sample`, because the sampling surface is a dependency they do
not need and their reference packages are not its reference packages.
`frmtmb.latent` is still a working name: nothing outside the monorepo
depends on it yet, so renaming it costs one DESCRIPTION, one workflow
file, the pointers in the core README and case-studies vignette, and
the extension path in the boundary test's positive control.

Tier context, from the split memo this protocol serves:

- Tier 1, ODE, splits first: two hook calls from frame.R, an off-CRAN
  hard dependency, and its own vignette. Removing it also removes
  Additional_repositories from the core DESCRIPTION.
- Tier 2, the draws surface (frm_sample, the draws methods, loo, the
  non-centered reparameterization, sampling default priors), splits
  in the same round as the ODE extraction (maintainer decision
  2026-09-03, revised from "at CRAN time": the sampling tests are the
  heaviest part of the suite, so the split buys iteration speed
  immediately; see dev/ode-extraction.md for the round's coordination
  rules). Prerequisite, DONE at v0.45.0: the CORE prior machinery in
  interop.R (resolve_prior_input, neg_log_prior_fn, resolve_bounds,
  the prior constructors) moved to priors.R; the sampling-only
  default-prior machinery deliberately stayed, because it leaves with
  this package.
  Note two corrections to the memo: `ncp_plan` and its family are
  sampling-only, not core (the fit never touches them), and
  `check_laplace()` ships WITH the draws package (it already needs
  tmbstan), which costs core its own approximation-verifier - core's
  docs then point at the draws package for it, rather than staying
  silent.
- Tier 3, hmm and lca, waited on this protocol and left with step 10
  for `frmtmb.latent`. Nothing in core names them any more.

## Order of work

1. [x] `frmtmb_structure()` constructor and validator.
   `frmtmb_family()` gains `structure =`. No behavior change.
2. [x] Frame: `blocks` slot, `frame_vars`, `keep_na`, `check_spec`,
   `frame_block`, `check_frame` call sites. Old slots still populated.
3. [x] Objective: the single `loglik` branch, with mixture and hmm
   forwarding to it. Old branches deleted.
4. [x] Predict, fitted, residuals, simulate, conditional effects: the
   `supports` checks and `fitted_mean` / `fitted_var` / `block$miss`.
5. [x] `mixture()` onto the protocol. `test-lca.R`,
   `test-mvn-mixture.R` and the mixture tests in `test-families.R`
   must pass unchanged. Delete `mix_g`.
6. [x] `hmm()` onto the protocol. `test-hmm.R` unchanged. Delete
   `hmm_g` and every hmm string in core. Boundary test goes green.
7. [x] `lca()` onto the protocol; drop its borrowed `fam$mix`.
8. [x] `latent_probs()` generic; aliases.
9. [x] Export the accessor list above, document `frmtmb_structure()`
   as the extension API with hmm as the worked example.
10. [x] Move hmm.R, lca.R and their tests to the extension package.
    Landed as `extensions/frmtmb.latent` (one package for both
    latent-state families). Four accessors joined the extension API to
    make it compile out of tree: `frame_block_of()`,
    `structure_supports_all()`, `mixture_posterior()` and
    `mixture_multimodal_refusals()`. Step 9's list was otherwise
    complete, and two of its guesses were wrong in the useful
    direction: neither `linpred_key()` nor `assemble_frame()` is
    reached at all. The written inventory is
    `extensions/frmtmb.latent/dev/out-of-tree-inventory.md`, and the
    extension's own suite re-runs the scan so a later edit cannot
    quietly reach for an internal again.

Four things landed with steps 6 through 9 that this document did not
anticipate, and its text above is the design as approved rather than as
built. The differences are: a `check_fit` slot, because hmm's
label-symmetry start warning is not a capability question and had
nowhere else to go; a `unit` slot and a `cluster_robust` flag, which
are what loo.R and sandwich.R read instead of naming a family; and
`loglik` becoming optional, which is what makes a structure usable as a
pure capability declaration by a family whose likelihood IS rowwise
(lca, and mixture() without groups). `fam$mix` stayed on `lca()`: it is
not a marker the core reads but `mixture()`'s component interface, and
it is why one posterior implementation serves all three mixture-type
families.

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

## The addition-term allow-list (added 2026-09-05)

Not a `frmtmb_structure()` slot, and it belongs here because it answers
the same question one file over: what a family declares about itself
instead of hand-writing a check.

`frmtmb_family(required_aterms =)` is a CONJUNCTION of what the density
cannot do without. There was no complementary allow-list, so a term the
density never reads was parsed, coerced, stored on the fitted object,
and then ignored. Measured, before the fix: `lba(3)` on
`rt | dec(two) + vint(choice) ~ 1` fitted with fixed effects
bit-identical to the model without `dec()` and no warning at any
point. Each family that wanted to refuse such a term had to write the
check itself, and the two sibling race families disagreed about it
for a release.

`wiener()` has the same defect and the allow-list does NOT close it,
which is worth stating here rather than leaving the doc claiming a
win it did not get. The earlier draft of this paragraph said
`wiener()` took a `vint()` it cannot use; that is wrong, and the
review measured it: `wiener()` reads `vint1` as the boundary when
`dec` is absent, so `rt | vint(upper)` and `rt | dec(upper)` give a
bit-identical logLik of -130.566406836. `vint()` is a working second
spelling. The real defect is the two supplied TOGETHER:
`rt | dec(upper) + vint(upper)` fits, `dec` wins, the `vint()` column
travels into the fit and changes nothing, silently. An allow-list
cannot catch that, because both terms are legitimately on the list;
it needs exclusivity among the alternatives of `required_aterms`,
which is the follow-up recorded under "left out" in
`dev/protocol-findings.md` and in the note at the end of this
section.

`accepts_aterms` is that allow-list. It names the terms a family reads
or lets the core act on, as a formula writes them and without
parentheses, and is read together with `required_aterms` (a required
term need not be repeated). `NULL`, the default, accepts every
registered term, so a custom family written before this argument
existed keeps its behavior. `character(0)` declares a family that takes
none.

Two decisions worth keeping:

- The check runs LAST of the addition-term guards, after `valid_y()`.
  Every earlier refusal says something specific about the pair it
  refuses (`se()` for a non-gaussian family, `cens()` without a CDF,
  `wiener_gng`'s "the indicator travels through dec(), not vint()"),
  and those messages are better than a generic one. What the allow-list
  catches is what nothing else did.
- It is spelled in TERMS, not in term values: `vint`, not `vint1` and
  `vint2`. So it closes "this family never reads that term" and does
  not close "this family reads one `vint()` value and you supplied
  two". The `wiener()` case the rdm/gng review found is the second
  kind - `dec()` and `vint()` are alternative spellings of one datum
  there, and supplying both leaves the second unread - and it needs an
  exclusivity rule on `required_aterms` alternative groups instead.
  `gddm()` is why that rule is not written here: it reads `dec()` and
  `vint()` together, with the condition index in `vint1` when `dec()`
  is present and in `vint2` when it is not.
