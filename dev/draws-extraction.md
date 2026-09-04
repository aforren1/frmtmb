# Extracting the sampling surface: frmtmb.sample

Status: inventory written 2026-09-03 against v0.46.0, before any code
moved. Sequencing and the coordination rules are in
`dev/ode-extraction.md`, section "The draws package". This document is
deliverable 1: the call graph, and one resolution per core internal the
sampling surface reaches. The move follows it.

## What "the sampling surface" is

The roots of the graph, as a set of top-level objects:

- `R/interop.R`, everything from `clamp_into_bounds()` to
  `check_tmbstan_build()`: `frm_sample()`, `as_tmbstan()`,
  `check_laplace()`, the mode inits and bounds, the non-centering plan
  and its transforms, the `frm_sample()` default priors,
  `sample_resolve_priors()`, the tmbstan build guard, and the
  `as.matrix()`/`print()` methods on `frmtmb_draws`. `interop.R`'s
  other three tenants stay: `check_custom_family()`, the
  emmeans/marginaleffects glue, and `getME.frmtmb_fit()`.
- `R/methods-draws.R` wholesale.
- `R/loo.R` wholesale.
- `conditional_effects.frmtmb_draws()` in `R/conditional-effects.R`.

`pp_check.frmtmb_draws()` and `prior_summary.frmtmb_draws()` are in
`methods-draws.R` and travel with it; their generics and their
`frmtmb_fit` methods stay in core.

## How the graph was taken

`codetools::findGlobals()` over the installed namespace, transitively
closed from those roots, intersected with the namespace's own objects.
The reverse direction was taken as well: for each internal reached,
which NON-moving core objects also reach it. An internal with no core
user is a candidate to leave rather than to be exported, and that
distinction is what keeps the export list from being the whole graph.
The script is `dev/` scratch, not committed; the two numbers it
produced are reproduced by re-running `findGlobals()` and are not
quoted here as counts.

## Resolutions

Three resolutions, as agreed:

- **(a)** already exported, nothing to do;
- **(b)** a NEW core export with a documented contract, grouped on a
  new `frmtmb-sampling-api` help page;
- **(c)** restructured so the dependency disappears - the extension
  carries its own copy of a context-free helper, or reads a documented
  field directly.

### (a) Already exported

`frm()`, `frmtmb_control()`, `set_prior()`, `ranef()`, `VarCorr()`,
`ngrps()`, `mixture_probs()`, `eval_dpars()`, `fit_extras()`,
`single_response()`. The last three are the structured-family
extension API, already public for the protocol the ODE and
structured-family lanes use; the sampling surface is their second
customer, which is the evidence they were the right shape.

### (b) New core exports

Grouped by the question each answers. Every one is reached from the
sampling surface and has at least one core caller as well, so none of
them is API invented for one consumer.

**The objective seam.** The sampling surface rebuilds the objective
three times over: with prior terms added, with blocks non-centered,
and with a MAP fit's own penalty removed.

| export | file | contract |
|---|---|---|
| `build_objective(frame)` | objective.R | The bare negative log-likelihood closure of an assembled frame. Honors `frame$ncp_blocks` (the non-centered block indices) and `frame$map`. Core's own fit path and cluster scores call it; the sampling package calls it to add priors or to non-center. |
| `row_lpdf(fam, y, yobs, dpars, aterms, extra)` | objective.R | One log-density per row, with `cens()` and `trunc()` folded in. The composition `log_lik()` must reproduce exactly, on numeric dpars instead of on the tape. |
| `with_cs_offsets(fit, newdata, dpars)` | predict.R | Category-specific offsets applied to evaluated dpars. |
| `us_chol_cor(theta, K)` | covstruct.R | The unstructured correlation matrix of a `thetar` segment. `log_lik()` needs it for a `set_rescor(TRUE)` model. |

**The prior seam.** The prior VOCABULARY (`set_prior()`,
`prior_normal()` and relatives, `prior_summary()`, `get_prior()`) is
already exported and stays in core. What was missing is the
RESOLUTION machinery, which turns a specification into taped terms.

| export | file | contract |
|---|---|---|
| `as_priorlist(x)` | priors.R | Coerce a `set_prior()` spec, a `brmsprior`, or a legacy named list to the internal spelling, or pass a named list through. |
| `resolve_prior_input(fit, pl)` | priors.R | A priorlist against a model, giving `list(entries, lower, upper)` on the internal parameter scale. |
| `neg_log_prior_fn(entries)` | priors.R | The tapeable negative-log-prior closure of resolved entries. |
| `resolve_bounds(fit, lower, upper)` | priors.R | User-spelled bounds to internal-scale `list(lower, upper)` over the outer parameters. |
| `spec_target(s)` | priors.R | The slot key one specification addresses, for supersession and for the message that names a dropped prior. |

**The covstruct block readers.** One block of `frame$re_blocks` in,
one structural fact out. These are what the non-centering plan and the
default priors ask, and they are the reason `covstruct_registry` and
`lkj_refusals` do NOT need to be exported: three new small readers
replace five direct registry and refusal-list peeks.

| export | file | contract |
|---|---|---|
| `ncp_eligible(bk)` | covstruct.R | Whether the block has a linear Cholesky factor and only sd/correlation parameters, so `b = L(theta) z` is a bijection. |
| `ncp_scale_b(bk, z, theta)` | covstruct.R | `b = L(theta) z`. Core's `build_objective()` already calls it for a non-centered frame. |
| `ncp_unscale_b(bk, b, theta)` | covstruct.R | `z = L(theta)^-1 b`, for mapping the ML mode onto the non-centered scale. |
| `covstruct_has_chol(bk)` | covstruct.R | NEW. Whether any Cholesky factor is registered for the block's structure. Replaces a `covstruct_registry[[cs]][["chol_sd"]]` peek. |
| `block_sd_idx(bk)` | covstruct.R | NEW. The block's theta positions that are standard deviations. Replaces a `covstruct_registry[[...]]$sd_idx(bk$dim)` peek. |
| `block_cor_prior(bk)` | covstruct.R | NEW. What correlation prior the block can carry: `"none"` (it has none), `"lkj"` (a density fits), `"unsupported"` (it has correlations but its parameterization has no correlation matrix over the whole of it). Three-valued because the two callers want different answers and they are not complements: choosing a default wants `"lkj"`, listing the slots left flat wants `"unsupported"` alone. Replaces two `names(lkj_refusals)` peeks and one `block_cor_spec()` reach. |
| `block_n_cor(bk)` | covstruct.R | How many correlation parameters the block has. |
| `is_student_block(bk)` | covstruct.R | Whether the latent is Student-t (a scale mixture, so not a linear factor). |

**The simulator seam**, for `posterior_predict()`.

| export | file | contract |
|---|---|---|
| `sim_can(fam)` / `sim_note(fam)` | families.R | Whether a family has a simulator, and the sentence naming what is missing. |
| `sim_context(fit, rspec, dpars, aterms, n, extra)` | families.R | The context object a draw is taken from. |
| `sim_draw(ctx)` | families.R | One simulated response from a context. |
| `sim_is_structured(ctx)` | families.R | Whether the draw walks a fitted structure (a state sequence, a group latent class, a correlated residual), which is what makes `newdata` and `re_formula` refusable. |

**Parameter labeling and the fitted-quantity readers.**

| export | file | contract |
|---|---|---|
| `par_name_bare(x)` | confint.R | The draws-side spelling of a parameter name: parentheses dropped, `Intercept` not `(Intercept)`. |
| `outer_par_names(fit)` | confint.R | Names of the outer parameter vector, in `confint()` row order. |
| `estimated_coef_names(fit)` | methods-fit.R | Names of the estimated coefficient vector, mapped `betad` entries excluded. |
| `log_sd_theta_index(fit)` | confint.R | Which theta positions are log standard deviations, so a boundary variance can be told from an AR(1) phi at its own boundary. |
| `sdr_of(fit)` | fit.R | The cached `sdreport()`. `check_laplace()` reads `cov.fixed` from it. |
| `require_fitted(fit, what)` | fit.R | The refusal a maximum-likelihood-only method owes an unfitted object. |

`all_par_labels()` is NOT here: it lives in `interop.R`, has no core
caller, and leaves with the extension.

**The hypothesis expression engine**, for
`hypothesis.frmtmb_draws()`. Five exports because the draws method
parses ONCE and evaluates PER DRAW, and core's fit method (which also
computes a Wald standard error the draws method must not) cannot be
reused whole.

| export | file | contract |
|---|---|---|
| `hyp_parse_all(hypothesis, names, class, group)` | confint.R | Hypothesis strings to expressions plus per-expression direction. |
| `hyp_vals_only(fit)` | confint.R | `list(vals, comp)`: the parameter values and their component map at one parameter vector. |
| `hyp_env_vals(fit, vals, comp)` | confint.R | The evaluation environment; its NAMES are the vocabulary `hyp_parse_all()` checks against. |
| `hyp_eval(fit, expr, vals, comp)` | confint.R | One expression at one parameter vector. |
| `hyp_tail_p(t, dir)` | confint.R | The posterior tail probability of a directional claim. |

**The conditional-effects engine**, for
`conditional_effects.frmtmb_draws()`. The stages core's own
bootstrap band already runs per refit, run here per draw.

| export | file | contract |
|---|---|---|
| `ce_grids_build(fit, rspec, lp, effects, resp, dpar, resolution, conditions, data)` | conditional-effects.R | The prediction grids, the effect list, the condition sets and the base values. |
| `ce_boot_one(fit, nd, categorical, resp, dpar, re_formula)` | conditional-effects.R | One grid evaluated at one parameter vector, flattened. |
| `ce_finalize(dfs_by_eff, effects, rspec, resp, dpar, band, base, categorical, cond_sets, groups)` | conditional-effects.R | The per-effect data frames to the returned object, with the attributes `plot()` reads. |
| `ce_cats_display(rspec, dpar)` | conditional-effects.R | Whether the display is per-category. |
| `ce_structure_check(rspec)` | conditional-effects.R | The refusal a structured likelihood owes a grid. |
| `ce_re_formula(re_formula, dots)` | conditional-effects.R | The `re_formula`/`re.form` resolution for this call. |

**The two-dialect argument seam.**

| export | file | contract |
|---|---|---|
| `arg_unset()` | utils.R | The "not supplied" marker a formal defaults to when `NULL` and `NA` are both real settings. |
| `re_form_arg(re_formula, re.form, caller)` | utils.R | One of the two spellings, or the refusal when both are given. |

**One more linear-predictor reader.**

| export | file | contract |
|---|---|---|
| `find_linpred(fit, resp, dpar)` | sugar.R | One linear predictor of the frame, with the multivariate "disambiguate with resp =" refusal. |

**The two registration seams.** `frmtmb_register_compat()` and
`compat_rule_builder()` were written for exactly this and marked
"Not exported yet: the two contributors are still in this package.
Moving them out is what makes it public API." A contributor now lives
outside the package, so they are exported, unchanged. The third,
`frmtmb_register_prior_defaults()`, is new and is described below.

### (c) Restructured away

| internal | file | resolution |
|---|---|---|
| `%||%` | utils.R | The extension defines its own. One line, and a universal idiom. |
| `check_count`, `check_named_list`, `check_probability` | utils.R | The extension carries its own copies. Context-free argument checks that know nothing about frmtmb; exporting them would make three messages into API for no gain. |
| `ce_pctl` | conditional-effects.R | Same: a five-line pointwise percentile. Core keeps its copy for the bootstrap band. |
| `fam_structure`, `structure_unit` | structure.R | The extension reads `fam[["structure"]]` and its `unit` field directly. Both are documented fields of `frmtmb_structure()`, so this is reading public protocol, not reaching into core. |
| `covstruct_registry`, `lkj_refusals` | covstruct.R | Replaced by the three new block readers above. Exporting the registry itself would have made every field of every structure's entry into API. |
| `refuse_retired_priors` | utils.R | No core caller: it refuses `frm_sample()`'s retired `priors =` spelling, which only that function ever accepted. Moves. |
| `VarCorr.frmtmb_fit` | methods-fit.R | Reached by name in `VarCorr.frmtmb_draws()`, which called the method directly to skip dispatch. The extension calls the exported generic `VarCorr()` instead. |

## The get_prior defaults seam

`get_prior()` prints `(flat)` in the `prior` column of every row,
which is true of `frm()` and false of `frm_sample()`. After the split
the sampling defaults live in the extension, so core cannot state them
and must not pretend to.

The seam is the compat pattern, at its smallest. Core gains a
one-function registry:

```r
frmtmb_register_prior_defaults(provider)
```

`provider(spec, frame)` returns a `frmtmb_priorlist` or `NULL`. Core
holds at most the providers registered with it, and `get_prior()`
consults them when filling the `prior` column: a row whose
class/dpar/nlpar/resp key matches a provided specification shows that
specification, and every other row still shows `(flat)`. With no
provider registered - core alone, extension absent - every row reads
`(flat)`, which is exactly what `frm()` does. With `frmtmb.sample`
loaded, `get_prior()` reports what `frm_sample()` would apply.

The column is therefore truthful in both worlds, which is the point,
and the registry is one environment and one lookup.

## The tmbstan refusal string

`interop.R`'s "the sampler returned no draws" message names `frm_ode`
and RTMBode. `dev/ode-extraction.md` asks that core stop naming the
extension, and prefers (a): a generic sentence, with the ODE package
documenting its own failure. That question resolves HERE, because the
whole message leaves core with `frm_sample()`.

The resolution taken: the message keeps a generic "a tape that calls
an external solver" sentence in `frmtmb.sample`, and does not name
`frm_ode` or RTMBode. Naming a package that `frmtmb.sample` does not
depend on, from a message a user sees only when they have already hit
the failure, buys a pointer that the ODE package's own documentation
gives better. Core is left naming neither.

## Test disposition, vignettes, NEWS

In the report, not here: this document is the export inventory, which
is the thing that has to be reviewed before the code moves.

## Consequence for the core-boundary ratchet

`tests/testthat/test-core-boundary.R` pinned three entries:
`interop.R|frm_ode`, `interop.R|RTMBode` and `methods-draws.R|lca`.
All three left with this extraction - the first two inside the tmbstan
refusal message, the third inside two `stop()` messages of
`methods-draws.R`, which is now a file of frmtmb.sample.

The pinned inventory is therefore EMPTY, which is the endgame
`dev/ode-extraction.md` predicted ("the inventory should then be empty
or compat-only, and the boundary test's comment flips from ratchet
language to zero-tolerance language"). The growth-only ratchet permits
a shrink, but the test does not survive an empty one: `boundary_scan()`
returns `integer(0)`, whose `names()` is `NULL`, and `order(NULL)`
errors before any assertion runs. Even past that, the deliberate
fail-closed guard `expect_gt(sum(found), 0L)` would fire, because a
scanner that matches nothing would otherwise pass vacuously.

This lane does NOT edit that file: re-pinning belongs to the
consolidator, after both extractions merge. What the consolidator needs
to do there is flip it to the zero-tolerance form - an empty `pinned`,
and the vacuity guard rewritten to assert that the scan found nothing
rather than that it found something.
