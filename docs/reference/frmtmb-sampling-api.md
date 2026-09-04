# Core internals for a sampling extension

These functions are the interface `frmtmb.sample` is written against.
They were internal until the sampling surface moved to its own package,
and they are documented here as one contract rather than one page each,
because what has to be reviewed is the set.

They are lower-level than the rest of frmtmb: they take assembled
frames, resolved prior entries and single random-effect blocks, not
formulas. Most users want
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md),
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
and
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.html)
instead. Nothing here validates its arguments the way the user-facing
surface does.

## Arguments

- frame:

  An assembled model frame (`fit$frame`).

- bk:

  One element of `frame$re_blocks`.

- fit:

  A `frmtmb_fit`.

- ...:

  Arguments of the individual functions; see the sections above and the
  source, which is the reference for these.

## Value

As described per function above.

## The objective seam

`build_objective(frame)` returns the bare negative log-likelihood
closure of an assembled frame, ready for `RTMB::MakeADFun()`. It honors
two fields a caller may set on a COPY of the frame: `frame$map`, and
`frame$ncp_blocks`, the indices of the random-effect blocks to build in
their non-centered form (`b = L(theta) z`). A prior-augmented objective
is this closure plus `neg_log_prior_fn()`'s.

`row_lpdf(fam, y, yobs, dpars, aterms, extra)` is the per-row
log-density composition the objective itself runs, with `cens()` and
[`trunc()`](https://rdrr.io/r/base/Round.html) folded in; it runs on
numeric dpar values as readily as on the tape, which is what makes a
pointwise `log_lik()` reproduce the fitted density exactly instead of
approximating it. `with_cs_offsets(fit, rspec, dpv)` takes one response
spec from `fit$spec$responses` and the dpar-value list
[`eval_dpars()`](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
returns for that response, and gives back the same list with the
category-specific (`cs()`) offsets applied; on a model without `cs()`
terms it returns `dpv` unchanged, so it is safe to call unconditionally
before `row_lpdf()`. `us_chol_cor(theta, K)` is the unstructured
correlation matrix of a `thetar` segment, which a `set_rescor(TRUE)`
model's joint row density needs.

## The prior seam

The prior VOCABULARY -
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md),
[`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
and its relatives,
[`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md),
[`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md) -
is ordinary exported API and is not part of this page. What is here is
the RESOLUTION machinery underneath it: `as_priorlist()` coerces the
accepted spellings (including a brms `brmsprior`) to one,
`resolve_prior_input()` turns a priorlist plus a model into
`list(entries, lower, upper)` on the internal parameter scale,
`neg_log_prior_fn()` turns resolved entries into a tapeable closure,
`resolve_bounds()` turns user-spelled bounds into internal-scale vectors
over the outer parameters, and `spec_target()` names the slot one
specification addresses.

`frmtmb_register_prior_defaults()` is the other direction: it lets a
package tell
[`get_prior()`](https://aforren1.github.io/frmtmb/reference/get_prior.md)
what defaults it would apply, so that the reported default is true of
the routes the session actually has.

## The covstruct block readers

Each takes ONE element of `frame$re_blocks` and answers one structural
question about it. `ncp_eligible()` says whether the block has a linear
Cholesky factor and only standard-deviation and correlation parameters,
so that `b = L(theta) z` is a bijection; `ncp_scale_b()` and
`ncp_unscale_b()` are that map and its inverse. `covstruct_has_chol()`
is the first half of eligibility on its own (a factor is registered at
all), `block_sd_idx()` gives the block's standard-deviation positions
within its theta segment, `block_n_cor()` counts its correlation
positions, `block_cor_prior()` says whether an LKJ density fits them
(`"none"`, `"lkj"` or `"unsupported"`), and `is_student_block()` reports
a Student-t latent, which is a scale mixture and so has no linear factor
at all.

The registry these read is deliberately NOT exported: exporting it would
make every field of every structure's entry into API, and these eight
questions are the ones anything outside covstruct.R has ever needed to
ask.

## The simulator seam

`sim_can()` and `sim_note()` say whether a family can be simulated from
and what is missing when it cannot; `sim_context()` builds the context
one draw is taken from and `sim_draw()` takes it. `sim_is_structured()`
reports whether the draw walks a fitted structure - a hidden state
sequence, a group-level latent class, a correlated residual - which is
what makes `newdata` and `re_formula` refusable rather than merely
unimplemented on a predictive method.

## Parameter labeling and fitted quantities

`par_name_bare()` is the draws-side spelling of a parameter name, with
parentheses dropped (`Intercept`, not `(Intercept)`);
`outer_par_names()` names the outer parameter vector in
[`confint()`](https://rdrr.io/r/stats/confint.html) row order;
`estimated_coef_names()` names the estimated coefficient vector with
mapped `betad` entries excluded; `log_sd_theta_index()` says which theta
positions are log standard deviations, so that a variance at its
boundary can be told from an AR(1) coefficient at its own. `sdr_of()` is
the cached
[`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html), and
`require_fitted()` is the refusal a maximum-likelihood-only method owes
an object that was assembled but never optimized.

## The hypothesis expression engine

`hyp_parse_all()` turns hypothesis strings into expressions and
per-expression directions, once. `hyp_vals_only(fit)` reads the fit's
current estimates into a flat named value vector plus a parallel
component vector (which template component each value came from);
`hyp_env_vals(fit, vals, comp)` takes exactly that pair and builds the
evaluation environment, and `hyp_eval()` evaluates one expression in
it - so a caller with many parameter vectors parses once, then per
vector swaps the estimates in and calls `hyp_vals_only()` and
`hyp_env_vals()` again. `hyp_tail_p()` is the tail probability of a
directional claim.

## The conditional-effects engine

`ce_grids_build()` builds the prediction grids, effect list, condition
sets and base values; `ce_boot_one()` evaluates one grid at one
parameter vector and flattens it; `ce_finalize()` assembles the
per-effect data frames into the returned object with the attributes
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) reads.
`ce_cats_display()` says whether the display is per-category,
`ce_structure_check()` is the refusal a structured likelihood owes a
grid, and `ce_re_formula()` resolves the random-effect argument of such
a call.

## The two-dialect argument seam

frmtmb answers to two argument dialects: a brms-named function takes
`re_formula`, frmtmb's own fit surface takes lme4's `re.form`, and the
brms-named ones accept both. `arg_unset()` is the "not supplied" marker
a formal defaults to when `NULL` and `NA` are both real settings and
neither can double as unset; `re_form_arg()` resolves the pair, refusing
rather than guessing when both are given.

## The compatibility registry

`frmtmb_register_compat()` contributes feature rows and compatibility
rules to
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
from another package's `.onLoad()`, and `compat_rule_builder()` is the
accumulator that makes a contributed rule read like a core one.
Contributions are appended, and rules of equal specificity resolve
later-wins, so a contributed rule may override a core default and a core
default can never silently override a contributed one.

## See also

[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the family-level accessors a structured family uses, and
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
for the protocol those serve.
[`vignette("compatibility")`](https://aforren1.github.io/frmtmb/articles/compatibility.md)
for
[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# the objective seam: the bare likelihood closure of the frame
nll <- build_objective(fit$frame)
nll(fit$estimates)
#> [1] 19.66275

# a block reader: this one has a diagonal factor and one sd, so it
# can be sampled non-centered
bk <- fit$frame$re_blocks[[1]]
c(eligible = ncp_eligible(bk), n_cor = block_n_cor(bk),
  cor_prior = block_cor_prior(bk))
#>  eligible     n_cor cor_prior 
#>    "TRUE"       "0"    "none" 

# parameter labeling, in the two spellings
head(outer_par_names(fit))
#> [1] "(Intercept)"       "x"                 "sigma_(Intercept)"
#> [4] "theta_1"          
head(par_name_bare(outer_par_names(fit)))
#> [1] "Intercept"       "x"               "sigma_Intercept" "theta_1"        
```
