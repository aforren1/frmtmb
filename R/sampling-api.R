# The core internals a sampling package needs, made public.
#
# The sampling surface (frm_sample(), check_laplace(), the frmtmb_draws
# methods, log_lik() and the leave-one-out estimators) lives in the
# frmtmb.sample package. dev/draws-extraction.md is the call graph it
# was taken from and carries one resolution per dependency; this file is
# the (b) column of that table, gathered in one place so the contract
# can be read as a whole rather than hunted through ten files.
#
# The exports are declared here with @rawNamespace rather than with an
# @export beside each definition, on purpose: what makes this API is the
# LIST, and a list that lives in one file can be reviewed, diffed and
# argued with. A tag scattered over ten files cannot.

#' Core internals for a sampling extension
#'
#' @description
#' These functions are the interface `frmtmb.sample` is written
#' against. They were internal until the sampling surface moved to its
#' own package, and they are documented here as one contract rather
#' than one page each, because what has to be reviewed is the set.
#'
#' They are lower-level than the rest of frmtmb: they take assembled
#' frames, resolved prior entries and single random-effect blocks, not
#' formulas. Most users want [frm()], [set_prior()] and
#' `frmtmb.sample::frm_sample()` instead. Nothing here validates its
#' arguments the way the user-facing surface does.
#'
#' @section The objective seam:
#' `build_objective(frame)` returns the bare negative log-likelihood
#' closure of an assembled frame, ready for `RTMB::MakeADFun()`. It
#' honors two fields a caller may set on a COPY of the frame:
#' `frame$map`, and `frame$ncp_blocks`, the indices of the
#' random-effect blocks to build in their non-centered form
#' (`b = L(theta) z`). A prior-augmented objective is this closure plus
#' `neg_log_prior_fn()`'s.
#'
#' `row_lpdf(fam, y, yobs, dpars, aterms, extra)` is the per-row
#' log-density composition the objective itself runs, with `cens()` and
#' `trunc()` folded in; it runs on numeric dpar values as readily as on
#' the tape, which is what makes a pointwise `log_lik()` reproduce the
#' fitted density exactly instead of approximating it.
#' `with_cs_offsets(fit, rspec, dpv)` takes one response spec from
#' `fit$spec$responses` and the dpar-value list `eval_dpars()` returns
#' for that response, and gives back the same list with the
#' category-specific (`cs()`) offsets applied; on a model without
#' `cs()` terms it returns `dpv` unchanged, so it is safe to call
#' unconditionally before `row_lpdf()`. `us_chol_cor(theta, K)` is the
#' unstructured correlation matrix of a `thetar` segment, which a
#' `set_rescor(TRUE)` model's joint row density needs.
#'
#' @section The prior seam:
#' The prior VOCABULARY - [set_prior()], [prior_normal()] and its
#' relatives, [get_prior()], [prior_summary()] - is ordinary exported
#' API and is not part of this page. What is here is the RESOLUTION
#' machinery underneath it: `as_priorlist()` coerces the accepted
#' spellings (including a brms `brmsprior`) to one,
#' `resolve_prior_input()` turns a priorlist plus a model into
#' `list(entries, lower, upper)` on the internal parameter scale,
#' `neg_log_prior_fn()` turns resolved entries into a tapeable closure,
#' `resolve_bounds()` turns user-spelled bounds into internal-scale
#' vectors over the outer parameters, and `spec_target()` names the slot
#' one specification addresses.
#'
#' `frmtmb_register_prior_defaults()` is the other direction: it lets a
#' package tell [get_prior()] what defaults it would apply, so that the
#' reported default is true of the routes the session actually has.
#'
#' @section The covstruct block readers:
#' Each takes ONE element of `frame$re_blocks` and answers one
#' structural question about it. `ncp_eligible()` says whether the
#' block has a linear Cholesky factor and only standard-deviation and
#' correlation parameters, so that `b = L(theta) z` is a bijection;
#' `ncp_scale_b()` and `ncp_unscale_b()` are that map and its inverse.
#' `covstruct_has_chol()` is the first half of eligibility on its own
#' (a factor is registered at all), `block_sd_idx()` gives the block's
#' standard-deviation positions within its theta segment,
#' `block_n_cor()` counts its correlation positions, `block_cor_prior()`
#' says whether an LKJ density fits them (`"none"`, `"lkj"` or
#' `"unsupported"`), and `is_student_block()` reports a Student-t
#' latent, which is a scale mixture and so has no linear factor at all.
#'
#' The registry these read is deliberately NOT exported: exporting it
#' would make every field of every structure's entry into API, and
#' these eight questions are the ones anything outside covstruct.R has
#' ever needed to ask.
#'
#' @section The simulator seam:
#' `sim_can()` and `sim_note()` say whether a family can be simulated
#' from and what is missing when it cannot; `sim_context()` builds the
#' context one draw is taken from and `sim_draw()` takes it.
#' `sim_is_structured()` reports whether the draw walks a fitted
#' structure - a hidden state sequence, a group-level latent class, a
#' correlated residual - which is what makes `newdata` and `re_formula`
#' refusable rather than merely unimplemented on a predictive method.
#'
#' @section Parameter labeling and fitted quantities:
#' `par_name_bare()` is the draws-side spelling of a parameter name,
#' with parentheses dropped (`Intercept`, not `(Intercept)`);
#' `outer_par_names()` names the outer parameter vector in
#' [confint()] row order; `estimated_coef_names()` names the estimated
#' coefficient vector with mapped `betad` entries excluded;
#' `log_sd_theta_index()` says which theta positions are log standard
#' deviations, so that a variance at its boundary can be told from an
#' AR(1) coefficient at its own. `sdr_of()` is the cached
#' `TMB::sdreport()`, and `require_fitted()` is the refusal a
#' maximum-likelihood-only method owes an object that was assembled but
#' never optimized.
#'
#' @section The hypothesis expression engine:
#' `hyp_parse_all()` turns hypothesis strings into expressions and
#' per-expression directions, once. `hyp_vals_only(fit)` reads the
#' fit's current estimates into a flat named value vector plus a
#' parallel component vector (which template component each value came
#' from); `hyp_env_vals(fit, vals, comp)` takes exactly that pair and
#' builds the evaluation environment, and `hyp_eval()` evaluates one
#' expression in it - so a caller with many parameter vectors parses
#' once, then per vector swaps the estimates in and calls
#' `hyp_vals_only()` and `hyp_env_vals()` again.
#' `hyp_tail_p()` is the tail probability of a directional claim.
#'
#' @section The conditional-effects engine:
#' `ce_grids_build()` builds the prediction grids, effect list,
#' condition sets and base values; `ce_boot_one()` evaluates one grid
#' at one parameter vector and flattens it; `ce_finalize()` assembles
#' the per-effect data frames into the returned object with the
#' attributes `plot()` reads. `ce_cats_display()` says whether the
#' display is per-category, `ce_structure_check()` is the refusal a
#' structured likelihood owes a grid, and `ce_re_formula()` resolves
#' the random-effect argument of such a call.
#'
#' @section The two-dialect argument seam:
#' frmtmb answers to two argument dialects: a brms-named function takes
#' `re_formula`, frmtmb's own fit surface takes lme4's `re.form`, and
#' the brms-named ones accept both. `arg_unset()` is the "not supplied"
#' marker a formal defaults to when `NULL` and `NA` are both real
#' settings and neither can double as unset; `re_form_arg()` resolves
#' the pair, refusing rather than guessing when both are given.
#'
#' @param frame An assembled model frame (`fit$frame`).
#' @param bk One element of `frame$re_blocks`.
#' @param fit A `frmtmb_fit`.
#' @param ... Arguments of the individual functions; see the sections
#'   above and the source, which is the reference for these.
#' @return As described per function above.
#' @seealso [frmtmb-extension-api] for the family-level accessors a
#'   structured family uses, [frmtmb_structure()] for the protocol
#'   those serve, and [frmtmb_register_compat()] for the compatibility
#'   registry, which was documented here while it had no page of its
#'   own.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
#' dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # the objective seam: the bare likelihood closure of the frame
#' nll <- build_objective(fit$frame)
#' nll(fit$estimates)
#'
#' # a block reader: this one has a diagonal factor and one sd, so it
#' # can be sampled non-centered
#' bk <- fit$frame$re_blocks[[1]]
#' c(eligible = ncp_eligible(bk), n_cor = block_n_cor(bk),
#'   cor_prior = block_cor_prior(bk))
#'
#' # parameter labeling, in the two spellings
#' head(outer_par_names(fit))
#' head(par_name_bare(outer_par_names(fit)))
#' @name frmtmb-sampling-api
# One @aliases tag per alias, not one wrapped tag: roxygen2 merges
# repeated tags into the same alias set, and pkgcheck's roxygen parse
# rejects a multi-line @aliases outright where roxygen2 only warns.
#' @aliases frmtmb-sampling-api
#' @aliases build_objective
#' @aliases row_lpdf
#' @aliases with_cs_offsets
#' @aliases us_chol_cor
#' @aliases aterms_for_newdata
#' @aliases has_trunc
#' @aliases as_priorlist
#' @aliases resolve_prior_input
#' @aliases neg_log_prior_fn
#' @aliases resolve_bounds
#' @aliases spec_target
#' @aliases frmtmb_register_prior_defaults
#' @aliases ncp_eligible
#' @aliases ncp_scale_b
#' @aliases ncp_unscale_b
#' @aliases covstruct_has_chol
#' @aliases block_sd_idx
#' @aliases block_cor_prior
#' @aliases block_n_cor
#' @aliases is_student_block
#' @aliases sim_can
#' @aliases sim_note
#' @aliases sim_context
#' @aliases sim_draw
#' @aliases sim_is_structured
#' @aliases par_name_bare
#' @aliases outer_par_names
#' @aliases estimated_coef_names
#' @aliases log_sd_theta_index
#' @aliases sdr_of
#' @aliases require_fitted
#' @aliases hyp_parse_all
#' @aliases hyp_vals_only
#' @aliases hyp_env_vals
#' @aliases hyp_eval
#' @aliases hyp_tail_p
#' @aliases ce_grids_build
#' @aliases ce_boot_one
#' @aliases ce_finalize
#' @aliases ce_cats_display
#' @aliases ce_structure_check
#' @aliases ce_re_formula
#' @aliases find_linpred
#' @aliases arg_unset
#' @aliases re_form_arg
#' @rawNamespace export(build_objective, row_lpdf, with_cs_offsets,
#'   us_chol_cor, aterms_for_newdata, has_trunc, as_priorlist,
#'   resolve_prior_input, neg_log_prior_fn, resolve_bounds, spec_target,
#'   frmtmb_register_prior_defaults, ncp_eligible, ncp_scale_b,
#'   ncp_unscale_b, covstruct_has_chol, block_sd_idx, block_cor_prior,
#'   block_n_cor, is_student_block, sim_can, sim_note, sim_context,
#'   sim_draw, sim_is_structured, par_name_bare, outer_par_names,
#'   estimated_coef_names, log_sd_theta_index, sdr_of, require_fitted,
#'   hyp_parse_all, hyp_vals_only, hyp_env_vals, hyp_eval, hyp_tail_p,
#'   ce_grids_build, ce_boot_one, ce_finalize, ce_cats_display,
#'   ce_structure_check, ce_re_formula, find_linpred, arg_unset,
#'   re_form_arg)
NULL

# ---- the prior-defaults registry -------------------------------------
#
# get_prior() reports the default in every targetable slot. Core's own
# default is flat everywhere - frm() is maximum likelihood until a prior
# is set - but frmtmb.sample's frm_sample() applies brms's
# weakly-informative defaults on both of its routes. After the split
# core cannot state those and must not pretend to, so the package that
# owns them registers a provider and get_prior() consults it.
#
# The same shape as the compat registry, and for the same reason: an
# environment filled from a contributor's .onLoad(), by which time every
# namespace is sealed and the collation-order question does not arise.
# With no provider registered, every row reads "(flat)", which is
# exactly what frm() does.

frmtmb_prior_defaults <- new.env(parent = emptyenv())
frmtmb_prior_defaults$providers <- list()

#' Register a source of default priors with [get_prior()].
#'
#' `provider(spec, frame)` is called with a parsed specification and an
#' assembled frame and returns a `frmtmb_priorlist` of the defaults it
#' would apply to that model, or `NULL`. [get_prior()] fills a row's
#' `prior` column from the registered providers when one of them
#' addresses that row's slot, and leaves it `(flat)` otherwise.
#'
#' @noRd
frmtmb_register_prior_defaults <- function(provider) {
  stopifnot(is.function(provider))
  frmtmb_prior_defaults$providers <-
    c(frmtmb_prior_defaults$providers, list(provider))
  invisible(NULL)
}

#' The slot key a default and a `get_prior()` row are matched on: the
#' class, plus every qualifier that makes one slot a different slot from
#' another of the same class.
#'
#' @noRd
prior_slot_key <- function(class, dpar = "", nlpar = "", resp = "") {
  paste0(class, "|", dpar, "|", nlpar, "|", resp)
}

#' The registered defaults for one model, as a key-to-printed-density
#' lookup. Empty when nothing is registered, which is the core-alone
#' case and the reason every row then reads "(flat)".
#'
#' @noRd
registered_prior_defaults <- function(spec, frame) {
  out <- list()
  for (p in frmtmb_prior_defaults$providers) {
    pl <- tryCatch(p(spec, frame), error = function(e) NULL)
    for (s in unclass(pl %||% list())) {
      k <- prior_slot_key(s$class, s$dpar %||% "", s$nlpar %||% "",
                          s$resp %||% "")
      out[[k]] <- format_prior_dist(s$dist)
    }
  }
  out
}

#' One prior density as `get_prior()` prints it.
#'
#' @noRd
format_prior_dist <- function(d) {
  if (is.null(d)) return(NA_character_)
  kind <- if (identical(d$kind, "t")) "student_t" else d$kind
  paste0(kind, "(", paste(unlist(d[-1L]), collapse = ", "), ")")
}
