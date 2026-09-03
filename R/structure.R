# The structured-family protocol (dev/structured-family-protocol.md).
#
# WHY THIS EXISTS. Three families have a likelihood that does not
# factorize over rows, and each one used to reach the core through its
# own named family slot, its own frame slot and its own branch in seven
# core files. A fourth would have added a fourth column to that table.
# One slot (`fam$structure`), one constructor and a fixed set of call
# sites replace it, so a structured family can be written outside this
# package.
#
# `[[ ]]`, NEVER `$`, on a structure, a block or a supports vector. `$`
# falls back to partial matching when no exact name matches, which is
# how `ctx$mix` once silently read `ctx$mix_g`. Every read below is
# bracketed so that renaming a slot breaks loudly instead of quietly.
#
# EVERYTHING A STRUCTURE HOLDS IS EITHER DATA OR A PURE FUNCTION. The
# block a structure builds is saved inside the fit and rebuilt by
# refit(), so it may not capture the model frame; `loglik` runs on the
# AD tape, so it may not branch on a parameter.

#' The capability flags a structure may declare, and their defaults.
#'
#' All `FALSE`: a new structured family starts fully refused and opts
#' in, so a capability nobody thought about fails with a sentence rather
#' than producing a number nobody checked.
#'
#' @noRd
frmtmb_structure_flags <- c(
  reml = FALSE, quadrature = FALSE, profile = FALSE,
  newdata_response = FALSE, se_fit_response = FALSE, re_form = FALSE,
  conditional_effects = FALSE, osa = FALSE, deviance = FALSE,
  multivariate = FALSE, cens_trunc = FALSE, mi = FALSE
)

#' Assert that a slot is a function or absent.
#'
#' @noRd
check_structure_fn <- function(f, arg) {
  if (!is.null(f) && !is.function(f)) {
    stop("frmtmb_structure(", arg, " =) must be a function or NULL, not ",
         arg_desc(f), call. = FALSE)
  }
  invisible(f)
}

#' The flag part of a `refusals` name: `"re_form"` and
#' `"re_form.simulate"` both name the `re_form` flag.
#'
#' @noRd
refusal_flag <- function(nm) sub("\\..*$", "", nm)

#' Declare a non-rowwise likelihood to the core
#'
#' A `frmtmb_structure()` is what a family carries when its likelihood
#' does not factorize over the rows of the data: a group-level
#' [mixture()], a hidden Markov chain, a latent class measurement model.
#' It is one object with one contract, so the core needs no branch per
#' family and a structured family can live in another package. Attach it
#' with `frmtmb_family(structure = )`.
#'
#' Only `loglik` is required. Every other slot defaults to the rowwise
#' behavior, which is what a family that changes the likelihood and
#' nothing else needs.
#'
#' @section The block:
#' `frame_block()` returns the **block**: a plain list of DATA, stored
#' at `fit$frame$blocks[[response]]` and handed back to `loglik()`,
#' `fitted_mean()`, `fitted_var()`, `latent_probs()` and `sim_ctx()`. It
#' is saved inside the fit and rebuilt by `refit()`, so it must hold no
#' AD values and no closure that captures the model frame. Read and
#' write it with `[[ ]]` only: `$` partial matching is how a `mix` read
#' once returned `mix_g`.
#'
#' Three names in it are reserved, because the core reads them:
#' \describe{
#'   \item{`y`}{If present, replaces the response vector for every later
#'     stage, which is how a family fills a placeholder in for an `NA`
#'     it keeps. Otherwise the response is unchanged.}
#'   \item{`miss`}{A logical `n`-vector. Residuals are `NA` at these
#'     rows. Optional.}
#'   \item{`mask`}{A numeric 0/1 `n`-vector the family multiplies into
#'     its own density. The core does not read it; the name is reserved
#'     so that every structured family spells it the same way.}
#' }
#' Everything else in the block belongs to the family.
#'
#' @section Capability flags:
#' `supports` is a named logical vector or list. Every name defaults to
#' `FALSE`, so a structured family starts fully refused and opts in:
#' \describe{
#'   \item{`reml`}{`REML = TRUE`.}
#'   \item{`quadrature`}{`quadrature =` other than the Laplace default.}
#'   \item{`profile`}{`frmtmb_control(profile = TRUE)`.}
#'   \item{`newdata_response`}{`predict(newdata =, type = "response")`.}
#'   \item{`se_fit_response`}{`predict(se.fit = TRUE, type =
#'     "response")`.}
#'   \item{`re_form`}{`re.form =` in `predict()` and `simulate()`.}
#'   \item{`conditional_effects`}{[conditional_effects()].}
#'   \item{`osa`}{`residuals(type = "osa")`.}
#'   \item{`deviance`}{`residuals(type = "deviance")`.}
#'   \item{`multivariate`}{[mvbf()] and `rescor = TRUE`.}
#'   \item{`cens_trunc`}{`cens()` and `trunc()`.}
#'   \item{`mi`}{`mi()` on the same response.}
#' }
#' `FALSE` refuses the capability outright. `TRUE` means only that the
#' structure does not stand in the way: the core's ordinary rules still
#' apply, so a mixture that sets `deviance = TRUE` is still refused a
#' deviance residual by the family's missing `dev_fn`, with the message
#' that refusal has always carried.
#'
#' Prediction on the link scale with `dpar =` is always available and is
#' not a flag: the linear predictors are rowwise and belong to the core.
#'
#' `reml`, `quadrature` and `profile` are read at fit time; the rest by
#' the method each one names. `multivariate`, `cens_trunc` and `mi`
#' describe model shapes a family will usually want to refuse from
#' `check_spec` instead, where it can say which addition term was the
#' problem.
#'
#' `refusals[[flag]]` is the WHOLE message a user sees when that
#' capability is asked for, so a family explains its own refusal in its
#' own words instead of through a generic sentence. A flag with no
#' string of its own gets a generic one naming the family. A name may
#' carry a context suffix, `"re_form.simulate"`, when one flag is
#' refused for two different reasons in two places; the core falls back
#' to the bare flag name.
#'
#' @param frame_vars `function(fam)` returning a list of language
#'   objects whose variables must be in the model frame but belong to no
#'   linear predictor: a grouping column, a time column, a sequence id.
#'   It runs before any data is seen, so it may read only the family.
#' @param keep_na `TRUE` means an `NA` in the response is data the
#'   family reads, so the row survives `na.action`. `NA`s in every other
#'   variable still drop the row. A family that keeps them must handle
#'   them in `frame_block` (a mask, a placeholder) or in `loglik`.
#' @param check_spec `function(resp, spec, av)` run before the generic
#'   addition-term guards, so a structured family refuses a term in its
#'   own words rather than through a missing CDF further down. It sees
#'   the response spec, the whole spec (for the univariate check) and
#'   the evaluated addition terms. It returns nothing and stops to
#'   refuse.
#' @param frame_block `function(resp, spec, av, mf, y, n)` run once at
#'   frame assembly, after `y` is coerced and before the random-effect
#'   blocks are built. It returns the block; see the Block section.
#' @param check_frame `function(spec, frame)` run after the predictors
#'   and random-effect blocks exist, for a refusal that depends on the
#'   design rather than on the response. It sees the assembled frame
#'   without the parameter template.
#' @param loglik `function(y, dpars, aterms, weights, block, extra)`
#'   returning the taped log-likelihood of the WHOLE response as one AD
#'   scalar: `sum(log-likelihood)`, not a negative. `y` is the response
#'   after any `block[["y"]]` replacement; `dpars` are the evaluated
#'   distributional parameters on the natural scale, one AD vector of
#'   length `n` or `1` each; `weights` are the effective row weights
#'   with cluster weights folded in, or `1`; `extra` are the family's
#'   extra parameters as an AD list, in the order `extra_pars` declared
#'   them. The family decides what a row weight means for a likelihood
#'   that is not rowwise, and may have refused weights in `check_spec`.
#'   It must not call `RTMB::OBS()`.
#' @param fitted_mean,fitted_var `function(fit, block)` giving the
#'   conditional mean and variance of each row GIVEN the whole observed
#'   response, for [fitted()], `predict(type = "response")` on the
#'   training data, and pearson residuals. `NULL` means "use the rowwise
#'   family mean", which is what a group-level mixture wants. A family
#'   with no mean supplies a function that stops. Both run at the
#'   estimates, outside the tape.
#' @param latent_probs `function(fit, block)` returning one matrix of
#'   posterior latent-state probabilities with column names, `n` rows or
#'   one row per group.
#' @param sim_ctx `function(ctx)` drawing the whole response, the
#'   structured-simulator contract of [frmtmb_family()] with
#'   `ctx[["block"]]` carrying the block. One implementation serves
#'   [simulate()], [posterior_predict()] and [frm_simulate()].
#' @param supports Named logical vector or list of capability flags; see
#'   Capability flags. Unnamed entries and unknown names are refused.
#' @param refusals Named list of one message per `FALSE` flag.
#' @return An object of class `frmtmb_structure`.
#' @seealso [frmtmb_family()], [mixture()]
#' @examples
#' # a structure whose likelihood is the rowwise one, written out: the
#' # smallest thing the protocol accepts
#' st <- frmtmb_structure(
#'   loglik = function(y, dpars, aterms, weights, block, extra) {
#'     sum(weights * RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE))
#'   },
#'   supports = list(conditional_effects = TRUE)
#' )
#' st
#' @srrstats {G2.0,G2.1} Every slot is asserted on type before the
#'   object is built: the function slots must be functions, `keep_na` a
#'   length-one non-missing logical, `supports` a fully named collection
#'   of length-one logicals drawn from a closed vocabulary, and
#'   `refusals` a fully named collection of length-one strings. An
#'   unknown `supports` name errors listing the known ones rather than
#'   being silently carried and never read.
#' @srrstats {G5.2,G5.2a} A refusal a family writes for itself is stored
#'   as data rather than raised from a `stop()` call of its own, so one
#'   capability refused in one place resolves to one message whatever
#'   family is asking.
#' @export
frmtmb_structure <- function(frame_vars = NULL, keep_na = FALSE,
                             check_spec = NULL, frame_block = NULL,
                             check_frame = NULL, loglik,
                             fitted_mean = NULL, fitted_var = NULL,
                             latent_probs = NULL, sim_ctx = NULL,
                             supports = list(), refusals = list()) {
  if (missing(loglik) || !is.function(loglik)) {
    stop("frmtmb_structure(loglik =) is required and must be a function ",
         "(y, dpars, aterms, weights, block, extra) returning the taped ",
         "log-likelihood of the whole response as one AD scalar",
         call. = FALSE)
  }
  check_structure_fn(frame_vars, "frame_vars")
  check_structure_fn(check_spec, "check_spec")
  check_structure_fn(frame_block, "frame_block")
  check_structure_fn(check_frame, "check_frame")
  check_structure_fn(fitted_mean, "fitted_mean")
  check_structure_fn(fitted_var, "fitted_var")
  check_structure_fn(latent_probs, "latent_probs")
  check_structure_fn(sim_ctx, "sim_ctx")
  check_flag(keep_na, "keep_na")
  # The core PASSES weights and lets the family decide what a row weight
  # means for a likelihood that is not rowwise. A family that ignores
  # them silently is the trap that costs a user their weighted fit, so
  # the omission is named at construction rather than discovered later.
  fm <- names(formals(loglik))
  if (!"..." %in% fm && !"weights" %in% fm) {
    warning("frmtmb_structure(loglik =) has no `weights` argument, so the ",
            "row weights the core passes are dropped: a weights() term on ",
            "this family would be accepted and then ignored. Take ",
            "`weights` and use it, or refuse weights() in check_spec",
            call. = FALSE)
  }
  supports <- validate_supports(supports)
  refusals <- validate_refusals(refusals, supports)
  # base::structure(), spelled out: `structure` is also this protocol's
  # own noun and frmtmb_family() takes it as an argument name
  base::structure(
    list(frame_vars = frame_vars, keep_na = isTRUE(keep_na),
         check_spec = check_spec, frame_block = frame_block,
         check_frame = check_frame, loglik = ad_overload_fn(loglik),
         fitted_mean = fitted_mean, fitted_var = fitted_var,
         latent_probs = latent_probs, sim_ctx = sim_ctx,
         supports = supports, refusals = refusals),
    class = "frmtmb_structure"
  )
}

#' Merge declared flags over the conservative defaults.
#'
#' @noRd
validate_supports <- function(supports) {
  out <- frmtmb_structure_flags
  if (is.null(supports) || !length(supports)) return(out)
  if (!is.list(supports) && !is.logical(supports)) {
    stop("frmtmb_structure(supports =) must be a named list or logical ",
         "vector of capability flags", call. = FALSE)
  }
  nms <- names(supports)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("frmtmb_structure(supports =) must name every flag; the known ",
         "flags are ", paste(names(out), collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("frmtmb_structure(supports =) names a flag twice: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "),
         call. = FALSE)
  }
  unknown <- setdiff(nms, names(out))
  if (length(unknown)) {
    stop("frmtmb_structure(supports =): unknown capability flag(s) ",
         paste0("'", unknown, "'", collapse = ", "), ". The flags are ",
         paste(names(out), collapse = ", "), call. = FALSE)
  }
  for (nm in nms) {
    v <- supports[[nm]]
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("frmtmb_structure(supports =): flag '", nm, "' must be TRUE or ",
           "FALSE, not ", arg_desc(v), call. = FALSE)
    }
    out[[nm]] <- v
  }
  out
}

#' A refusal must name a flag that is actually refused, or it is text
#' nobody will ever read.
#'
#' @noRd
validate_refusals <- function(refusals, supports) {
  if (is.null(refusals) || !length(refusals)) return(list())
  if (!is.list(refusals) && !is.character(refusals)) {
    stop("frmtmb_structure(refusals =) must be a named list of one ",
         "message per refused capability", call. = FALSE)
  }
  refusals <- as.list(refusals)
  nms <- names(refusals)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("frmtmb_structure(refusals =) must name the flag each message ",
         "explains", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("frmtmb_structure(refusals =) gives two messages for ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "),
         call. = FALSE)
  }
  for (nm in nms) {
    flag <- refusal_flag(nm)
    if (!flag %in% names(supports)) {
      stop("frmtmb_structure(refusals =): '", nm, "' explains no known ",
           "capability flag. The flags are ",
           paste(names(supports), collapse = ", "), call. = FALSE)
    }
    if (isTRUE(supports[[flag]])) {
      stop("frmtmb_structure(refusals =): '", nm, "' explains a refusal ",
           "of '", flag, "', which supports = declares SUPPORTED, so the ",
           "message could never be shown", call. = FALSE)
    }
    v <- refusals[[nm]]
    if (!is.character(v) || length(v) != 1L || is.na(v) || !nzchar(v)) {
      stop("frmtmb_structure(refusals =): the message for '", nm,
           "' must be one non-empty string", call. = FALSE)
    }
  }
  refusals
}

#' @export
print.frmtmb_structure <- function(x, ...) {
  sup <- x[["supports"]]
  cat("<frmtmb_structure>\n")
  cat("  slots:    ",
      paste(names(Filter(Negate(is.null),
                         x[c("frame_vars", "check_spec", "frame_block",
                             "check_frame", "loglik", "fitted_mean",
                             "fitted_var", "latent_probs", "sim_ctx")])),
            collapse = ", "),
      "\n", sep = "")
  cat("  keeps NA: ", x[["keep_na"]], "\n", sep = "")
  cat("  supports: ",
      if (any(sup)) paste(names(sup)[sup], collapse = ", ") else "(nothing)",
      "\n", sep = "")
  invisible(x)
}

## ---- what the core reads --------------------------------------------

#' The structure of a family, or `NULL`. One spelling, so that every
#' call site is greppable and none of them names a family.
#'
#' @noRd
fam_structure <- function(fam) {
  if (is.null(fam)) NULL else fam[["structure"]]
}

#' The frame block of one response, or `NULL`.
#'
#' @noRd
frame_block_of <- function(frame, resp) {
  (frame[["blocks"]] %||% list())[[resp]]
}

#' Whether a family declares a capability. A family with no structure is
#' rowwise and supports everything the core offers.
#'
#' @noRd
structure_allows <- function(st, flag) {
  is.null(st) || isTRUE(st[["supports"]][[flag]])
}

#' The generic refusal, for a structure that declared a capability
#' unsupported without writing a sentence of its own.
#'
#' @noRd
structure_generic <- function(fam, what) {
  paste0(what, " is not available for a '", fam[["family"]], "' family: ",
         "its likelihood does not factorize over the rows of the data, ",
         "so the quantity this needs per row is not defined")
}

#' The fitting options a structured family may refuse. All three
#' integrate something out with a Laplace approximation about a single
#' inner mode.
#'
#' @noRd
check_structure_fit <- function(spec, REML, quadrature, control) {
  for (resp in spec$responses) {
    fam <- resp$family
    st <- fam_structure(fam)
    if (is.null(st)) next
    if (isTRUE(REML)) {
      structure_gate(st, "reml", structure_generic(fam, "REML = TRUE"))
    }
    if (isTRUE(quadrature)) {
      structure_gate(st, "quadrature",
                     structure_generic(fam, "quadrature = TRUE"))
    }
    if (isTRUE(control$profile)) {
      structure_gate(st, "profile",
                     structure_generic(fam,
                                       "frmtmb_control(profile = TRUE)"))
    }
  }
  invisible(NULL)
}

#' Refuse a capability in the family's own words.
#'
#' `context` names the call site when one flag is refused for two
#' different reasons (`re.form` in `predict()` and in `simulate()`);
#' the bare flag name is the fallback, and `generic` the fallback for a
#' family that declared the refusal without explaining it.
#'
#' @noRd
structure_gate <- function(st, flag, generic, context = NULL) {
  if (structure_allows(st, flag)) return(invisible(NULL))
  ref <- st[["refusals"]]
  msg <- NULL
  if (!is.null(context)) msg <- ref[[paste0(flag, ".", context)]]
  msg <- msg %||% ref[[flag]] %||% generic
  stop(msg, call. = FALSE)
}
