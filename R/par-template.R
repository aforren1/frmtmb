# Discovering the parameter vocabulary. Two questions meet here: what
# names does this model's parameter vector carry, and what values does a
# fit start from. Both are answered by the same object, so a template
# read off one model can be edited and handed straight back as
# `frm(start =)` or `frm_simulate(newparams =)`.

# ---------------------------------------------------------------------
# Component names
# ---------------------------------------------------------------------

#' Display and matching names for one template component. `beta` and
#' `betad` carry design-matrix column names; the covariance components
#' do not, and fall back to the `theta_1` spelling `outer_par_map()`
#' already uses, so `confint(parm =)` rows and a template read the same.
#'
#' @noRd
par_template_names <- function(v, comp) {
  nm <- names(v)
  if (!is.null(nm) && all(nzchar(nm))) return(nm)
  paste0(comp, "_", seq_along(v))
}

#' The template as an editable object: every component a NAMED numeric
#' vector, whatever the frame stored.
#'
#' @noRd
new_par_template <- function(tpl, frame, fitted) {
  for (nm in names(tpl)) {
    tpl[[nm]] <- stats::setNames(as.numeric(tpl[[nm]]),
                                 par_template_names(tpl[[nm]], nm))
  }
  fixed <- list()
  if (length(frame[["betad_fixed_idx"]])) {
    fixed$betad <- as.integer(frame[["betad_fixed_idx"]])
  }
  structure(tpl, fitted = isTRUE(fitted), fixed = fixed,
            class = c("frmtmb_par_template", "list"))
}

# ---------------------------------------------------------------------
# The user-facing verb
# ---------------------------------------------------------------------

#' Parameter names and starting values
#'
#' The parameter vector a model optimizes, as an editable named list:
#' one component per parameter block (`beta`, `betad`, `theta`, ...),
#' each a named numeric vector. It answers "what do I call these?"
#' before there is a fit to ask, which is what `frm(start =)`,
#' [frm_simulate()]'s `newparams =`, and every box constraint all need:
#' these are the names a bound resolves to, and the names
#' [set_prior()]'s `coef` takes for a raw covariance parameter
#' (`class = "theta"`, `coef = "thetaac_1"`).
#'
#' On a fitted model the values are the estimates. On a formula and
#' data the frame is assembled but nothing is fitted, and the values are
#' the ones `frm()` would start from, `start =` and `prior =` included.
#' Either way the result can be edited and passed straight back:
#'
#' ```
#' st <- par_template(bf(y ~ x) + gaussian(), data = dd)
#' st$beta["x"] <- 2
#' frm(bf(y ~ x) + gaussian(), dd, start = st)
#' ```
#'
#' @section The names:
#' Fixed-effect components (`beta`, `betad`) carry the design-matrix
#' column names, so they agree with [fixef()] and `summary()`, and the
#' intercept is spelled `"(Intercept)"`. Covariance components
#' (`theta`, `thetaac`, `thetar`) have no design of their own and are
#' spelled `theta_1`, `theta_2`, ... - the same names `confint(parm =)`
#' reports. When you supply names BACK, the parentheses may be dropped:
#' `"Intercept"` and `"(Intercept)"` address the same coefficient.
#'
#' @section Two vocabularies, one model:
#' This is the flat internal parameterization. [get_prior()] describes
#' the same model in the structured `class`/`coef`/`group`/`dpar`/`nlpar`
#' vocabulary that `set_prior()` addresses, and `frm_simulate()` also
#' accepts a natural-scale spelling (`sigma`, `sd_g__Intercept`). Use
#' `get_prior()` to write priors and `par_template()` to write starting
#' values and box constraints.
#'
#' brms has no counterpart for starting values: `brm(init =)` takes a
#' list of Stan program names and its own documentation calls it mainly
#' an internal testing facility. brms does not need one, because it
#' initializes at random in a bounded range and its priors carry
#' location information. `frm()` optimizes rather than samples, so it
#' evaluates the objective AT the starting values; this is the discovery
#' step that makes that survivable.
#'
#' @param object A `frmtmb_fit`, or a [bf()] formula or plain formula.
#' @param data Model data. Required when `object` is a formula.
#' @param family Family, when `object` does not carry one.
#' @param start Optional `start` list, applied exactly as `frm()` would,
#'   so the result shows the starting values a fit with that `start`
#'   uses. Formula method only.
#' @param prior Optional [set_prior()] specification, so the result
#'   shows the nonlinear starting values its locations place. Formula
#'   method only.
#' @param na.action,data2 As in [frm()]. Formula method only.
#' @param ... Unused.
#' @return A `frmtmb_par_template`: a named list of named numeric
#'   vectors, accepted directly as `frm(start =)` and as
#'   `frm_simulate(newparams =)`.
#' @seealso [frm()] for `start =`, [get_prior()] for the prior
#'   addressing vocabulary, [frm_simulate()] for `newparams =`.
#' @examples
#' dd <- data.frame(x = rnorm(30), g = factor(rep(1:5, 6)))
#' dd$y <- 1 + 2 * dd$x + rnorm(30)
#'
#' # before fitting: the names and the cold starting values
#' par_template(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
#'
#' # after fitting: the same layout, holding the estimates
#' fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd)
#' par_template(fit)
#'
#' # discover, edit, fit
#' st <- par_template(bf(y ~ x) + gaussian(), data = dd)
#' st$beta["x"] <- 2
#' frm(bf(y ~ x) + gaussian(), dd, start = st)
#' @export
par_template <- function(object, ...) UseMethod("par_template")

#' @rdname par_template
#' @export
par_template.frmtmb_fit <- function(object, ...) {
  new_par_template(object$estimates, object$frame,
                   fitted = !inherits(object, "frmtmb_unfitted"))
}

#' @rdname par_template
#' @export
par_template.default <- function(object, data, family = NULL,
                                 start = NULL, prior = NULL,
                                 na.action = stats::na.omit,
                                 data2 = list(), ...) {
  if (missing(data)) {
    stop("par_template() needs `data` to assemble the design: the ",
         "parameter vector's length and names are properties of the ",
         "model matrices, not of the formula alone", call. = FALSE)
  }
  if (!is.null(start)) {
    check_named_list(start, "start", "start = list(beta = c(0, 1))")
  }
  prior <- as_priorlist(prior)
  bform <- resolve_deferred_families(as_bform(object, family), data)
  spec <- parse_spec(bform)
  frame <- assemble_frame(spec, data, na.action = na.action,
                          data2 = validate_data2(data2))
  entries <- if (!is.null(prior)) {
    resolve_prior_input(list(frame = frame, spec = spec), prior)$entries
  }
  new_par_template(make_start(frame, start, entries), frame, fitted = FALSE)
}

#' @export
print.frmtmb_par_template <- function(x, n = 10L, ...) {
  cat("<frmtmb parameter template> ",
      if (attr(x, "fitted")) "estimates" else "starting values", "\n",
      sep = "")
  fixed <- attr(x, "fixed") %||% list()
  for (nm in names(x)) {
    v <- x[[nm]]
    cat("$", nm, "  (", length(v), ")\n", sep = "")
    show <- if (length(v) > n) v[seq_len(n)] else v
    print(show)
    if (length(v) > n) {
      cat("... ", length(v) - n, " more\n", sep = "")
    }
    if (!is.null(fixed[[nm]])) {
      cat("(held at a constant, not optimized: ",
          paste(names(v)[fixed[[nm]]], collapse = ", "), ")\n", sep = "")
    }
  }
  cat("Edit and pass back as frm(start = ) or ",
      "frm_simulate(newparams = ). Parentheses are optional in names ",
      "you supply.\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------
# Writing start values into a component
# ---------------------------------------------------------------------

#' One `start` component resolved against its template slot.
#'
#' An UNNAMED vector is positional and must be full length, which is the
#' contract every existing call relies on. A NAMED vector addresses
#' entries by name through the same paren-tolerant matcher `confint()`
#' uses, overrides those and leaves the rest at their
#' defaults. Mixing the two in one vector is refused rather than guessed
#' at: neither reading is safe when half the entries are labeled.
#'
#' @noRd
resolve_start_component <- function(cur, val, comp) {
  if (!is.numeric(val)) {
    stop("start$", comp, " must be a numeric vector, not ", arg_desc(val),
         call. = FALSE)
  }
  nms <- names(val)
  if (is.null(nms)) {
    if (length(val) != length(cur)) {
      stop("start$", comp, " must have length ", length(cur),
           call. = FALSE)
    }
    cur[] <- val
    return(cur)
  }
  if (!all(nzchar(nms))) {
    stop("start$", comp, " mixes named and unnamed entries; name every ",
         "entry or none. par_template() lists the names", call. = FALSE)
  }
  target <- par_template_names(cur, comp)
  # match_par_name() raises the ambiguity error itself, in the one
  # wording every parm-style argument already uses
  idx <- match_par_name(nms, target)
  if (anyNA(idx)) {
    # a full-length vector whose names resolve to NOTHING is usually a
    # positional start carrying names from somewhere else (coef() of a
    # different model, a covariance matrix's dimnames), which the
    # positional-only contract of earlier releases ignored
    stop("start$", comp, " names no parameter of this model: ",
         paste(nms[is.na(idx)], collapse = ", "),
         ". It has ", paste(target, collapse = ", "),
         " (see par_template())",
         if (all(is.na(idx)) && length(val) == length(cur)) {
           ". To set it positionally instead, unname() the vector"
         }, call. = FALSE)
  }
  if (anyDuplicated(idx)) {
    stop("start$", comp, " addresses one parameter more than once: ",
         paste(unique(target[idx[duplicated(idx)]]), collapse = ", "),
         call. = FALSE)
  }
  cur[idx] <- as.numeric(val)
  cur
}

#' Positions of one component that a `start` entry sets. Called only
#' after resolve_start_component() accepted the same vector, so the name
#' match cannot fail here.
#'
#' @noRd
start_claimed_idx <- function(cur, val, comp) {
  if (is.null(val)) return(integer(0))
  nms <- names(val)
  if (is.null(nms)) return(seq_along(cur))
  as.integer(match_par_name(nms, par_template_names(cur, comp)))
}

# ---------------------------------------------------------------------
# Prior-placed nonlinear starts
# ---------------------------------------------------------------------

#' Positions in `beta` that belong to a nonlinear parameter's own
#' sub-formula. The nl body itself has no design matrix, so it
#' contributes none.
#'
#' @noRd
nl_beta_positions <- function(frame) {
  out <- integer(0)
  for (lp in frame[["linpreds"]]) {
    if (!identical(lp[["par"]], "beta")) next
    if (!is.null(lp[["nl_body"]]) || !is.null(lp[["constant"]])) next
    nl <- frame[["spec"]]$responses[[lp[["resp"]]]]$nlpars %||% character(0)
    if (lp[["dpar"]] %in% nl) out <- c(out, as.integer(lp[["idx"]]))
  }
  out
}

#' Starting values a prior's locations place on NONLINEAR parameters.
#'
#' A nonlinear body is usually flat, singular or undefined at zero, and
#' zero is where make_start() otherwise leaves it: a family's
#' `init_dpars` can only seed an intercept of its own linear predictor,
#' and an nl mu has no design. That is the documented brms#734 failure.
#' A prior that says `normal(5000, 1000)` on `ult` has already named the
#' region, and brms uses exactly that information to initialize, so read
#' it here too.
#'
#' Scoped to nonlinear coefficients on purpose. Every other parameter
#' keeps the start it has always had: a prior is a penalty, not a claim
#' about where the optimizer should begin, and moving every start would
#' change optimizer paths throughout for no failure it fixes.
#'
#' @noRd
prior_nl_starts <- function(frame, entries) {
  pos <- nl_beta_positions(frame)
  out <- list()
  if (!length(pos) || !length(entries)) return(out)
  nms <- par_template_names(frame[["par_template"]][["beta"]], "beta")
  for (e in entries) {
    if (!identical(e$comp, "beta") || length(e$idx) != 1L) next
    if (!identical(e$scale %||% "internal", "internal")) next
    if (!e$idx %in% pos) next
    loc <- prior_dist_location(e$dist)
    if (is.na(loc)) next
    out[[nms[e$idx]]] <- list(idx = as.integer(e$idx), value = loc)
  }
  out
}
