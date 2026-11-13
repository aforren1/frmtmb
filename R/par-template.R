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
  paste0(comp, "_", seq_along(v), recycle0 = TRUE)
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
#' [set_prior()]'s `coef` takes for a covariance parameter that no brms
#' class names (`class = "theta"`, `coef = "theta_12"`).
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
#' @section The default starting values:
#' For most location-scale families the coefficients of `mu` start at a
#' least-squares fit on its fixed design, and the scale (`sigma`,
#' `shape` or `phi`) starts at the spread of that fit's residuals. Thus
#' a covariate that explains most of the response does not start the
#' scale at the spread of the response. The families are gaussian(),
#' student(), skew_normal(), exgaussian(), asym_laplace(),
#' zero_inflated_asym_laplace(), huber(), lognormal(),
#' shifted_lognormal() and hurdle_lognormal() with the identity link of
#' `mu`, Gamma() and weibull() and negbinomial() with the log link, and
#' Beta() with the logit link. Every other fit, and every fit with a
#' nonlinear or intercept-only `mu`, a smooth or `gp()` term in `mu`, an
#' ill-conditioned design, a censored, truncated or `se()` response, or a
#' scale held at a constant, starts each intercept at the family's
#' initializer. A nonlinear coefficient starts at 0, or at the location
#' of a prior on it. Where the body is defined at those values and the
#' link cannot evaluate it, as under the default links of Gamma() and
#' inverse.gaussian(), one intercept the body reads moves to where the
#' body equals the response mean, and a message says so. The entries you
#' name in `start` replace these values. The random-effect standard
#' deviations start at the response's unit.
#'
#' A body can be constant across the rows at these values and stationary
#' there, as `b1 * exp(b2 * x)` is at 0. [frm()] then restarts from moved
#' coefficients of every such body and keeps a higher point, and it does the
#' same for a fit from the default start that runs along a flat direction. A
#' message says so. If no restart leaves a constant body, a warning names it:
#' the data may simply have no effect there, and a `start` can test that. A
#' message names a constant body that is a local maximum along its coefficients.
#'
#' A body can also have more than one true maximum, and a default start finds
#' one of them. Nothing in that fit marks it. So from the default start
#' [frm()] also fits from each sign pattern of the coefficients that the body
#' does not read as an added term, and returns the highest point. With four
#' such coefficients or fewer that is every pattern; with more, every
#' pattern of the four the body moves most by, so the order of the terms
#' does not change the answer. A coefficient inside a curve starts where it
#' moves the curve's argument by about 3 (`b2` near 1e-3 in
#' `b1 * exp(b2 * x)` with `x` in the thousands), and one the body is linear
#' in where it moves the body by about two of the response's units, each
#' only where that is far from 1, else at 1. A message (class
#' `frmtmb_nl_multistart`) says so when that point is higher than the first
#' fit by more than 0.01. `b1 * exp(b2 * x)` under [bernoulli()] on one data
#' set reads -193.4967 from the default start alone and -189.999 with the
#' patterns. This makes a fit 1.2 to 2.8 times as slow with two or three such
#' coefficients, 3.5 to 5.4 times with four, and up to 47 times with five. A
#' `start`, or `options(frmtmb.nl_sign_probe = FALSE)`, skips it; a located
#' prior does not. The patterns are not a proof: a `start` elsewhere, or a
#' located prior, can still show a higher maximum. [frm_allfit()] cannot, as
#' its optimizers start from one point.
#'
#' A nonlinear coefficient that the body reads as an added term (`a` and `b`
#' of `a + b * x`) is stepped in the units by which it moves the body, where
#' these are 16 or more times from 1: `a * x` with `x` near `1e-20` reaches
#' the maximum of `y ~ 0 + x`. The curvature of a normal, Student t or
#' logistic prior on it is added; under another prior with a curved density,
#' such as a gamma or lognormal prior, it keeps its old steps. In the sign
#' patterns' runs and the refit of their winner, a coefficient inside a curve
#' is stepped in its start size where that is 16 or more times from 1 and no
#' prior reaches it, so `b2` of `b1 * (1 - exp(-b2 * x))` fits alike with `x`
#' near 1 and near 1e-6.
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
#' @param ... Refused: an argument the method does not have is an
#'   error naming it, rather than silently changing nothing.
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
  frm_check_dots(...)
  new_par_template(object$estimates, object$frame,
                   fitted = !inherits(object, "frmtmb_unfitted"))
}

#' @rdname par_template
#' @export
par_template.default <- function(object, data, family = NULL,
                                 start = NULL, prior = NULL,
                                 na.action = stats::na.omit,
                                 data2 = list(), ...) {
  frm_check_dots(...)
  if (missing(data)) {
    frm_stop("par_template() needs `data` to assemble the design: the ",
             "parameter vector's length and names are properties of the ",
             "model matrices, not of the formula alone", call. = FALSE)
  }
  if (!is.null(start)) {
    check_named_list(start, "start", "start = list(beta = c(0, 1))")
  }
  prior <- as_priorlist(prior)
  check_prior_slots(prior)
  bform <- expand_dot_bform(as_bform(object, family), data)
  bform <- resolve_deferred_families(bform, data)
  spec <- parse_spec(bform)
  frame <- assemble_frame(spec, data, na.action = na.action,
                          data2 = validate_data2(data2))
  # resolve_prior_input() reads the spec's families, so it needs the
  # finalized ones: a family that derived its dpars from the response
  # would otherwise be priced against the vocabulary it was written
  # with. No in-repo family does that today.
  spec <- carry_finalized_responses(spec, frame)
  entries <- if (!is.null(prior)) {
    resolve_prior_input(list(frame = frame, spec = spec), prior)$entries
  }
  new_par_template(make_start(frame, start, entries), frame, fitted = FALSE)
}

#' @export
print.frmtmb_par_template <- function(x, n = 10L, ...) {
  frm_check_dots(...)
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
    frm_stop("start$", comp, " must be a numeric vector, not ", arg_desc(val),
             call. = FALSE)
  }
  nms <- names(val)
  if (is.null(nms)) {
    if (length(val) != length(cur)) {
      frm_stop("start$", comp, " must have length ", length(cur),
               call. = FALSE)
    }
    cur[] <- val
    return(cur)
  }
  if (!all(nzchar(nms))) {
    frm_stop("start$", comp, " mixes named and unnamed entries; name every ",
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
    frm_stop("start$", comp, " names no parameter of this model: ",
             paste(nms[is.na(idx)], collapse = ", "),
             ". It has ", paste(target, collapse = ", "),
             " (see par_template())",
             if (all(is.na(idx)) && length(val) == length(cur)) {
               ". To set it positionally instead, unname() the vector"
             }, call. = FALSE)
  }
  if (anyDuplicated(idx)) {
    frm_stop("start$", comp, " addresses one parameter more than once: ",
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

#' The value of a body at a template's coefficients, with every
#' random effect at its start of 0: each parameter it reads through its
#' design, then the bodies in dependency order. `NULL` for a body that
#' reads a `ps()` block, which has no value before its closure exists.
#'
#' @noRd
nl_body_at <- function(frame, lp, tpl) {
  vals <- list()
  for (q in frame[["linpreds"]]) {
    if (!identical(q[["resp"]], lp[["resp"]])) next
    if (!is.null(q[["nl_body"]])) {
      if (length(q[["ps_terms"]] %||% list())) return(NULL)
      eta <- eval(q[["nl_body"]],
                  c(vals[c(q[["nl_pars"]], q[["nl_dpar_refs"]])],
                    q[["data_list"]]),
                  ad_overload_env(q[["nl_env"]], q[["nl_body"]]))
    } else {
      b <- tpl[[q[["par"]]]][q[["idx"]]]
      eta <- if (length(b)) {
        as.numeric(as.matrix(q[["X"]] %*% b))
      } else {
        numeric(nrow(q[["X"]]))
      }
      if (!is.null(q[["offset"]])) eta <- eta + q[["offset"]]
    }
    if (identical(q[["dpar"]], lp[["dpar"]])) return(as.numeric(eta))
    vals[[q[["dpar"]]]] <- q[["link"]]$linkinv(eta)
  }
  NULL
}

#' The template positions of the coefficients a body reads, directly or
#' through further bodies: `list(<par> = positions)`, with `seen` the
#' parameters walked and `hollow` one whose formula has terms and whose
#' design has no column.
#'
#' @noRd
nl_body_coefs <- function(frame, lp) {
  lps <- frame[["linpreds"]]
  rn <- lp[["resp"]]
  seen <- character(0)
  front <- lp[["nl_pars"]] %||% character(0)
  coefs <- list()
  hollow <- character(0)
  while (length(front)) {
    np <- front[1L]
    front <- front[-1L]
    if (np %in% seen) next
    seen <- c(seen, np)
    q <- lps[[linpred_key(rn, np)]]
    if (is.null(q)) next
    if (!is.null(q[["nl_body"]])) {
      front <- c(front, q[["nl_pars"]] %||% character(0))
      next
    }
    if (length(q[["idx"]])) {
      coefs[[q[["par"]]]] <- c(coefs[[q[["par"]]]], q[["idx"]])
      next
    }
    tt <- q[["terms"]]
    has_terms <- !is.null(tt) && (attr(tt, "intercept") == 1L ||
                                    length(attr(tt, "term.labels")))
    if (has_terms && is.null(q[["Z"]]) && is.null(q[["offset"]]) &&
          is.null(q[["constant"]]) && !length(q[["smooths"]])) {
      hollow <- c(hollow, np)
    }
  }
  list(coefs = coefs, seen = seen, hollow = hollow)
}

#' A template with a body's coefficients moved, each by its own irrational
#' multiple of `step`, so that no linear body (`2 * a - b`, `a - 2 * b + c`)
#' has the move as a null direction.
#'
#' @noRd
nl_moved_tpl <- function(tpl, coefs, step) {
  k <- 0L
  for (par in names(coefs)) {
    at <- coefs[[par]]
    tpl[[par]][at] <- tpl[[par]][at] +
      step * sqrt(2 + (k + seq_along(at)) * 0.6180339887498949)
    k <- k + length(at)
  }
  tpl
}

#' The taped gradient at the start, by template component, where a body's
#' coefficients can be read from it: `NULL` with random effects (a gradient
#' there runs the inner optimizer), a whitened or profiled `beta`, or a frame
#' with no body.
#'
#' @noRd
nl_guard_gradient <- function(obj, frame, bw) {
  has_body <- any(vapply(frame[["linpreds"]], function(lp) {
    !is.null(lp[["nl_body"]])
  }, NA))
  if (!has_body || !is.null(bw) || length(obj$env$random)) return(NULL)
  g <- tryCatch(suppressWarnings(as.numeric(obj$gr(obj$par))),
                error = function(e) NULL)
  if (is.null(g) || length(g) != length(obj$par)) return(NULL)
  out <- list()
  for (par in unique(names(obj$par))) {
    tp <- frame[["par_template"]][[par]]
    gp <- g[names(obj$par) == par]
    # a mapped component has fewer entries than the template; its
    # positions cannot be read off, and the guard then asks the objective
    if (length(tp) == length(gp)) out[[par]] <- gp
  }
  out
}

#' Refuse a nonlinear body the model does not read.
#'
#' Through 0.78.0 an ordinal family's rule that thresholds replace the
#' intercept also took the intercept of `b1 ~ 1`, so the body was 0 on
#' every row and frm() fitted the thresholds alone without a word. The
#' cause was in the design and is fixed there; this guards the class: a
#' parameter whose formula has terms and whose design came out with no
#' column, or a family's density that does not read the dpar the body
#' computes.
#'
#' The second is one question with two halves, and the refusal needs both:
#' the body's VALUE moved, by more than the response's rounding could hide,
#' and the objective did not. A body that did not move says nothing about
#' the family: `step(x - c)` flat around a start, or `exp(-k * x)` saturated,
#' is a start problem, and `a * x` with x about 1e-17 is a scale problem.
#' A taped gradient that is not zero answers without any of this, and the
#' objective on doubles runs only for what is left. A non-finite value, or
#' an objective that cannot run on doubles, passes.
#'
#' @noRd
nl_check_body_read <- function(frame, nll, tpl, grad = NULL) {
  lps <- frame[["linpreds"]]
  for (lp in lps) {
    if (is.null(lp[["nl_body"]])) next
    rn <- lp[["resp"]]
    bc <- nl_body_coefs(frame, lp)
    if (length(bc$hollow)) {
      frm_stop("The nonlinear parameter '", bc$hollow[1L], "' is read by the ",
               "body of '", lp[["dpar"]], "', but its formula left it with no ",
               "coefficient, so the body would be fitted without it. This ",
               "is a fault in frmtmb, not in the model; please report it",
               call. = FALSE)
    }
    coefs <- bc$coefs
    if (!length(coefs)) next
    if (!is.null(grad) &&
          any(unlist(lapply(names(coefs), function(par) {
            g <- grad[[par]]
            !is.null(g) && any(is.finite(g[coefs[[par]]]) &
                                 g[coefs[[par]]] != 0)
          })))) {
      next
    }
    y <- frame[["y"]][[rn]]
    sy <- if (is.numeric(y)) stats::sd(as.numeric(y), na.rm = TRUE) else 1
    if (!is.finite(sy) || sy <= 0) sy <- 1
    moved_at <- function(step) {
      t1 <- nl_moved_tpl(tpl, coefs, step)
      e0 <- tryCatch(suppressWarnings(nl_body_at(frame, lp, tpl)),
                     error = function(e) NULL)
      e1 <- tryCatch(suppressWarnings(nl_body_at(frame, lp, t1)),
                     error = function(e) NULL)
      if (is.null(e0) || is.null(e1)) return(NULL)
      m0 <- suppressWarnings(lp[["link"]]$linkinv(e0))
      m1 <- suppressWarnings(lp[["link"]]$linkinv(e1))
      d <- abs(m1 - m0)
      d <- d[is.finite(d)]
      if (length(d) && max(d) > 1e-4 * sy) t1
    }
    v0 <- NULL
    unmoved <- 0L
    for (step in c(0.1, -0.15)) {
      t1 <- moved_at(step)
      if (is.null(t1)) next
      v0 <- v0 %||% tryCatch(suppressWarnings(as.numeric(nll(tpl))),
                             error = function(e) NA_real_)
      v1 <- tryCatch(suppressWarnings(as.numeric(nll(t1))),
                     error = function(e) NA_real_)
      if (!is.finite(v0) || !is.finite(v1) || v0 != v1) {
        unmoved <- -1L
        break
      }
      unmoved <- unmoved + 1L
    }
    if (unmoved > 0L) {
      frm_stop("The nonlinear body of '", lp[["dpar"]], "' does not reach ",
               "the likelihood of family '",
               frame[["spec"]]$responses[[rn]]$family[["family"]],
               "': moving the coefficients of ",
               paste0("'", bc$seen, "'", collapse = ", "),
               " moves the body and leaves the objective unchanged. The fit ",
               "would ignore the body, so it is refused. A custom family's ",
               "lpdf has to read dpars[[\"", lp[["dpar"]], "\"]]",
               call. = FALSE)
    }
  }
  invisible(NULL)
}

#' The bodies sitting on a stationary point where they are constant across
#' rows.
#'
#' `b1 * exp(b2 * x)` at b1 = b2 = 0 is constant, and on an ordinal family
#' the thresholds start at their exact maximum for a constant predictor, so
#' every gradient of a body coefficient there is 0 and the optimizer stops at
#' once, at the log-likelihood of the model without the body. The same point
#' is reached under any family from a body whose slopes multiply,
#' `a + b1 * b2 * x`. It is usually a saddle: not a maximum, and not a sign
#' that the data do not determine the coefficients.
#'
#' Every body is looked at, because a model has more than one and they sit
#' there independently: moving only the first left `sigma` or `mu2` on the
#' point, 16.8 and 65.7 nats short.
#'
#' Returns `saddle` and `local`, each a list of bodies. `local` holds the
#' points that are a true local maximum along the body's coefficients, which
#' the curvature says (`a + b^2 * x` with a negative slope in the data). A
#' body that varies across the rows costs one evaluation and is in neither,
#' nor is a body that cannot vary at all (`a`, or `a + b * x` with x constant).
#'
#' @noRd
nl_saddle_find <- function(fit, use_curvature = TRUE) {
  frame <- fit$frame
  est <- fit$estimates
  saddle <- list()
  local <- list()
  for (lp in frame[["linpreds"]]) {
    if (is.null(lp[["nl_body"]])) next
    bc <- nl_body_coefs(frame, lp)
    if (!length(bc$coefs)) next
    e0 <- tryCatch(suppressWarnings(nl_body_at(frame, lp, est)),
                   error = function(e) NULL)
    if (is.null(e0) || length(e0) < 2L || !all(is.finite(e0))) next
    if (diff(range(e0)) > 1e-8 * max(1, abs(mean(e0)))) next
    e1 <- tryCatch(
      suppressWarnings(nl_body_at(frame, lp, nl_moved_tpl(est, bc$coefs, 0.1))),
      error = function(e) NULL)
    if (is.null(e1) || !any(is.finite(e1))) next
    if (diff(range(e1, finite = TRUE)) <= 1e-6 * max(1, abs(mean(e1)))) next
    one <- list(lp = lp, coefs = bc$coefs, seen = bc$seen)
    if (use_curvature && nl_coef_curved(frame, est, bc$coefs)) {
      local[[length(local) + 1L]] <- one
    } else {
      saddle[[length(saddle) + 1L]] <- one
    }
  }
  list(saddle = saddle, local = local)
}

#' Is the negative log-likelihood strictly convex along a body's
#' coefficients, the rest of the model held? Second differences on doubles,
#' at most eight coefficients; a value that will not evaluate answers no.
#'
#' @noRd
nl_coef_curved <- function(frame, est, coefs) {
  pos <- do.call(rbind, lapply(names(coefs), function(par) {
    data.frame(par = par, at = coefs[[par]], stringsAsFactors = FALSE)
  }))
  k <- nrow(pos)
  if (k > 8L) return(FALSE)
  nll <- tryCatch(build_objective(frame), error = function(e) NULL)
  if (is.null(nll)) return(FALSE)
  h <- 1e-3 * vapply(seq_len(k), function(i) {
    max(1, abs(est[[pos$par[i]]][pos$at[i]]))
  }, 0)
  val <- function(s) {
    t1 <- est
    for (i in which(s != 0)) {
      t1[[pos$par[i]]][pos$at[i]] <- t1[[pos$par[i]]][pos$at[i]] + s[i] * h[i]
    }
    tryCatch(suppressWarnings(as.numeric(nll(t1))),
             error = function(e) NA_real_)
  }
  f0 <- val(rep(0, k))
  H <- matrix(0, k, k)
  for (i in seq_len(k)) {
    ei <- replace(rep(0, k), i, 1)
    H[i, i] <- (val(ei) + val(-ei) - 2 * f0) / h[i]^2
    for (j in seq_len(i - 1L)) {
      ej <- replace(rep(0, k), j, 1)
      H[i, j] <- H[j, i] <- (val(ei + ej) - val(ei - ej) - val(ej - ei) +
                               val(-ei - ej)) / (4 * h[i] * h[j])
    }
  }
  if (!all(is.finite(H))) return(FALSE)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > 1e-6 * max(abs(ev))
}

#' The deterministic restarts from a saddle: every coefficient of every named
#' body moved off the point in one start, a slope by less where its column is
#' wide, and the sign of one parameter of body `flip` at a time flipped, at most
#' four starts in all. They read neither the session's RNG nor the response.
#'
#' A body before `flip` starts where `kept` has it. Flips are per body because
#' the bodies need their signs independently: `mvbf()` with two `exp()` bodies
#' needs a negative slope on the second, and four flips shared over the two
#' never reached it, 3.6 nats short.
#'
#' @noRd
nl_saddle_starts <- function(frame, bodies, start, flip = 1L, kept = NULL) {
  frame_names <- function(par) {
    par_template_names(frame[["par_template"]][[par]], par)
  }
  lps <- frame[["linpreds"]]
  pars <- list()
  owner <- integer(0)
  taken <- character(0)
  for (k in seq_along(bodies)) {
    rn <- bodies[[k]]$lp[["resp"]]
    for (p in bodies[[k]]$seen) {
      key <- linpred_key(rn, p)
      q <- lps[[key]]
      if (key %in% taken || is.null(q) || !is.null(q[["nl_body"]]) ||
            !length(q[["idx"]])) {
        next
      }
      taken <- c(taken, key)
      pars[[length(pars) + 1L]] <- q
      owner <- c(owner, k)
    }
  }
  signs <- list(rep(1, length(pars)))
  for (i in which(owner == flip)) {
    if (length(signs) >= 4L) break
    s <- rep(1, length(pars))
    s[i] <- -1
    signs[[length(signs) + 1L]] <- s
  }
  lapply(signs, function(s) {
    st <- start %||% list()
    for (i in seq_along(pars)) {
      q <- pars[[i]]
      par <- q[["par"]]
      nm <- frame_names(par)
      X <- as.matrix(q[["X"]])
      w <- apply(X, 2L, function(v) {
        r <- diff(range(v))
        if (is.finite(r) && r > 0) max(1, r / 2) else 1
      })
      v <- if (!is.null(kept) && owner[i] < flip) {
        kept$estimates[[par]][q[["idx"]]]
      } else {
        s[i] * 0.1 * sqrt(2 + i * 0.6180339887498949) /
          w[seq_along(q[["idx"]])]
      }
      names(v) <- nm[q[["idx"]]]
      old <- st[[par]]
      if (!is.null(old) && is.null(names(old))) names(old) <- nm
      st[[par]] <- c(old[setdiff(names(old), names(v))], v)
    }
    st
  })
}

#' "the body 'x'" or "the bodies 'x' and 'y'", with the verb's number.
#'
#' @noRd
nl_bodies_words <- function(bodies) {
  txt <- vapply(bodies, function(b) deparse1(b$lp[["nl_body"]]), "")
  txt <- paste0("'", txt, "'")
  if (length(txt) == 1L) {
    return(list(noun = paste0("the body ", txt), is = "is", its = "its"))
  }
  list(noun = paste0("the bodies ", paste(txt[-length(txt)], collapse = ", "),
                     " and ", txt[length(txt)]),
       is = "are", its = "their")
}

#' What the model is at these constant bodies. "Thresholds-only" is true only
#' where nothing on the response still varies: a `nlf(disc ~ ...)` constant
#' beside a `mu` body that moves is not that model.
#'
#' @noRd
nl_model_words <- function(frame, bodies) {
  plain <- "the model without those covariates"
  rn <- unique(vapply(bodies, function(b) b$lp[["resp"]], ""))
  if (length(rn) != 1L) return(plain)
  rs <- frame[["spec"]]$responses[[rn]]
  if (!isTRUE(rs$family[["drop_intercept"]])) return(plain)
  named <- vapply(bodies, function(b) b$lp[["dpar"]] %||% "", "")
  if (!"mu" %in% named) return(plain)
  under <- unlist(lapply(bodies, function(b) b$seen))
  for (lp in frame[["linpreds"]]) {
    if (!identical(lp[["resp"]], rn)) next
    dp <- lp[["dpar"]] %||% ""
    if (dp %in% c(named, under)) next
    if (!is.null(lp[["nl_body"]]) || NCOL(lp[["X"]]) > 1L ||
          !is.null(lp[["Z"]]) || length(lp[["smooths"]])) {
      return(plain)
    }
  }
  "the thresholds-only model"
}

#' Two log-likelihoods, printed with the digits that tell them apart.
#'
#' @noRd
nl_ll_pair <- function(a, b) {
  d <- 7L
  while (d < 15L && identical(format(a, digits = d), format(b, digits = d))) {
    d <- d + 1L
  }
  c(format(a, digits = d), format(b, digits = d))
}

#' Do these held warnings report a flat direction in the standard errors?
#' A family parameter at the end of its range is one too (its message,
#' class `frmtmb_limit_fit`, is held beside the warnings): a student()
#' `nu` that the default start ran to exp(22.7) on a `pow_a` body with
#' gp() stopped at a log-likelihood of -4.95, and the restart reaches
#' 55.67 (lane gpstop's sweep, `pow_a_student_n120_s2`).
#'
#' @noRd
nl_flat_warning <- function(held) {
  any(vapply(held, function(w) {
    inherits(w, "frmtmb_limit_fit") ||
      grepl("Standard errors are not available|are not identified",
            conditionMessage(w))
  }, NA))
}

#' Do these held warnings say the stop may not be the maximum: the optimizer's
#' code, a gradient left at the stop, or a gradient that is not finite there
#' (check_convergence())? By class, so a reworded message cannot turn the
#' ridge pass off. "Hessian is not positive definite" is a standard-error
#' verdict and is left out.
#'
#' @noRd
nl_stop_warning <- function(held) {
  any(vapply(held, function(w) {
    inherits(w, c("frmtmb_conv_code", "frmtmb_conv_gradient",
                  "frmtmb_conv_gradient_nonfinite"))
  }, NA))
}

#' The restarts of one round: the highest fit above `floor`, with its held
#' warnings, or `NULL` in `best`.
#'
#' Two restarts that agree have found the point this start pattern finds, and
#' the third of a set can cost more than the first two together (12 s against
#' 5 s at 3e4 rows). They agree only if their estimates do as well as their
#' values. Equal values alone can be two trips up one ridge: `b1 * (1 - exp(-b2 *
#' x))` under `cratio()` reads -377.0010 from the first two starts, at b1 = 21293
#' and 22565, and -373.2498 only from the third, which gives b2 the sign it needs.
#' A standard error need not report that ridge. Agreeing on the saddle's own value
#' says nothing either (`saddle_floor`). A restart from a flat direction that
#' agrees with the first fit has its answer at once (`from_flat`).
#'
#' @noRd
nl_restart_round <- function(run, starts, floor, saddle_floor, from_flat) {
  same <- function(u, v) abs(u - v) <= 1e-7 * max(1, abs(u))
  best <- NULL
  best_ll <- floor
  best_held <- list()
  for (st in starts) {
    h1 <- new.env(parent = emptyenv())
    h1$w <- list()
    r <- tryCatch(
      withCallingHandlers(suppressMessages(run(st, FALSE)),
                          warning = function(w) {
                            h1$w[[length(h1$w) + 1L]] <- w
                            invokeRestart("muffleWarning")
                          }),
      error = function(e) NULL)
    if (is.null(r)) next
    # the restarts are judged on what they maximize, which under a prior
    # holds the log prior
    ll <- suppressWarnings(as.numeric(logLik_penalized(r)))
    if (!is.finite(ll)) next
    if (ll > best_ll + 1e-8 * max(1, abs(best_ll))) {
      best <- r
      best_ll <- ll
      best_held <- h1$w
    }
    if (from_flat && same(ll, floor)) break
  }
  list(best = best, ll = best_ll, held = best_held)
}

#' `frm()`'s fit of a model with a nonlinear body, kept off a saddle.
#'
#' The warnings of a fit are held until it is the one returned, because a
#' saddle's warning says "the data do not determine them", which is false,
#' and a restart replaces that fit. A fit that converged with its standard
#' errors and is not on a saddle comes back as `run()` returned it, its
#' warnings after it, and restarts cost nothing there. Any fit from a default
#' start that a restart does not beat comes back the same way, digit for
#' digit.
#'
#' A kept restart is looked at again, because a body it did not move can sit on
#' the point in its turn: the next pass moves that body as well. What is left
#' is a warning and never a stop. The restarts have measured the point, and a
#' body with no effect in the data rests there correctly: an exact zero slope on
#' `a + b1 * b2 * x`, or 1 of 400 null data sets at n = 40, where a stop would
#' end a simulation loop.
#'
#' @noRd
nl_fit_off_saddle <- function(run, frame, start, prior) {
  hq <- new.env(parent = emptyenv())
  hq$w <- list()
  hold <- function(w) {
    hq$w[[length(hq$w) + 1L]] <- w
    invokeRestart("muffleWarning")
  }
  # a family limit's message is held as the warnings are: it marks the
  # flat direction a restart is for, and it is about this fit only
  hold_limit <- function(m) {
    hq$w[[length(hq$w) + 1L]] <- m
    invokeRestart("muffleMessage")
  }
  release <- function() {
    for (w in hq$w) if (inherits(w, "message")) message(w) else warning(w)
  }
  fit <- tryCatch(withCallingHandlers(run(start, TRUE), warning = hold,
                                      frmtmb_limit_fit = hold_limit),
                  error = function(e) {
                    release()
                    stop(e)
                  })
  held <- hq$w
  # under a prior the curvature on the likelihood alone says nothing about
  # the mode, so a restart decides
  use_curv <- is.null(prior) || NROW(prior) == 0L
  look <- nl_saddle_find(fit, use_curvature = use_curv)
  on_saddle <- length(look$saddle) > 0L
  # A default start the optimizer ran along a ridge is the other case worth a
  # restart: b1 * exp(b2 * x) under acat() takes b1 up and b2 down to the
  # linear limit, 5.9 nats below the maximum. The standard errors report that
  # flat direction under both BLASes. The return code is not read: rounding
  # sets it (1 under R's BLAS, 0 under OpenBLAS), and a cap on the iterations
  # sets it too, where a restart costs a factor of 4 and gains nothing
  flagged <- nl_flat_warning(held)
  bodies <- look$saddle
  if (!on_saddle && is.null(start) && flagged) {
    for (lp in frame[["linpreds"]]) {
      if (is.null(lp[["nl_body"]])) next
      bc <- nl_body_coefs(frame, lp)
      if (length(bc$coefs)) {
        bodies[[length(bodies) + 1L]] <-
          list(lp = lp, coefs = bc$coefs, seen = bc$seen)
      }
    }
  }
  if (!length(bodies)) {
    nl_local_message(frame, fit, look$local)
    release()
    return(fit)
  }
  key <- function(b) linpred_key(b$lp[["resp"]], b$lp[["dpar"]] %||% "")
  ll0 <- as.numeric(logLik_penalized(fit))
  # the values compared and printed below are what the fits maximized,
  # which under a prior with a density is not logLik()
  llw <- if (fit_penalized(fit)) "penalized log-likelihood" else
    "log-likelihood"
  first <- bodies
  best <- NULL
  best_ll <- ll0
  best_held <- list()
  left <- look$saddle
  # on a flat direction the first body's patterns move them all and settle it
  todo <- if (on_saddle) seq_along(bodies) else 1L
  pass_first <- TRUE
  ridge <- FALSE
  repeat {
    for (k in todo) {
      starts <- nl_saddle_starts(frame, bodies, start, k,
                                 if (ridge) NULL else best)
      # the first body's base pattern has already moved this one
      if (pass_first && k > 1L) starts <- starts[-1L]
      res <- nl_restart_round(run, starts, best_ll,
                              saddle_floor = on_saddle && is.null(best),
                              from_flat = !on_saddle && any(vapply(held, function(w) {
                                grepl("are not identified", conditionMessage(w))
                              }, NA)))
      if (is.null(res$best)) next
      best <- res$best
      best_ll <- res$ll
      best_held <- res$held
    }
    ridge <- FALSE
    if (is.null(best)) break
    look <- nl_saddle_find(best, use_curvature = use_curv)
    left <- look$saddle
    fresh <- Filter(function(b) {
      !key(b) %in% vapply(bodies, key, "")
    }, left)
    # a kept restart that still warns can have a body it did not move on a
    # ridge: `b1 * exp(b2 * x)` beside a saddle on `disc` ran to the linear
    # limit from its zeros, -377.03 at code 1 against -371.54 from a move
    # (dev/fitedge-findings.md, item 2). Every body moves then, since the kept
    # ones sit where the ridge left them: disc at d2 = 1499. Only a warning
    # that the stop is not the maximum asks for it: a product body always holds
    # "Standard errors are not available", and on two gaussian responses with
    # four bodies the pass cost 1.5 times for nothing
    # (dev/reviews/2026-11-06-fitedge.md, m2)
    if (!length(fresh) && nl_stop_warning(best_held)) {
      fresh <- Filter(function(b) {
        !key(b) %in% vapply(bodies, key, "")
      }, nl_all_bodies(frame))
      ridge <- length(fresh) > 0L
    }
    if (!length(fresh)) break
    todo <- length(bodies) + seq_along(fresh)
    bodies <- c(bodies, fresh)
    on_saddle <- TRUE
    pass_first <- FALSE
  }
  if (!is.null(best)) {
    wd <- nl_bodies_words(first)
    ll <- nl_ll_pair(ll0, best_ll)
    frm_message(
      if (length(nl_saddle_named(first, fit))) {
        paste0("The fit stopped where ", wd$noun, " ", wd$is, " constant ",
               "across the rows, ", llw, " ", ll[1L], ", so frm() ")
      } else {
        paste0("The fit of ", wd$noun, " from its default start ran along a ",
               "flat direction to ", llw, " ", ll[1L], ", so frm() ")
      },
      "restarted from moved coefficients and reached a higher point, ",
      llw, " ", ll[2L], ", optimizer return code ",
      format(best$opt$convergence), ". Give `start` to choose your own",
      class = "frmtmb_nl_restart")
    nl_local_message(frame, best, look$local)
    for (w in best_held) warning(w)
    if (length(left)) nl_saddle_warning(frame, left, best_ll, llw)
    return(best)
  }
  if (!length(left)) {
    release()
    return(fit)
  }
  nl_local_message(frame, fit, look$local)
  release()
  nl_saddle_warning(frame, left, ll0, llw)
  fit
}

#' Every body with coefficients, in nl_saddle_find()'s form.
#'
#' @noRd
nl_all_bodies <- function(frame) {
  out <- list()
  for (lp in frame[["linpreds"]]) {
    if (is.null(lp[["nl_body"]])) next
    bc <- nl_body_coefs(frame, lp)
    if (length(bc$coefs)) {
      out[[length(out) + 1L]] <- list(lp = lp, coefs = bc$coefs, seen = bc$seen)
    }
  }
  out
}

#' Which of `bodies` are constant at `fit`'s estimates: a restart's message
#' says so only of those, and not of a body it moved from a flat direction.
#'
#' @noRd
nl_saddle_named <- function(bodies, fit) {
  Filter(function(b) {
    e0 <- tryCatch(suppressWarnings(nl_body_at(fit$frame, b$lp, fit$estimates)),
                   error = function(e) NULL)
    !is.null(e0) && all(is.finite(e0)) &&
      diff(range(e0)) <= 1e-8 * max(1, abs(mean(e0)))
  }, bodies)
}

#' The warning on a fit that a moved start could not take off a constant body.
#'
#' @noRd
nl_saddle_warning <- function(frame, bodies, ll, llw = "log-likelihood") {
  w <- nl_bodies_words(bodies)
  seen <- unique(unlist(lapply(bodies, function(b) b$seen)))
  frm_warning(
    "The fit stops where ", w$noun, " ", w$is, " constant across the rows, ",
    "which is ", nl_model_words(frame, bodies), ": ", llw, " ",
    format(ll, digits = 7), ". Moved starts of ",
    paste0("'", seen, "'", collapse = ", "), " reached nothing higher, so the ",
    "data may simply have no effect there. A `start` in a region where the ",
    "body varies can test it",
    call. = FALSE, class = "frmtmb_nl_saddle")
}

#' The message for a constant body that the curvature calls a local maximum:
#' correct, and worth saying, since a thresholds-only fit passes as one with
#' code 0 and nothing else to mark it (`b^2 * x * (1 - b^2)`, 28.8 nats below
#' brms's mode on one data set).
#'
#' @noRd
nl_local_message <- function(frame, fit, bodies) {
  if (!length(bodies)) return(invisible(NULL))
  w <- nl_bodies_words(bodies)
  frm_message(
    "The fit stops where ", w$noun, " ", w$is, " constant across the rows, ",
    "which is ", nl_model_words(frame, bodies), ": log-likelihood ",
    format(as.numeric(stats::logLik(fit)), digits = 7), ". The point is a ",
    "local maximum along ", w$its, " coefficients. A `start` elsewhere can ",
    "test for a higher one",
    class = "frmtmb_nl_local")
  invisible(NULL)
}

#' The natural unit of each nonlinear coefficient, aligned with
#' `par_names` (`names(obj$par)`), or `NULL` when every unit is 1.
#'
#' A linear coefficient is stepped in the response's unit over its
#' column's spread (predictor_par_units()). A body has no columns, so the
#' spread is read from the body itself: the root mean square over the
#' rows of the body's derivative in the coefficient, at the start `tpl`.
#' Without it `a * x` with x about 1e-20 stopped at its start, 26.8 below
#' `y ~ 0 + x`, at code 0 with no warning, because the gradient was 1e-18.
#' The response's unit enters on the same terms as for a linear `mu`
#' (predictor_unit_input()), with the response's own spread, as a body
#' has no design to take out.
#'
#' Only a coefficient the body reads as an added term (nl_body_additive(),
#' `a` and `b` of `a + b * x`, `a` of `a * x`) gets a unit, because only
#' there is the derivative at the start the derivative everywhere. A unit
#' is a power of two, and one within a factor of 16 of 1 is not applied,
#' so a body whose derivatives are of order 1 keeps its path exactly. A
#' position `have` already scales (autoscale_units(), smooth_fx_units())
#' keeps that.
#'
#' A prior on the coefficient adds its curvature (nl_prior_precision()) to
#' the likelihood's, per row: the unit is `1 / sqrt(r^2 / resp_u^2 + p / n)`
#' for derivative spread `r`, prior precision `p` and `n` rows, which is
#' `resp_u / r` with no prior. Read from the likelihood alone, the unit of
#' `a` in `a * x` with x about 1e-20 under `normal(0, 1e-3)` was 1.5e20, the
#' prior's curvature dominated the optimizer's coordinates by about 1e40,
#' and 15 of 63 prior rows stopped at code 1 below 0.79.1
#' (dev/reviews/2026-11-06-nlstart.md, blocker 1).
#'
#' @noRd
nl_par_units <- function(frame, tpl, par_names, have = NULL,
                         prior_entries = NULL, curve = FALSE) {
  unit <- lapply(frame[["par_template"]], function(v) rep(NA_real_, length(v)))
  hit <- FALSE
  for (lp in frame[["linpreds"]]) {
    if (is.null(lp[["nl_body"]])) next
    coefs <- nl_body_coefs(frame, lp)$coefs
    if (!length(coefs)) next
    resp_u <- nl_resp_unit(frame, lp)
    body_at <- function(t) {
      tryCatch(suppressWarnings(nl_body_at(frame, lp, t)),
               error = function(e) NULL)
    }
    # only a coefficient the body reads as an added term has one derivative
    # everywhere; elsewhere the start's derivative is local, and a unit read
    # from it sent an ODE rate's late-phase shift to 9.1 against log(2),
    # where the late phase is 0 and flat (frmtmb.ode's test-ode-tv.R)
    added <- nl_body_additive(frame, lp, tpl, coefs)
    k <- 0L
    for (par in names(coefs)) {
      for (at in coefs[[par]]) {
        k <- k + 1L
        if (!added[k]) next
        h <- 1e-4 * max(1, abs(tpl[[par]][at]))
        up <- tpl
        dn <- tpl
        up[[par]][at] <- up[[par]][at] + h
        dn[[par]][at] <- dn[[par]][at] - h
        e1 <- body_at(up)
        e0 <- body_at(dn)
        if (is.null(e1) || is.null(e0) || length(e1) != length(e0)) next
        r <- sqrt(mean(((e1 - e0) / (2 * h))^2))
        if (!is.finite(r) || r <= 0) next
        p <- nl_prior_precision(prior_entries, par, at)
        if (is.na(p)) next
        u <- 2^round(log2(1 / sqrt(r^2 / resp_u^2 + p / length(e1))))
        # an order of magnitude is the failure this cures; a factor of 2
        # only changes the path, and on a body with more than one maximum
        # that can change which one it reaches
        if (!is.finite(u) || u <= 0 || abs(log2(u)) < 4) next
        # a coefficient two bodies read is stepped for the one it moves most
        unit[[par]][at] <- min(unit[[par]][at], u, na.rm = TRUE)
        hit <- TRUE
      }
    }
    # in the refit of a probe's winner (`curve`), a coefficient inside a
    # curve is stepped in the size the probe started it at
    # (nl_probe_scale()): with x about 1e-6, b2 near 1e6 stepped in absolute
    # units left 17 of 40 fits 0.67 to 58.6 below the same model at x about
    # 1 (short check 2, R3). Not on a first fit: there the size is read away
    # from the start, and frmtmb.ode's test-ode-tv.R (lk's late shift ran to
    # a flat phase) and test-open-issues.R's nlme check (xmid, a gradient of
    # 59380) failed with it. A prior on it leaves it unscaled
    if (!curve) next
    probed <- list()
    k <- 0L
    for (par in names(coefs)) {
      for (at in coefs[[par]]) {
        k <- k + 1L
        if (added[k]) next
        if (!identical(nl_prior_precision(prior_entries, par, at), 0)) next
        probed[[length(probed) + 1L]] <- list(par = par, at = at)
      }
    }
    if (!length(probed)) next
    sc <- nl_probe_scale(frame, lp, tpl, probed)
    for (i in seq_along(probed)) {
      u <- sc["mag", i]
      if (sc["curve", i] != 1 || is.na(u)) next
      q <- probed[[i]]
      unit[[q$par]][q$at] <- min(unit[[q$par]][q$at], u, na.rm = TRUE)
      hit <- TRUE
    }
  }
  if (!hit) return(NULL)
  out <- numeric(0)
  for (cp in unique(par_names)) {
    u <- unit[[cp]]
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      u <- u[-frame[["betad_fixed_idx"]]]
    }
    out <- c(out, u)
  }
  if (length(out) != length(par_names)) return(NULL)
  out[is.na(out)] <- 1
  if (length(have) == length(out)) out[have != 1] <- 1
  if (all(out == 1)) return(NULL)
  out
}

#' The response's unit for a body, on the terms predictor_units() sets for
#' a linear `mu`: the sd of the response on `mu`'s scale for an
#' identity-link location-scale family, else 1.
#'
#' @noRd
nl_resp_unit <- function(frame, lp) {
  plain <- lp
  plain[["nl_body"]] <- NULL
  inp <- tryCatch(predictor_unit_input(frame, plain),
                  error = function(e) NULL)
  if (is.null(inp)) return(1)
  s <- stats::sd(inp$z[inp$ok])
  if (is.finite(s) && s > 0) s else 1
}

#' The curvature at its centre of the log prior density on template
#' position `at` of component `par`, summed over the prior entries that
#' reach it on the internal scale: `1 / s^2` for a normal, `(df + 1) /
#' (df s^2)` for a Student t and `1 / (2 s^2)` for a logistic of scale `s`,
#' and 0 for an exponential, a uniform or a flat prior, whose log density
#' has no curvature. 0 where no prior reaches it.
#'
#' `NA` where any other prior reaches it (gamma, lognormal, beta, ...), or
#' one on another scale: its curvature depends on where it is read, and the
#' caller then gives no unit. Read as 0, a `gamma(2, 1000)` or
#' `lognormal(0, 0.001)` prior on `a` in `a * x` with x about 1e-20 left the
#' likelihood's unit of 1.5e20, and 8 rows stopped at code 1 below 0.79.1, by
#' up to 992 (dev/reviews/2026-11-06-nlstart.md, short check, R1).
#'
#' @noRd
nl_prior_precision <- function(prior_entries, par, at) {
  p <- 0
  for (e in prior_entries %||% list()) {
    if (!identical(e$comp, par) || !at %in% e$idx) next
    if (!identical(e$scale %||% "internal", "internal")) return(NA_real_)
    d <- e$dist
    kind <- d$kind %||% ""
    if (kind %in% c("exponential", "uniform", "flat")) next
    if (!kind %in% c("normal", "t", "logistic")) return(NA_real_)
    s <- d$scale
    if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s <= 0) {
      return(NA_real_)
    }
    if (kind == "t" && !(is.numeric(d$df) && isTRUE(d$df > 0))) {
      return(NA_real_)
    }
    p <- p + switch(kind,
                    normal = 1 / s^2,
                    t = (d$df + 1) / (d$df * s^2),
                    logistic = 1 / (2 * s^2))
  }
  p
}

#' Which coefficients a body reads as an added term, `a` or `b3 * z`: the
#' body's derivative in the coefficient is the same at two points with
#' every coefficient nonzero and different, so that `b1` of `b1 * exp(b2 *
#' x)`, linear alone but scaled by `b2`, is not one. The sign a start gives
#' such a term cannot decide which maximum a fit reaches. A logical vector
#' over `coefs` in its order; a body that will not evaluate answers `FALSE`
#' throughout.
#'
#' @noRd
nl_body_additive <- function(frame, lp, tpl, coefs) {
  n <- sum(lengths(coefs))
  base <- tpl
  for (par in names(coefs)) base[[par]][coefs[[par]]] <- 0.5
  pts <- list(base, nl_moved_tpl(base, coefs, 0.3))
  at_body <- function(t) {
    v <- tryCatch(suppressWarnings(nl_body_at(frame, lp, t)),
                  error = function(e) NULL)
    if (is.null(v) || !all(is.finite(v))) NULL else v
  }
  out <- logical(0)
  for (par in names(coefs)) {
    for (at in coefs[[par]]) {
      d <- lapply(pts, function(t) {
        up <- t
        dn <- t
        up[[par]][at] <- up[[par]][at] + 1e-3
        dn[[par]][at] <- dn[[par]][at] - 1e-3
        e1 <- at_body(up)
        e0 <- at_body(dn)
        if (is.null(e1) || is.null(e0)) NULL else (e1 - e0) / 2e-3
      })
      if (is.null(d[[1L]]) || is.null(d[[2L]])) return(rep(FALSE, n))
      out <- c(out, max(abs(d[[2L]] - d[[1L]])) <=
                 1e-6 * max(abs(d[[1L]]), 1e-300))
    }
  }
  out
}

#' The sign patterns frm() starts a nonlinear body's coefficients from:
#' every one of the coefficients at positions `lead` (at most four), the
#' others positive. With four coefficients or fewer that is every pattern,
#' and the result does not depend on their order. A cap of four patterns
#' picked by position gave one model two answers by the order of its terms
#' (bernoulli, two saturating terms: -163.5320 and -162.0191,
#' dev/reviews/2026-11-06-nlstart.md, major 1); nl_sign_probe() picks
#' `lead` by what the body does, not by position.
#'
#' @noRd
nl_sign_patterns <- function(k, lead = seq_len(min(k, 4L))) {
  m <- length(lead)
  lapply(seq_len(2^m), function(i) {
    s <- rep(1, k)
    s[lead] <- ifelse(bitwAnd(i - 1L, 2^(seq_len(m) - 1L)) > 0, -1, 1)
    s
  })
}

#' Each probed coefficient's start size and its ranking key, from the body.
#' `mag` is a power of two where it is 16 or more times from 1, else `NA`,
#' and the caller keeps 1 over the column's half range. `key` is the root
#' mean square of the body's derivative in the coefficient, `-Inf` where it
#' does not evaluate. Both are set by the body and the data alone, with every
#' coefficient sized at once, so neither the order of the terms nor the
#' names change them.
#'
#' Two kinds of coefficient need two sizes. One the body reads inside a
#' curve (`b2` of `b1 * exp(b2 * x)`) is sized so that the curve's argument
#' moves by 3 on the row it moves most: 3 over the largest absolute
#' derivative where it is 0 and the other probed coefficients are 1. On `x`
#' in 0 to 3000 that is about 1e-3; at 1 the body overflowed and every
#' pattern failed, 12.8 to 33.0 short on 5 of 5 seeds
#' (dev/reviews/2026-11-06-nlstart.md, major 2). One the body is linear in
#' (`b1`) is then sized so that it moves the body by two of the response's
#' units on the average row, read with the curve's coefficients at their
#' own sizes. Read with them at 1, `b1 * (1 - exp(-b2 * x))` with x about
#' 1e-4 gave `b1` the size of `b2`, 32768, and the fit stayed 74 to 86 short
#' on 5 of 5 gaussian seeds (short check, R2). Sizing both kinds alike, a
#' pass that takes each size from the others' does not settle: `b1 b2`
#' moves the body as one product, so the two sizes trade places.
#'
#' @noRd
nl_probe_scale <- function(frame, lp, tpl, probed) {
  resp_u <- nl_resp_unit(frame, lp)
  n <- length(probed)
  body_with <- function(vals, q = NULL, h = NULL) {
    t <- tpl
    for (i in seq_len(n)) t[[probed[[i]]$par]][probed[[i]]$at] <- vals[i]
    if (!is.null(q)) t[[probed[[q]]$par]][probed[[q]]$at] <- h
    v <- tryCatch(suppressWarnings(nl_body_at(frame, lp, t)),
                  error = function(e) NULL)
    if (is.null(v) || !all(is.finite(v))) NULL else v
  }
  deriv <- function(vals, q) {
    e1 <- body_with(vals, q, 1e-6)
    e0 <- body_with(vals, q, -1e-6)
    if (is.null(e1) || is.null(e0) || length(e1) != length(e0)) return(NULL)
    (e1 - e0) / 2e-6
  }
  # the 16-fold zone is read before rounding: rounded first, a location
  # inside a logistic curve sized 12 became 16 and started at +-16, and 5 of
  # 5 gamma fits of the review's logistic stayed 9.7 to 23.4 short, silent
  zone <- function(m) {
    if (!is.finite(m) || m <= 0 || abs(log2(m)) < 4) return(NA_real_)
    2^round(log2(m))
  }
  ones <- rep(1, n)
  d1 <- lapply(seq_len(n), function(q) deriv(ones, q))
  key <- vapply(d1, function(d) {
    r <- if (is.null(d)) NA_real_ else sqrt(mean(d^2))
    if (is.finite(r) && r > 0) r else -Inf
  }, 0)
  curve <- vapply(d1, function(d) {
    if (is.null(d)) return(NA_real_)
    zone(3 / max(abs(d)))
  }, 0)
  s1 <- ifelse(is.na(curve), 1, curve)
  # linear in the coefficient: its second difference is 0 at the sizes
  linear <- vapply(seq_len(n), function(q) {
    e <- lapply(c(0, 1, 2), function(k) body_with(s1, q, k * s1[q]))
    if (any(vapply(e, is.null, NA))) return(FALSE)
    d1q <- max(abs(e[[2L]] - e[[1L]]))
    d1q > 0 &&
      max(abs(e[[3L]] - 2 * e[[2L]] + e[[1L]])) <= 1e-6 * d1q
  }, NA)
  mag <- ifelse(linear, NA_real_, curve)
  st <- ifelse(is.na(mag), 1, mag)
  for (q in which(linear)) {
    d <- deriv(st, q)
    r <- if (is.null(d)) NA_real_ else sqrt(mean(d^2))
    if (!is.finite(r) || r <= 0) next
    mag[q] <- zone(2 * resp_u / r)
    key[q] <- r
  }
  rbind(mag = mag, key = key, curve = as.numeric(!linear))
}

#' Optimize a fitted nonlinear model again from every sign pattern of its
#' bodies' coefficients, on the fit's own objective, and return a full fit
#' from the best point when that is higher than `fit`; else `NULL`.
#'
#' A body that is not affine in its coefficients can have more than one
#' maximum, and a fit from the default start reaches the one its path
#' finds. On 450 nonlinear rows (9 families, 5 bodies, 2 slopes, 5 seeds;
#' dev/nlstart-scan.R) 16 defaults were short of the best of 8 random
#' starts by 0.43 to 74.8, 11 of them with no warning: `b1 * (1 - exp(-b2
#' * x))` under gaussian ran to its linear limit (b1 near 5e4), 70 below
#' the maximum at b1 and b2 both negative, which no restart of
#' nl_fit_off_saddle() tries; a bernoulli `b1 * exp(b2 * x)` stopped at
#' the other sign of b2. Nothing in such a fit marks it, so the only test
#' is to look. The patterns vary the coefficients a body does not read as
#' an added term (nl_body_additive()), all of them up to four, else the four
#' with the largest derivative (nl_sign_patterns()). Each coefficient starts
#' at the size nl_probe_scale() reads from the body, or else at 1 over its
#' column's half range, the rest of the model where the fit started. Each
#' pattern runs the optimizer on the tape the fit already has, so it costs
#' an optimization and no tape or standard error. Only a winner by more
#' than nl_probe_gain() is fitted again from its point, by `run()`.
#'
#' Not run where a coefficient is outside the outer vector (a map, a
#' profiled or whitened `beta`), under quadrature or importance
#' sampling, with a custom optimizer or one asked to stop early
#' (optctrl_stops_earlier()), when every coefficient is an added term
#' (`a + b * x`, whose only maxima are the family's), or under
#' `options(frmtmb.nl_sign_probe = FALSE)`.
#'
#' @noRd
nl_sign_probe <- function(run, frame, fit) {
  obj <- fit$obj
  ctl <- fit$control
  if (!isTRUE(getOption("frmtmb.nl_sign_probe", TRUE)) ||
        is.null(obj) || is.null(fit$opt) || !is.finite(fit$opt$objective) ||
        isTRUE(fit$quadrature) || !is.null(fit$importance) ||
        !is.null(fit$beta_whiten) || isTRUE(ctl$profile) ||
        is.function(ctl$optimizer) ||
        !identical(ctl$optimizer %||% "nlminb", "nlminb") ||
        optctrl_stops_earlier(ctl$optCtrl)) {
    return(NULL)
  }
  pn <- names(obj$par)
  nl_named <- unique(unlist(lapply(frame[["spec"]]$responses,
                                   function(r) r$nlpars %||% character(0))))
  if (any(unique(pn) %in% nl_named)) return(NULL)
  est <- fit$estimates
  lps <- frame[["linpreds"]]
  at_obj <- integer(0)
  width <- numeric(0)
  mag_b <- numeric(0)
  key <- numeric(0)
  val <- numeric(0)
  nm <- character(0)
  for (lp in lps) {
    if (is.null(lp[["nl_body"]])) next
    bc <- nl_body_coefs(frame, lp)
    if (!length(bc$coefs)) next
    added <- nl_body_additive(frame, lp, est, bc$coefs)
    probed <- list()
    k <- 0L
    for (par in names(bc$coefs)) {
      pos <- which(pn == par)
      if (length(pos) != length(frame[["par_template"]][[par]])) return(NULL)
      for (at in bc$coefs[[par]]) {
        k <- k + 1L
        if (!added[k]) probed[[length(probed) + 1L]] <- list(par = par,
                                                             at = at,
                                                             j = pos[at])
      }
    }
    if (!length(probed)) next
    sc <- nl_probe_scale(frame, lp, est, probed)
    for (i in seq_along(probed)) {
      q <- probed[[i]]
      hit <- match(q$j, at_obj)
      # a coefficient two bodies read starts at the smaller size
      if (!is.na(hit)) {
        two <- c(mag_b[hit], sc["mag", i])
        mag_b[hit] <- if (all(is.na(two))) NA_real_ else min(two, na.rm = TRUE)
        key[hit] <- max(key[hit], sc["key", i])
        next
      }
      w <- 1
      for (qq in lps) {
        col <- match(q$at, qq[["idx"]])
        if (!identical(qq[["par"]], q$par) || is.na(col) ||
              !is.null(qq[["nl_body"]])) next
        r <- diff(range(as.matrix(qq[["X"]])[, col]))
        if (is.finite(r) && r > 0) w <- max(1, r / 2)
        break
      }
      at_obj <- c(at_obj, q$j)
      width <- c(width, w)
      mag_b <- c(mag_b, sc["mag", i])
      key <- c(key, sc["key", i])
      val <- c(val, est[[q$par]][q$at])
      nm <- c(nm, paste(q$par, par_template_names(
        frame[["par_template"]][[q$par]], q$par)[q$at]))
    }
  }
  if (!length(at_obj)) return(NULL)
  pu <- fit$par_units
  u <- if (length(pu) == length(pn)) pu[at_obj] else rep(1, length(at_obj))
  mag <- ifelse(is.na(mag_b), u / width, mag_b)
  # with more than four, the four the body moves most by; a tie (terms that
  # read one covariate alike) goes by the fitted values and then the name,
  # not by the order of the formula, which changed 8 of 36 fits by up to
  # 0.147 when the terms were reversed (short check, minor)
  lead <- utils::head(order(-key, -mag, -abs(val), -val, nm), 4L)
  pats <- nl_sign_patterns(length(at_obj), lead)
  # the runs step a curve's coefficient in its start size, as the winner's
  # refit does (nl_par_units(curve = TRUE)); the ratio of the two calls is
  # that factor alone, everything the fit already scales staying as it was
  pe <- if (!is.null(fit$prior) && NROW(fit$prior)) {
    tryCatch(resolve_prior_input(list(frame = frame, spec = fit$spec),
                                 fit$prior)$entries,
             error = function(e) NULL)
  }
  cu <- nl_par_units(frame, est, pn, pu, pe, curve = TRUE) %||% 1
  c0 <- nl_par_units(frame, est, pn, pu, pe, curve = FALSE) %||% 1
  pu_run <- if (all(cu / c0 == 1)) pu else (pu %||% 1) * (cu / c0)
  bounds <- fit_outer_box(fit)
  saved <- obj_state_save(obj)
  # a basis from nl_basis_descent() whitens the body's Jacobian at the
  # first fit's stop; a pattern starts at other signs, where it does not
  # hold, and with it no pattern won on rows the probe exists for
  # (dev/rel080-nlrows.R). The runs step in the design's coordinates.
  ob_fit <- obj$env$.frmtmb_opt_basis
  if (isTRUE(ob_fit$nl)) {
    obj$env$.frmtmb_opt_basis <- opt_basis(frame, obj, ctl$optimizer)
  }
  best <- fit$opt
  for (s in pats) {
    p <- obj$par
    p[at_obj] <- s * mag
    p <- pmin(pmax(p, bounds$lower), bounds$upper)
    o <- tryCatch(suppressWarnings(optimize_obj(obj, ctl, bounds, pu_run,
                                                start_par = p)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$objective) &&
          o$objective < best$objective) {
      best <- o
    }
  }
  obj$env$.frmtmb_opt_basis <- ob_fit
  obj_settle(obj, fit$opt, fit$opt, saved)
  gain <- fit$opt$objective - best$objective
  if (!(gain > nl_probe_gain(fit$opt$objective))) return(NULL)
  st <- list()
  for (cp in unique(pn)) {
    v <- est[[cp]]
    pos <- seq_along(v)
    if (cp == "betad" && length(frame[["betad_fixed_idx"]])) {
      pos <- setdiff(pos, frame[["betad_fixed_idx"]])
    }
    v[pos] <- unname(best$par[pn == cp])
    st[[cp]] <- unname(v)
  }
  # read by fit_assembled(): step the curves' coefficients as the run did
  attr(st, "frmtmb_nl_probe") <- TRUE
  held <- new.env(parent = emptyenv())
  held$w <- list()
  r <- tryCatch(
    withCallingHandlers(suppressMessages(run(st, FALSE)),
                        warning = function(w) {
                          held$w[[length(held$w) + 1L]] <- w
                          invokeRestart("muffleWarning")
                        }),
    error = function(e) NULL)
  if (is.null(r)) return(NULL)
  ll0 <- as.numeric(logLik_penalized(fit))
  ll1 <- suppressWarnings(as.numeric(logLik_penalized(r)))
  if (!is.finite(ll1) || !(ll1 - ll0 > nl_probe_gain(ll0))) {
    return(NULL)
  }
  # did the refit settle where the pattern's run stopped, at code 0? A refit
  # that warns says so itself; one that is silent and did not settle is
  # told by nl_fit_probed(), so that a fit never comes back quieter than it
  # is (short check 2, R3)
  settled <- length(held$w) > 0L ||
    (identical(as.integer(r$opt$convergence), 0L) &&
       is.finite(r$opt$objective) &&
       r$opt$objective <= best$objective + nl_probe_gain(best$objective))
  list(fit = r, held = held$w, ll0 = ll0, ll1 = ll1, n = length(pats),
       settled = settled, ll_tape = -best$objective)
}

#' The gain a probe's point needs over the first fit before frm() refits
#' from it and says so: 0.01 log-likelihood units, or 1e-9 of the value
#' when that is larger. On the 450 rows of dev/nlstart-scan.R every gain
#' was 0.43 or more; on the review's 335 fresh rows three were 9.6e-4 to
#' 2.0e-3, each a slide along a ridge to the linear limit that left the fit
#' 3.0 to 11.2 short, which a refit and a message do not change.
#'
#' @noRd
nl_probe_gain <- function(v) max(1e-2, 1e-9 * abs(v))

#' `frm()`'s fit of a model with a nonlinear body: nl_fit_off_saddle(),
#' then, from the default start only, nl_sign_probe(). What the first
#' says is held, and said only if its fit is the one returned; a message
#' about where the fit started stays true either way and is kept.
#'
#' @noRd
nl_fit_probed <- function(run, frame, start, prior) {
  held <- new.env(parent = emptyenv())
  held$c <- list()
  keep <- function(cnd) {
    held$c[[length(held$c) + 1L]] <- cnd
    if (inherits(cnd, "warning")) invokeRestart("muffleWarning")
    invokeRestart("muffleMessage")
  }
  say <- function(cnds) {
    for (cnd in cnds) {
      if (inherits(cnd, "warning")) warning(cnd) else message(cnd)
    }
  }
  fit <- tryCatch(
    withCallingHandlers(nl_fit_off_saddle(run, frame, start, prior),
                        warning = keep, message = keep),
    error = function(e) {
      say(held$c)
      stop(e)
    })
  pr <- if (is.null(start)) nl_sign_probe(run, frame, fit)
  if (is.null(pr)) {
    say(held$c)
    return(fit)
  }
  # what the first fit's stop was, or a restart it kept, is no longer the
  # answer
  say(Filter(function(cnd) {
    !inherits(cnd, c("warning", "frmtmb_nl_local", "frmtmb_nl_restart"))
  }, held$c))
  ll <- nl_ll_pair(pr$ll0, pr$ll1)
  # the values are what the fits maximized, which under a prior with a
  # density is not logLik()
  llw <- if (fit_penalized(pr$fit)) "penalized log-likelihood" else
    "log-likelihood"
  frm_message(
    "The fit of the nonlinear model from its default start reached ",
    llw, " ", ll[1L], ". Of ", pr$n, " starts with the signs of the ",
    "body's coefficients varied, the best reached ", llw, " ", ll[2L],
    ", optimizer return code ", format(pr$fit$opt$convergence), ", and frm() ",
    "returns that fit. Give `start` to choose your own",
    class = "frmtmb_nl_multistart")
  say(pr$held)
  if (!pr$settled) {
    lt <- nl_ll_pair(pr$ll_tape, pr$ll1)
    frm_warning(
      "The refit from the best start pattern did not settle where that ",
      "pattern's run stopped: it ended at ", llw, " ", lt[2L],
      ", optimizer return code ", format(pr$fit$opt$convergence),
      ", against ", lt[1L], " on the first fit's objective. The fit may ",
      "not be at a maximum; a `start` near the estimates can test it",
      call. = FALSE, class = "frmtmb_nl_multistart_unsettled")
  }
  pr$fit
}

#' Move a body the link cannot evaluate at its start to the response mean.
#'
#' A nonlinear coefficient starts at 0, so a body such as `b0 + b1 * x`
#' starts at 0 too, and the inverse link of Gamma() and the 1/mu^2 link of
#' inverse.gaussian() are undefined there (as a link with a bound at 0
#' is, on a family that needs a positive mean). The fit then stopped at
#' its first evaluation where the linear fit of the same model starts at
#' the family's initializer. frm() has the data, so one intercept of a
#' parameter the body reads moves to where the body's mean equals
#' link(mean(y)). The parameters are tried in the order the body names
#' them, each with the others at their starts, and the first that makes
#' the objective finite is kept: `b0 + b1 * x` takes b0, since b1 x
#' changes sign on these rows.
#'
#' This is one linear solve per parameter and not a search of the
#' objective, which `nl_start_message()` explains the danger of. It runs
#' only where the objective at the start is not finite and the body
#' itself is, so every fit that could evaluate its start keeps it, and a
#' body undefined at 0, such as `log(a)`, still gets the error that asks
#' for a start.
#'
#' `taken` holds the `beta` positions `start` or a located prior set.
#'
#' @noRd
nl_domain_start <- function(frame, tpl, taken, announce, nll = NULL) {
  bodies <- Filter(function(lp) {
    if (is.null(lp[["nl_body"]])) return(FALSE)
    fam <- frame[["spec"]]$responses[[lp[["resp"]]]]$family
    # a latent or category predictor has the whole real line for its domain
    lp[["dpar"]] %in% (fam[["primary_dpars"]] %||% "mu") &&
      !isTRUE(fam[["type"]] %in% c("ordinal", "categorical")) &&
      is.numeric(frame[["y"]][[lp[["resp"]]]]) &&
      is.null(dim(frame[["y"]][[lp[["resp"]]]]))
  }, frame[["linpreds"]])
  if (!length(bodies)) return(tpl)
  nll <- nll %||% tryCatch(build_objective(frame), error = function(e) NULL)
  if (is.null(nll)) return(tpl)
  value <- function(t) {
    tryCatch(suppressWarnings(as.numeric(nll(t))), error = function(e) NA_real_)
  }
  if (is.finite(value(tpl))) return(tpl)
  # per body, the moves to try; a body left alone is the first of each, so
  # a multivariate model moves only the bodies that need it
  moves <- lapply(bodies, function(lp) {
    e0 <- tryCatch(suppressWarnings(nl_body_at(frame, lp, tpl)),
                   error = function(e) NULL)
    if (is.null(e0) || !all(is.finite(e0))) return(list(NULL))
    y <- frame[["y"]][[lp[["resp"]]]]
    target <- suppressWarnings(lp[["link"]]$linkfun(mean(y, na.rm = TRUE)))
    if (!is.finite(target)) return(list(NULL))
    out <- list(NULL)
    for (np in lp[["nl_pars"]]) {
      q <- frame[["linpreds"]][[linpred_key(lp[["resp"]], np)]]
      if (!identical(q[["par"]], "beta")) next
      at <- q[["idx"]][match("(Intercept)", colnames(q[["X"]]))]
      if (!length(at) || is.na(at) || at %in% taken) next
      t1 <- tpl
      t1[["beta"]][at] <- t1[["beta"]][at] + 1
      e1 <- tryCatch(suppressWarnings(nl_body_at(frame, lp, t1)),
                     error = function(e) NULL)
      if (is.null(e1)) next
      slope <- mean(e1 - e0)
      if (!is.finite(slope) || slope == 0) next
      v <- tpl[["beta"]][at] + (target - mean(e0)) / slope
      out[[length(out) + 1L]] <- list(at = at, value = v, lp = lp)
    }
    out
  })
  ways <- as.matrix(expand.grid(lapply(moves, seq_along)))
  for (i in seq_len(min(nrow(ways), 16L))[-1L]) {
    cand <- tpl
    done <- list()
    for (j in seq_along(moves)) {
      m <- moves[[j]][[ways[i, j]]]
      if (is.null(m)) next
      cand[["beta"]][m$at] <- m$value
      done[[length(done) + 1L]] <- m
    }
    if (length(done) && is.finite(value(cand))) {
      if (announce) {
        nms <- par_template_names(frame[["par_template"]][["beta"]], "beta")
        frm_message(
          "Nonlinear starting values moved to where the body equals the ",
          "response mean, because the link cannot evaluate it at 0: ",
          paste(vapply(done, function(m) {
            paste0(nms[m$at], " = ", format(m$value, digits = 4))
          }, ""), collapse = ", "), ". Give `start` to choose your own.")
      }
      return(cand)
    }
  }
  tpl
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
