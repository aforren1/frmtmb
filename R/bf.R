#' Set up a model formula
#'
#' Specifies a model with brms-compatible syntax. Distributional
#' parameters (dpars) can get their own formulas with the full predictor
#' grammar, or be fixed to constants:
#' `bf(y ~ x + (1 | g), sigma ~ z + (1 | g))`, `bf(y ~ x, sigma = 1)`.
#' Nonlinear formulas (`nl = TRUE`) and multivariate models (see
#' [mvbf()]) use the same grammar.
#'
#' The left-hand side accepts addition terms after `|`:
#' `y | weights(w) ~ ...` and `y | trials(n) ~ ...`.
#' Every linear predictor accepts lme4-style random effects `(1 | g)`,
#' `(1 + x | g)`, `(x || g)`, and explicit covariance-structure wrappers
#' `us(x | g)` and `diag(x | g)`.
#'
#' Attach a family with `+`, for example `bf(y ~ x) + gaussian()`.
#'
#' @param formula The model formula for `mu`.
#' @param ... Two-sided formulas for other dpars (the left-hand side names
#'   the dpar, e.g. `sigma ~ z`), or named scalars fixing a dpar to a
#'   constant on the response scale (e.g. `sigma = 1`).
#' @param family Optional family; can also be attached with `+` or passed
#'   to [frm()].
#' @param nl Nonlinear-formula flag: the main formula becomes a
#'   nonlinear expression of named parameters, each given its own
#'   `...` formula with the full predictor grammar.
#' @return An object of class `frmtmb_formula`.
#' @examples
#' # brms-style model formulas: attach a family with `+`
#' bf(y ~ x + (1 | g)) + gaussian()
#' # distributional parameters get their own formulas or constants
#' bf(y ~ x, sigma ~ x)
#' bf(y ~ x, shape = 2) + Gamma()
#' # nonlinear models declare parameter formulas and nl = TRUE
#' bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1 + (1 | g), nl = TRUE)
#' @export
bf <- function(formula, ..., family = NULL, nl = FALSE) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula", call. = FALSE)
  }
  # mvbind(y1, y2) ~ rhs: shared predictors, one bf per response
  if (length(formula) == 3L && is.call(formula[[2]]) &&
      identical(formula[[2]][[1]], as.name("mvbind"))) {
    resps <- as.list(formula[[2]])[-1]
    forms <- lapply(resps, function(r) {
      f1 <- formula
      f1[[2]] <- r
      bf(f1, ..., family = family)
    })
    return(do.call(mvbf, forms))
  }
  dots <- list(...)
  pforms <- list()
  pfix <- list()
  for (i in seq_along(dots)) {
    d <- dots[[i]]
    nm <- names(dots)[i] %||% ""
    if (inherits(d, "formula")) {
      if (length(d) != 3L) {
        stop("dpar formulas must be two-sided, naming the dpar on the ",
             "left: e.g. sigma ~ x", call. = FALSE)
      }
      dpar <- deparse1(d[[2]])
      if (!grepl("^[a-zA-Z][a-zA-Z0-9]*$", dpar)) {
        stop("Invalid parameter name '", dpar, "': names must be ",
             "alphanumeric without dots or underscores (they collide ",
             "with coefficient naming)", call. = FALSE)
      }
      if (dpar %in% c(names(pforms), names(pfix))) {
        stop("Duplicated dpar: '", dpar, "'", call. = FALSE)
      }
      pforms[[dpar]] <- d
    } else if (is.numeric(d) && length(d) == 1L) {
      if (nm == "") {
        stop("Constant dpar values must be named, e.g. bf(y ~ x, sigma = 1)",
             call. = FALSE)
      }
      if (nm %in% c(names(pforms), names(pfix))) {
        stop("Duplicated dpar: '", nm, "'", call. = FALSE)
      }
      pfix[[nm]] <- d
    } else {
      stop("Cannot interpret bf() argument ",
           if (nm != "") paste0("'", nm, "'") else i,
           ": expected a dpar formula or a named numeric constant",
           call. = FALSE)
    }
  }
  if (isTRUE(nl) && !length(pforms)) {
    stop("nl = TRUE needs at least one parameter formula, e.g. ",
         "bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)", call. = FALSE)
  }
  structure(
    list(formula = formula, pforms = pforms, pfix = pfix, nl = isTRUE(nl),
         family = if (!is.null(family)) as_frmtmb_family(family)),
    class = "frmtmb_formula"
  )
}

#' @export
"+.frmtmb_formula" <- function(e1, e2) {
  if (missing(e2)) return(e1)
  if (inherits(e2, "frmtmb_formula")) {
    return(mvbf(e1, e2))
  }
  if (inherits(e2, "frmtmb_rescor")) {
    stop("set_rescor() applies to multivariate formulas; combine ",
         "responses with mvbf() or `bf() + bf()` first", call. = FALSE)
  }
  if (inherits(e2, "frmtmb_family") || inherits(e2, "family") ||
      is.function(e2)) {
    e1$family <- as_frmtmb_family(e2)
    return(e1)
  }
  stop("Cannot add an object of class ", paste(class(e2), collapse = "/"),
       " to a 'frmtmb_formula'", call. = FALSE)
}

#' Combine formulas into a multivariate model
#'
#' Each response keeps its own formula, family, dpar formulas, and
#' addition terms. Residual correlation between gaussian responses is
#' requested with `rescor = TRUE` or [set_rescor()]. Random-effect
#' correlation across responses uses the brms `|ID|` syntax, e.g.
#' `(1 | p | g)` in several formulas correlates their `g` effects.
#'
#' @param ... `bf()` formulas, each with a family attached (or supply one
#'   `family` to [frm()] for all of them).
#' @param rescor Model residual correlation between the responses
#'   (gaussian only).
#' @return An object of class `frmtmb_mvformula`.
#' @export
mvbf <- function(..., rescor = FALSE) {
  forms <- list(...)
  flat <- list()
  for (f in forms) {
    if (inherits(f, "frmtmb_mvformula")) {
      flat <- c(flat, f$forms)
      if (isTRUE(f$rescor)) rescor <- TRUE
    } else if (inherits(f, "frmtmb_formula")) {
      flat <- c(flat, list(f))
    } else {
      stop("mvbf() takes bf() formulas", call. = FALSE)
    }
  }
  if (length(flat) < 2) {
    stop("mvbf() needs at least two responses", call. = FALSE)
  }
  structure(list(forms = flat, rescor = isTRUE(rescor)),
            class = "frmtmb_mvformula")
}

#' @rdname mvbf
#' @param rescor_value For `set_rescor()`: turn residual correlation on
#'   or off.
#' @export
set_rescor <- function(rescor_value = TRUE) {
  structure(list(rescor = isTRUE(rescor_value)), class = "frmtmb_rescor")
}

#' @export
"+.frmtmb_mvformula" <- function(e1, e2) {
  if (missing(e2)) return(e1)
  if (inherits(e2, "frmtmb_rescor")) {
    e1$rescor <- e2$rescor
    return(e1)
  }
  if (inherits(e2, "frmtmb_mvformula")) {
    return(mvbf(e1, e2, rescor = e1$rescor))
  }
  if (inherits(e2, "frmtmb_formula")) {
    return(mvbf(e1, e2, rescor = e1$rescor))
  }
  if (inherits(e2, "frmtmb_family") || inherits(e2, "family") ||
      is.function(e2)) {
    fam <- as_frmtmb_family(e2)
    e1$forms <- lapply(e1$forms, function(f) {
      if (is.null(f$family)) f$family <- fam
      f
    })
    return(e1)
  }
  stop("Cannot add an object of class ", paste(class(e2), collapse = "/"),
       " to a 'frmtmb_mvformula'", call. = FALSE)
}

#' @export
print.frmtmb_mvformula <- function(x, ...) {
  for (f in x$forms) print(f)
  cat("rescor:", x$rescor, "\n")
  invisible(x)
}

#' @export
print.frmtmb_formula <- function(x, ...) {
  cat(deparse1(x$formula), "\n")
  for (f in x$pforms) cat(deparse1(f), "\n")
  for (nm in names(x$pfix)) cat(nm, "=", x$pfix[[nm]], "\n")
  if (!is.null(x$family)) {
    cat("Family:", x$family$family, "\n")
  }
  invisible(x)
}

# Normalize a plain formula or bf()/mvbf() object plus an optional
# family argument into a bform with families attached: the argument
# fills empty per-response slots of a multivariate form and overrides
# a univariate one. Shared by frm(), get_prior(), and frm_simulate()
# so the coercion cannot drift between entry points.
as_bform <- function(formula, family = NULL) {
  bform <- if (inherits(formula,
                        c("frmtmb_formula", "frmtmb_mvformula"))) {
    formula
  } else {
    bf(formula)
  }
  if (!is.null(family)) {
    fam <- as_frmtmb_family(family)
    if (inherits(bform, "frmtmb_mvformula")) {
      bform$forms <- lapply(bform$forms, function(f) {
        if (is.null(f$family)) f$family <- fam
        f
      })
    } else {
      bform$family <- fam
    }
  }
  bform
}
