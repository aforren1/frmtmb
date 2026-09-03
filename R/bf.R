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
#' Attach a family with `+`, for example `bf(y ~ x) + gaussian()`, or
#' pass one to [frm()]. A model that names no family is gaussian.
#'
#' @param formula The model formula for `mu`.
#' @param ... Two-sided formulas for other dpars (the left-hand side
#'   names the dpar, e.g. `sigma ~ z`, or several sharing one
#'   right-hand side, e.g. `b1 + b2 ~ 1`), or named scalars fixing a
#'   dpar to a constant on the response scale (e.g. `sigma = 1`).
#' @param family Optional family; can also be attached with `+` or
#'   passed to [frm()], which uses `gaussian()` when nothing names one.
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
#' @srrstats {G2.0} Inputs expected to be single-valued are asserted to be
#'   so. A distributional parameter fixed to a constant must satisfy
#'   `is.numeric(d) && length(d) == 1L`; the tuning arguments of the
#'   special terms (`gp()`, `rr()`, `car()`, `mm()`, `gr()`) go through
#'   one validator that errors with the argument name and the length it
#'   received when a scalar was expected.
#' @srrstats {G2.1} Inputs are asserted to be of the expected type.
#'   `formula` must inherit from `"formula"` and errors otherwise; the
#'   special-term validator requires finite numeric tuning values; and
#'   multipliers in a `mo()` or `mi()` interaction must be numeric, with
#'   character and factor columns both refused rather than coerced to
#'   all-`NA` and surfacing later as an optimizer failure.
#' @srrstats {G2.14c} Missing predictor values can be replaced by imputed
#'   ones inside the model rather than dropped. `bf(x | mi() ~ ...)`
#'   declares an imputation model for a partially observed variable, and
#'   `mi(x)` uses it in another formula, so the missing values become
#'   latent parameters estimated jointly with everything else.
#'   `bf(x | mi(sdx) ~ ...)` supplies known measurement standard
#'   deviations the same way. For imputation performed outside the
#'   model, [frm_multiple()] fits each completed data set and pools by
#'   Rubin's rules.
#' @srrstats {RE2.2} Missing values in the response and in the predictors
#'   are processed differently and separately. Rows with a missing
#'   predictor are removed by `na.action` before fitting, so a model can
#'   be fitted on complete predictor data and used to generate values for
#'   every associated response point; missing values in a variable
#'   declared with `mi()` are not removed but imputed in-model. This is
#'   why frame assembly builds the `mi()` branch with `stats::na.pass`
#'   and drops rows only on the non-`mi()` columns.
#'
#' @export
bf <- function(formula, ..., family = NULL, nl = FALSE) {
  if (inherits(formula, c("brmsformula", "bform"))) {
    stop("this formula was built by brms::bf(): attaching brms after ",
         "frmtmb masks frmtmb's bf(), so a bare bf() call now reaches ",
         "brms. Call frmtmb::bf() explicitly, or attach brms before ",
         "frmtmb", call. = FALSE)
  }
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
      for (dpar in lhs_dpar_names(d[[2]])) {
        if (dpar %in% c(names(pforms), names(pfix))) {
          stop("Duplicated dpar formula: '", dpar, "'", call. = FALSE)
        }
        di <- d
        di[[2]] <- as.name(dpar)
        pforms[[dpar]] <- di
      }
    } else if (is.numeric(d) && length(d) == 1L) {
      if (nm == "") {
        stop("Constant dpar values must be named, e.g. bf(y ~ x, sigma = 1)",
             call. = FALSE)
      }
      if (nm %in% c(names(pforms), names(pfix))) {
        stop("Duplicated dpar constant: '", nm, "'", call. = FALSE)
      }
      pfix[[nm]] <- d
    } else {
      stop("Cannot interpret bf() argument ",
           if (nm != "") paste0("'", nm, "'") else i,
           ": expected a dpar formula or a named numeric constant",
           call. = FALSE)
    }
  }
  # A body that is one bare name (`bf(y ~ a, nl = TRUE)`) is the brms
  # nlf() spelling, where the parameter formulas arrive afterwards with
  # `+ nlf(a ~ ...)`; anything else with no formula here is the usual
  # slip of forgetting them, and is worth catching at the call.
  if (isTRUE(nl) && !length(pforms) &&
      !is.name(reformulas::RHSForm(formula))) {
    stop("nl = TRUE needs at least one parameter formula, e.g. ",
         "bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE). Formulas ",
         "added afterwards with lf() or nlf() are not visible here, so ",
         "give bf() at least one of them", call. = FALSE)
  }
  structure(
    list(formula = formula, pforms = pforms, pfix = pfix, nl = isTRUE(nl),
         nlforms = list(),
         family = if (!is.null(family)) as_frmtmb_family(family)),
    class = "frmtmb_formula"
  )
}

#' The left-hand side of a dpar or nonlinear-parameter formula, checked
#' once for `bf()` and `lf()`: it becomes part of a coefficient name, so
#' dots and underscores are refused rather than allowed to collide with
#' the `dpar_term` separator.
#'
#' @noRd
check_dpar_name <- function(dpar) {
  if (!grepl("^[a-zA-Z][a-zA-Z0-9]*$", dpar)) {
    stop("Invalid parameter name '", dpar, "': names must be ",
         "alphanumeric without dots or underscores (they collide ",
         "with coefficient naming)", call. = FALSE)
  }
  dpar
}

#' The parameter names on the left of one dpar formula. brms lets a
#' single formula name several nonlinear parameters (`b1 + b2 ~ 1`);
#' each name then gets its own copy of the right-hand side, which is
#' the same model as writing the formulas out.
#'
#' @noRd
lhs_dpar_names <- function(lhs) {
  if (is.call(lhs) && identical(lhs[[1L]], as.name("+")) &&
      length(lhs) == 3L) {
    return(c(lhs_dpar_names(lhs[[2L]]), lhs_dpar_names(lhs[[3L]])))
  }
  check_dpar_name(deparse1(lhs))
}

#' Add parameter formulas to a model formula
#'
#' brms's `lf()`: one or more two-sided formulas for distributional or
#' nonlinear parameters, added to a [bf()] with `+`. It is sugar for
#' passing the same formulas to `bf()` directly, and is useful when the
#' parameter formulas are built somewhere else than the response
#' formula.
#'
#' `bf(y ~ x) + lf(sigma ~ z)` and `bf(y ~ x, sigma ~ z)` give the same
#' model. In a multivariate model an `lf()` must be added to the `bf()`
#' of the response it belongs to, before the responses are combined.
#'
#' @param ... Two-sided formulas naming the parameter on the left, e.g.
#'   `sigma ~ x` or (with `nl = TRUE` on the `bf()`) a nonlinear
#'   parameter's formula `a ~ 1 + (1 | g)`.
#' @return An object of class `frmtmb_lf`, to be added to a [bf()].
#' @examples
#' # the two spellings are the same model
#' bf(y ~ x) + lf(sigma ~ z)
#' bf(y ~ x, sigma ~ z)
#'
#' # nonlinear parameter formulas can arrive the same way
#' bf(y ~ a * exp(-b * x), a ~ 1, nl = TRUE) + lf(b ~ 1 + (1 | g))
#' @export
lf <- function(...) {
  dots <- list(...)
  pforms <- list()
  for (d in dots) {
    if (!inherits(d, "formula") || length(d) != 3L) {
      stop("lf() takes two-sided formulas naming the parameter on the ",
           "left: e.g. lf(sigma ~ x)", call. = FALSE)
    }
    for (dpar in lhs_dpar_names(d[[2]])) {
      if (dpar %in% names(pforms)) {
        stop("Duplicated parameter formula in lf(): '", dpar, "'",
             call. = FALSE)
      }
      di <- d
      di[[2]] <- as.name(dpar)
      pforms[[dpar]] <- di
    }
  }
  if (!length(pforms)) {
    stop("lf() needs at least one parameter formula, e.g. lf(sigma ~ x)",
         call. = FALSE)
  }
  structure(list(pforms = pforms), class = "frmtmb_lf")
}

#' @export
print.frmtmb_lf <- function(x, ...) {
  for (f in x$pforms) cat(deparse1(f), "\n")
  invisible(x)
}

#' Add a nonlinear parameter formula to a model formula
#'
#' brms's `nlf()`: it declares that one parameter is computed by a
#' NONLINEAR expression rather than by a linear predictor. The names in
#' that expression are either other model parameters, each of which
#' needs its own formula, or columns of the data.
#'
#' `bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)` is the same model as
#' `bf(y ~ exp(b * x), b ~ 1, nl = TRUE)`, written the other way round.
#' Where `nl = TRUE` makes the response formula the nonlinear body,
#' `nlf()` names the parameter it belongs to, so any parameter can have
#' one: `bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1)` is a
#' nonlinear model for the residual standard deviation with a linear
#' `mu`, which `nl = TRUE` cannot spell.
#'
#' Bodies may be chained: an `nlf()` body can name a parameter that
#' another `nlf()` defines, to any depth. The parameters are evaluated
#' in dependency order, and a cycle is refused by name. Each
#' parameter's link is applied to the body's value, so
#' `nlf(sigma ~ a + b * z)` gives `sigma = exp(a + b * z)` under the
#' default log link, exactly as in brms.
#'
#' A body may also name another distributional parameter of the same
#' response, which reads that parameter's per-row VALUE. This is the one
#' place frmtmb goes beyond brms, where such a name is a data column and
#' the model is refused when no column has it. It buys the variance
#' function of the model's own mean that nlme writes as
#' `varPower(form = ~ fitted(.))`: with the default log link on `sigma`,
#' `nlf(sigma ~ ls + th * log(abs(mu)))` is
#' `sd = exp(ls) * |mu|^th`. A column of the data still wins over the
#' parameter name, so a body ported from brms keeps its meaning.
#'
#' Like [lf()], an `nlf()` in a multivariate model must be added to the
#' `bf()` of the response it belongs to, before the responses are
#' combined.
#'
#' @param formula A two-sided formula naming the parameter on the left
#'   and its nonlinear body on the right, e.g. `sigma ~ a * exp(b * x)`.
#' @param ... Further two-sided formulas, treated as LINEAR parameter
#'   formulas exactly as if passed to [lf()] - the brms convention.
#' @param loop Accepted for brms source compatibility and ignored.
#'   frmtmb evaluates a nonlinear body once over whole vectors, which is
#'   brms's `loop = FALSE`; a body built from elementwise operations has
#'   the same value either way.
#' @return An object of class `frmtmb_nlf`, to be added to a [bf()].
#' @examples
#' # the composed spelling of a nonlinear model
#' bf(y ~ a) + nlf(a ~ exp(b * x)) + lf(b ~ 1)
#'
#' # a nonlinear sigma with a linear mu
#' bf(y ~ x) + nlf(sigma ~ a + b * z) + lf(a ~ 1, b ~ 1)
#'
#' # bodies chain: cc feeds a, a feeds mu
#' bf(y ~ a, nl = TRUE) + nlf(a ~ cc * x) + nlf(cc ~ exp(b)) + lf(b ~ 1)
#'
#' # a variance function of the fitted mean: sd = exp(ls) * |mu|^th
#' bf(y ~ x) + nlf(sigma ~ ls + th * log(abs(mu))) + lf(ls ~ 1, th ~ 1)
#' @export
nlf <- function(formula, ..., loop = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("nlf() takes a two-sided formula naming the parameter on the ",
         "left: e.g. nlf(sigma ~ a * exp(b * x))", call. = FALSE)
  }
  lhs <- formula[[2L]]
  if (is.call(lhs) && identical(lhs[[1L]], as.name("+"))) {
    stop("nlf() declares one parameter at a time, and '", deparse1(lhs),
         "' names several. Sharing one nonlinear body would make them ",
         "the same function of the data, which leaves the model ",
         "aliased; write one nlf() per parameter", call. = FALSE)
  }
  dpar <- check_dpar_name(deparse1(lhs))
  nlforms <- list()
  nlforms[[dpar]] <- formula
  pforms <- if (length(list(...))) lf(...)$pforms else list()
  if (dpar %in% names(pforms)) {
    stop("nlf() gives '", dpar, "' both a nonlinear body and a linear ",
         "formula; it can have one or the other", call. = FALSE)
  }
  structure(list(nlforms = nlforms, pforms = pforms),
            class = "frmtmb_nlf")
}

#' @export
print.frmtmb_nlf <- function(x, ...) {
  for (f in x$nlforms) cat(deparse1(f), " (nonlinear)\n", sep = "")
  for (f in x$pforms) cat(deparse1(f), "\n")
  invisible(x)
}

#' @export
"+.frmtmb_formula" <- function(e1, e2) {
  if (missing(e2)) return(e1)
  if (inherits(e2, "frmtmb_formula")) {
    return(mvbf(e1, e2))
  }
  if (inherits(e2, "frmtmb_lf")) {
    for (nm in names(e2$pforms)) {
      if (nm %in% c(names(e1$pforms), names(e1$pfix),
                    names(e1$nlforms))) {
        stop("lf() sets '", nm, "', which the bf() it is added to ",
             "already sets", call. = FALSE)
      }
      e1$pforms[[nm]] <- e2$pforms[[nm]]
    }
    return(e1)
  }
  if (inherits(e2, "frmtmb_nlf")) {
    set <- function(nm) {
      nm %in% c(names(e1$pforms), names(e1$pfix), names(e1$nlforms))
    }
    for (nm in names(e2$nlforms)) {
      if (set(nm)) {
        stop("nlf() sets '", nm, "', which the bf() it is added to ",
             "already sets", call. = FALSE)
      }
      e1$nlforms[[nm]] <- e2$nlforms[[nm]]
    }
    for (nm in names(e2$pforms)) {
      if (set(nm)) {
        stop("The linear parameter formulas passed to nlf() set '", nm,
             "', which the bf() it is added to already sets",
             call. = FALSE)
      }
      e1$pforms[[nm]] <- e2$pforms[[nm]]
    }
    return(e1)
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
#' The linked terms merge into one covariance block, so they must all
#' name the same grouping specification. When they all write
#' `gr(g, cov = A)` (or all `gr(g, prec = Q)`) with the same matrix, the
#' merged block keeps it: its covariance is `A (x) Sigma`, with `Sigma`
#' unstructured across the merged coefficients. A two-trait animal model
#' is therefore the same fit whether written across two responses with
#' `(1 | q | gr(id, cov = A))` or in long format as a single
#' `(0 + trait | gr(id, cov = A))`. Mixing structures under one key -
#' a plain `g` in one formula and `gr(g, cov = A)` in another, or
#' `cov =` against `prec =` - is refused, because a merged block has
#' room for one structure.
#'
#' @param ... `bf()` formulas, each with a family attached (or supply one
#'   `family` to [frm()] for all of them).
#' @param rescor Model residual correlation between the responses
#'   (gaussian only).
#' @return An object of class `frmtmb_mvformula`.
#' @examples
#' set.seed(2)
#' n <- 160
#' dd <- data.frame(x = rnorm(n), g = factor(rep(1:16, 10)))
#' u <- cbind(rnorm(16, 0, 0.8), rnorm(16, 0, 0.8))
#' e <- rnorm(n)                      # a disturbance both responses see
#' dd$y1 <- 1 + 0.5 * dd$x + u[dd$g, 1] + e + rnorm(n, 0, 0.5)
#' dd$y2 <- 2 - 0.3 * dd$x + u[dd$g, 2] + e + rnorm(n, 0, 0.5)
#'
#' # each response keeps its own formula and family
#' fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
#' fixef(fit)
#'
#' # rescor estimates the correlation of the residuals
#' fit_rc <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
#'               data = dd)
#' rescor_matrix(fit_rc)
#'
#' # set_rescor() turns it on after the fact, and `+` also combines bf()s
#' mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE)
#' bf(y1 ~ x) + bf(y2 ~ x)
#'
#' # |ID| correlates the random effects of the two responses
#' fit_id <- frm(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
#'                 gaussian(), data = dd)
#' VarCorr(fit_id)
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
  if (inherits(e2, "frmtmb_lf")) {
    stop("lf() does not say which response it belongs to. Put it ",
         "directly after the bf() it modifies, e.g. ",
         "bf(y1 ~ x) + lf(sigma ~ z) + bf(y2 ~ x)", call. = FALSE)
  }
  if (inherits(e2, "frmtmb_nlf")) {
    stop("nlf() does not say which response it belongs to. Put it ",
         "directly after the bf() it modifies, e.g. ",
         "bf(y1 ~ x) + nlf(sigma ~ a * z) + bf(y2 ~ x)", call. = FALSE)
  }
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
  cat(deparse1(x$formula), if (isTRUE(x$nl)) " (nonlinear)" else "", "\n",
      sep = "")
  for (f in x$nlforms) cat(deparse1(f), " (nonlinear)\n", sep = "")
  for (f in x$pforms) cat(deparse1(f), "\n")
  for (nm in names(x$pfix)) cat(nm, "=", x$pfix[[nm]], "\n")
  if (!is.null(x$family)) {
    cat("Family:", x$family$family, "\n")
  }
  invisible(x)
}

#' Normalize a plain formula or `bf()`/`mvbf()` object plus an optional
#' family argument into a bform with families attached: the argument
#' fills empty per-response slots of a multivariate form and overrides
#' a univariate one. A response left without a family after that gets
#' `gaussian()`, the brms/lme4/glmmTMB convention. Shared by `frm()`,
#' `get_prior()`, and `frm_simulate()` so the coercion cannot drift
#' between entry points.
#'
#' @noRd
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
  # Silent, like glmmTMB and lme4: a message here would fire on every
  # linear mixed model, which is the most common fit there is.
  if (inherits(bform, "frmtmb_mvformula")) {
    bform$forms <- lapply(bform$forms, function(f) {
      if (is.null(f$family)) f$family <- default_family()
      f
    })
  } else if (is.null(bform$family)) {
    bform$family <- default_family()
  }
  bform
}

#' The family a model gets when neither the `family` argument nor a `+`
#' attachment names one.
#'
#' @noRd
default_family <- function() as_frmtmb_family(stats::gaussian())
