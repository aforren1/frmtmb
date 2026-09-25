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
#' @param formula The model formula for `mu`, or a formula `bf()`
#'   already built. Given one of those and nothing else, `bf()` returns
#'   it unchanged, as brms's `bf()` does; further arguments add to it.
#' @param ... Two-sided formulas for other dpars (the left-hand side
#'   names the dpar, e.g. `sigma ~ z`, or several sharing one
#'   right-hand side, e.g. `b1 + b2 ~ 1`), one-sided formulas named by
#'   their dpar (`sigma = ~ z`, the same formula), or named scalars
#'   fixing a dpar to a constant on the response scale (e.g.
#'   `sigma = 1`).
#' @param family Optional family; can also be attached with `+` or
#'   passed to [frm()], which uses `gaussian()` when nothing names one.
#' @param nl Nonlinear-formula flag: the main formula becomes a
#'   nonlinear expression of named parameters, each given its own
#'   `...` formula with the full predictor grammar. `NULL`, the default
#'   as in brms, means `FALSE` for a new formula and keeps the setting
#'   of a formula `bf()` already built.
#' @param center Whether the location formula's intercept is brms's
#'   class `"Intercept"`, a density on the intercept at the means of the
#'   predictors. `NULL`, the default, means `TRUE` for a new formula and
#'   keeps the setting of a formula `bf()` already built. `FALSE` makes
#'   the intercept an ordinary coefficient, class `"b"` with coef
#'   `"Intercept"`, as `0 + Intercept` in the formula does; see the
#'   section on brms's reserved `Intercept` below. It changes priors
#'   only: the likelihood and the maximum likelihood fit are the same.
#'   It applies to the location formula alone, as in brms; give a
#'   parameter formula its own with [lf()].
#' @return An object of class `frmtmb_formula`.
#' @section brms's reserved `Intercept`:
#' In a formula without an intercept, `Intercept` is a reserved name, as
#' in brms: `y ~ 0 + Intercept + x` is the model `y ~ 1 + x`, with the
#' intercept as an ordinary population-level coefficient. Factors get
#' treatment contrasts, as they do beside an intercept. The likelihood
#' is the same, so a maximum likelihood fit is the same; what changes
#' is the prior. brms places a class `"Intercept"` prior at the means
#' of the predictors, and a class `"b"` prior on this coefficient at
#' zero, and so does frmtmb: `set_prior(..., class = "b")` reaches it,
#' `class = "Intercept"` does not. `bf(y ~ x, center = FALSE)` gives
#' the same model. The spelling works in every linear formula: the
#' location, a distributional parameter (`sigma ~ 0 + Intercept + z`)
#' and a nonlinear parameter (`a ~ 0 + Intercept + x`, where it changes
#' nothing, because brms never centers a nonlinear parameter).
#'
#' `Intercept` must be a term of its own; `Intercept:x` is refused.
#' An ordinal family refuses `0 + Intercept`, as brms does, because its
#' thresholds take the intercept's place. A data column named
#' `Intercept` must hold only ones in such a model, as brms requires;
#' in a formula with an intercept, `Intercept` is an ordinary variable.
#' Prediction needs no `Intercept` column in `newdata`.
#' @examples
#' # brms-style model formulas: attach a family with `+`
#' bf(y ~ x + (1 | g)) + gaussian()
#' # distributional parameters get their own formulas or constants
#' bf(y ~ x, sigma ~ x)
#' bf(y ~ x, shape = 2) + Gamma()
#' # nonlinear models declare parameter formulas and nl = TRUE
#' bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1 + (1 | g), nl = TRUE)
#' # an intercept that a class "b" prior reaches; the two are one model
#' bf(y ~ 0 + Intercept + x)
#' bf(y ~ x, center = FALSE)
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
bf <- function(formula, ..., family = NULL, nl = NULL, center = NULL) {
  if (inherits(formula, c("brmsformula", "bform"))) {
    frm_stop("this formula was built by brms::bf(): attaching brms after ",
             "frmtmb masks frmtmb's bf(), so a bare bf() call now reaches ",
             "brms. Call frmtmb::bf() explicitly, or attach brms before ",
             "frmtmb", call. = FALSE)
  }
  # brms's bf() takes its own output and returns it unchanged when
  # nothing else is given, so code that normalizes a formula by calling
  # bf() on it must not change a model already built
  existing <- inherits(formula, "frmtmb_formula")
  if (!existing && !inherits(formula, "formula")) {
    frm_stop("`formula` must be a formula", call. = FALSE)
  }
  # nl reaches isTRUE() at the end of this function, which reads "yes",
  # NA and c(TRUE, FALSE) as FALSE: bf(..., nl = "yes") used to build a
  # LINEAR model and say nothing. NULL keeps what an existing bf() says
  if (!is.null(nl)) check_flag(nl, "nl")
  if (!is.null(center)) check_flag(center, "center")
  if (existing) {
    return(bf_update(formula, ..., family = family, nl = nl,
                     center = center))
  }
  refuse_nested_formula(formula)
  # mvbind(y1, y2) ~ rhs: shared predictors, one bf per response
  if (length(formula) == 3L && is.call(formula[[2]]) &&
      identical(formula[[2]][[1]], as.name("mvbind"))) {
    resps <- as.list(formula[[2]])[-1]
    forms <- lapply(resps, function(r) {
      f1 <- formula
      f1[[2]] <- r
      bf(f1, ..., family = family, nl = nl, center = center)
    })
    return(do.call(mvbf, forms))
  }
  parsed <- bf_dots(list(...))
  pforms <- parsed[["pforms"]]
  pfix <- parsed[["pfix"]]
  # A body that is one bare name (`bf(y ~ a, nl = TRUE)`) is the brms
  # nlf() spelling, where the parameter formulas arrive afterwards with
  # `+ nlf(a ~ ...)`; anything else with no formula here is the usual
  # slip of forgetting them, and is worth catching at the call.
  if (isTRUE(nl) && !length(pforms) &&
      !is.name(reformulas::RHSForm(formula))) {
    frm_stop("nl = TRUE needs at least one parameter formula, e.g. ",
             "bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE). Formulas ",
             "added afterwards with lf() or nlf() are not visible here, so ",
             "give bf() at least one of them", call. = FALSE)
  }
  out <- structure(
    list(formula = formula, pforms = pforms, pfix = pfix, nl = isTRUE(nl),
         nlforms = list(),
         family = if (!is.null(family)) as_frmtmb_family(family)),
    class = c("frmtmb_formula", "frmtmb_bform")
  )
  if (!is.null(center)) out$center <- center
  out
}

#' Read `bf()`'s dots into dpar formulas and constants, onto whatever an
#' existing formula already sets. `taken` names the parameters a
#' nonlinear body already defines, which a new entry may not reuse.
#'
#' @noRd
bf_dots <- function(dots, pforms = list(), pfix = list(),
                    taken = character(0)) {
  for (i in seq_along(dots)) {
    d <- dots[[i]]
    nm <- names(dots)[i] %||% ""
    if (inherits(d, "formula")) {
      d <- named_par_formula(d, nm)
      refuse_nested_formula(d)
      for (dpar in lhs_dpar_names(d[[2]])) {
        if (dpar %in% c(names(pforms), names(pfix), taken)) {
          frm_stop("Duplicated dpar formula: '", dpar, "'", call. = FALSE)
        }
        di <- d
        di[[2]] <- as.name(dpar)
        pforms[[dpar]] <- di
      }
    } else if (is.numeric(d) && length(d) == 1L) {
      if (nm == "") {
        frm_stop(
          "Constant dpar values must be named, e.g. bf(y ~ x, sigma = 1)",
          call. = FALSE)
      }
      if (nm %in% c(names(pforms), names(pfix), taken)) {
        frm_stop("Duplicated dpar constant: '", nm, "'", call. = FALSE)
      }
      pfix[[nm]] <- d
    } else {
      frm_stop("Cannot interpret bf() argument ",
               if (nm != "") paste0("'", nm, "'") else i,
               ": expected a dpar formula or a named numeric constant",
               call. = FALSE)
    }
  }
  list(pforms = pforms, pfix = pfix)
}

#' `bf()` on a formula `bf()` already built. With nothing else given the
#' object comes back untouched, which is brms's contract
#' (`identical(form, bf(form))`). Further dpar formulas and constants are
#' added; one naming a parameter the formula already sets is refused, as
#' `+ lf()` refuses it, rather than replacing it in silence.
#'
#' @noRd
bf_update <- function(formula, ..., family = NULL, nl = NULL,
                      center = NULL) {
  dots <- list(...)
  if (!length(dots) && is.null(family) && is.null(nl) && is.null(center)) {
    return(formula)
  }
  if (length(dots)) {
    parsed <- bf_dots(dots, formula[["pforms"]], formula[["pfix"]],
                      taken = names(formula[["nlforms"]]))
    formula[["pforms"]] <- parsed[["pforms"]]
    formula[["pfix"]] <- parsed[["pfix"]]
  }
  if (!is.null(nl)) formula[["nl"]] <- isTRUE(nl)
  if (!is.null(center)) formula[["center"]] <- center
  if (!is.null(family)) formula[["family"]] <- as_frmtmb_family(family)
  formula
}

#' A parameter formula in its two-sided form. brms takes the parameter
#' either on the left, `sigma ~ x`, or as the argument name of a
#' one-sided formula, `sigma = ~ x`; both are the same formula. A
#' one-sided formula with no name belongs to no parameter.
#'
#' @noRd
named_par_formula <- function(d, nm) {
  if (length(d) == 3L) return(d)
  if (!nzchar(nm)) {
    frm_stop("Additional formulas must be named: '", deparse1(d), "' does ",
             "not say which parameter it belongs to. Write it two-sided, ",
             "as sigma ~ x, or name it, as sigma = ~ x", call. = FALSE)
  }
  structure(call("~", as.name(nm), d[[2L]]), class = "formula",
            .Environment = environment(d))
}

#' Refuse a formula whose right-hand side is itself a formula.
#'
#' `y ~ ~ x` parses, and until this refusal it fitted `y ~ x` without a
#' word, so a doubled tilde that lost a term in editing went unnoticed.
#' brms refuses the same shape by name (its issue #749) and checks the
#' top level only, which is what is checked here: a `~` nested inside
#' another call is left to whatever reads that call.
#'
#' @noRd
refuse_nested_formula <- function(f) {
  rhs <- f[[length(f)]]
  if (is.call(rhs) && identical(rhs[[1L]], as.name("~"))) {
    frm_stop("Nested formulas are not allowed: the right-hand side of '",
             deparse1(f), "' is itself a formula. Did you use '~~' ",
             "somewhere?", call. = FALSE)
  }
  invisible(f)
}

#' The left-hand side of a dpar or nonlinear-parameter formula, checked
#' once for `bf()` and `lf()`: it becomes part of a coefficient name, so
#' dots and underscores are refused rather than allowed to collide with
#' the `dpar_term` separator.
#'
#' @noRd
check_dpar_name <- function(dpar) {
  if (!grepl("^[a-zA-Z][a-zA-Z0-9]*$", dpar)) {
    frm_stop("Invalid parameter name '", dpar, "': names must be ",
             "alphanumeric and must not contain dots or underscores (they ",
             "collide with coefficient naming)", call. = FALSE)
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
#' model. In a multivariate model, add an `lf()` to the `bf()` of the
#' response it belongs to, or name that response with `resp =`:
#' `bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + lf(sigma ~ z, resp = "y3")`.
#'
#' @param ... Two-sided formulas naming the parameter on the left, e.g.
#'   `sigma ~ x` or (with `nl = TRUE` on the `bf()`) a nonlinear
#'   parameter's formula `a ~ 1 + (1 | g)`, or one-sided formulas named
#'   by their parameter, `sigma = ~ x`.
#' @param resp The response the formulas belong to, when the `lf()` is
#'   added to a multivariate formula. `NULL` (the default) adds them to
#'   the `bf()` on the left of the `+`, which must then be a single
#'   formula.
#' @param center `FALSE` makes the intercept of each of these formulas an
#'   ordinary coefficient, class `"b"` with coef `"Intercept"`, instead
#'   of brms's class `"Intercept"`: the same as `0 + Intercept` in the
#'   formula. See [bf()]. `NULL`, the default, leaves the intercept as
#'   class `"Intercept"`.
#' @return An object of class `frmtmb_lf`, to be added to a [bf()].
#' @examples
#' # the two spellings are the same model
#' bf(y ~ x) + lf(sigma ~ z)
#' bf(y ~ x, sigma ~ z)
#'
#' # nonlinear parameter formulas can arrive the same way
#' bf(y ~ a * exp(-b * x), a ~ 1, nl = TRUE) + lf(b ~ 1 + (1 | g))
#'
#' # in a multivariate formula, resp = says which response it modifies
#' bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + lf(sigma ~ z, resp = "y3")
#' @export
lf <- function(..., resp = NULL, center = NULL) {
  check_lf_resp(resp, "lf()")
  if (!is.null(center)) check_flag(center, "center")
  dots <- list(...)
  pforms <- list()
  for (i in seq_along(dots)) {
    d <- dots[[i]]
    if (!inherits(d, "formula")) {
      frm_stop("lf() takes two-sided formulas naming the parameter on the ",
               "left: e.g. lf(sigma ~ x)", call. = FALSE)
    }
    d <- named_par_formula(d, names(dots)[i] %||% "")
    refuse_nested_formula(d)
    for (dpar in lhs_dpar_names(d[[2]])) {
      if (dpar %in% names(pforms)) {
        frm_stop("Duplicated parameter formula in lf(): '", dpar, "'",
                 call. = FALSE)
      }
      di <- d
      di[[2]] <- as.name(dpar)
      # carried on the formula itself, which is what reaches the parser
      # through bf() and mvbf() unchanged
      if (!is.null(center)) attr(di, "center") <- center
      pforms[[dpar]] <- di
    }
  }
  if (!length(pforms)) {
    frm_stop("lf() needs at least one parameter formula, e.g. lf(sigma ~ x)",
             call. = FALSE)
  }
  structure(list(pforms = pforms, resp = resp), class = "frmtmb_lf")
}

#' `resp =` of lf() and nlf() is one response name.
#'
#' @noRd
check_lf_resp <- function(resp, fn) {
  if (is.null(resp)) return(invisible(NULL))
  if (!is.character(resp) || length(resp) != 1L || is.na(resp) ||
      !nzchar(resp)) {
    frm_stop(fn, ": resp must be a single response name, not ",
             arg_desc(resp), call. = FALSE)
  }
  invisible(NULL)
}

#' @export
print.frmtmb_lf <- function(x, ...) {
  frm_check_dots(...)
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
#' A function the body CALLS is looked up in RTMB before anywhere else,
#' so a bare `pnorm()` or `qgamma()` is the tape-capable version and
#' needs no `RTMB::` prefix. See the "What a nonlinear body sees"
#' section of [frm()] for the whole rule, including the `stats::` escape
#' hatch.
#'
#' Like [lf()], an `nlf()` in a multivariate model is added to the
#' `bf()` of the response it belongs to, or names that response with
#' `resp =`.
#'
#' @param formula A two-sided formula naming the parameter on the left
#'   and its nonlinear body on the right, e.g. `sigma ~ a * exp(b * x)`.
#' @param ... Further two-sided formulas, treated as LINEAR parameter
#'   formulas exactly as if passed to [lf()] - the brms convention.
#' @param resp The response the formulas belong to, when the `nlf()` is
#'   added to a multivariate formula, as in [lf()].
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
nlf <- function(formula, ..., resp = NULL, loop = NULL) {
  check_lf_resp(resp, "nlf()")
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    frm_stop("nlf() takes a two-sided formula naming the parameter on the ",
             "left: e.g. nlf(sigma ~ a * exp(b * x))", call. = FALSE)
  }
  refuse_nested_formula(formula)
  lhs <- formula[[2L]]
  if (is.call(lhs) && identical(lhs[[1L]], as.name("+"))) {
    frm_stop("nlf() declares one parameter at a time, and '", deparse1(lhs),
             "' names several. Sharing one nonlinear body would make them ",
             "the same function of the data, which leaves the model ",
             "aliased; write one nlf() per parameter", call. = FALSE)
  }
  dpar <- check_dpar_name(deparse1(lhs))
  nlforms <- list()
  nlforms[[dpar]] <- formula
  pforms <- if (length(list(...))) lf(...)$pforms else list()
  if (dpar %in% names(pforms)) {
    frm_stop("nlf() gives '", dpar, "' both a nonlinear body and a linear ",
             "formula; it can have one or the other", call. = FALSE)
  }
  structure(list(nlforms = nlforms, pforms = pforms, resp = resp),
            class = "frmtmb_nlf")
}

#' @export
print.frmtmb_nlf <- function(x, ...) {
  frm_check_dots(...)
  for (f in x$nlforms) cat(deparse1(f), " (nonlinear)\n", sep = "")
  for (f in x$pforms) cat(deparse1(f), "\n")
  invisible(x)
}

#' Add to a model formula
#'
#' One `+` method serves [bf()] and [mvbf()] objects alike, so a sum of
#' any number of formulas reads left to right. The sum of `bf(y1 ~ x)`,
#' `bf(y2 ~ x)` and `bf(y3 ~ x)` is
#' `mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x))`.
#'
#' The right-hand side can be another formula, a family, [lf()],
#' [nlf()] or [set_rescor()]. A family added to a single formula sets
#' its family. A family added to a multivariate formula goes to every
#' response that has no family yet, so in the sum of `bf(o ~ x)`,
#' `cumulative()`, `bf(y ~ x)` and `gaussian()` the response `o` stays
#' ordinal. brms instead gives the last family to every response. An
#' [lf()] or [nlf()] added to a multivariate formula names its response
#' with `resp =`.
#'
#' @param e1 A `bf()` or `mvbf()` object.
#' @param e2 The object to add.
#' @return A `frmtmb_formula` or a `frmtmb_mvformula`.
#' @examples
#' # three responses, one family for all of them
#' bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + gaussian()
#'
#' # a family per response, and a dpar formula for the third response
#' bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian() + bf(y2 ~ x) +
#'   lf(sigma ~ x, resp = "y2")
#' @name plus-bform
#' @export
"+.frmtmb_bform" <- function(e1, e2) {
  if (missing(e2)) return(e1)
  # R's Ops group dispatch calls this method when EITHER operand is a
  # bf() or mvbf() object. A single method for the shared parent class
  # is what lets a formula and a multivariate formula meet: with two
  # different methods R warns "Incompatible methods" and falls back to
  # the internal `+`. chooseOpsMethod() would settle that too, but only
  # from R 4.3.
  if (!inherits(e1, "frmtmb_bform")) {
    frm_stop("Cannot add a 'frmtmb_formula' to an object of class ",
             paste(class(e1), collapse = "/"), ". Start the sum with ",
             "the bf() formula, e.g. bf(y ~ x) + gaussian()", call. = FALSE)
  }
  if (inherits(e1, "frmtmb_mvformula")) plus_mvbf(e1, e2) else plus_bf(e1, e2)
}

#' The response label of one `bf()`, as `resp =` addresses it: the
#' left-hand side without its addition terms.
#'
#' @noRd
bform_resp_label <- function(f) {
  lhs <- f$formula[[2L]]
  if (is.call(lhs) && identical(lhs[[1L]], as.name("|"))) lhs <- lhs[[2L]]
  deparse1(lhs)
}

#' `resp =` of an lf() or nlf() matches a response by its name or by the
#' name brms makes of it (`y_a` is `ya`).
#'
#' @noRd
bform_resp_matches <- function(f, resp) {
  lab <- bform_resp_label(f)
  resp %in% c(lab, brms_stan_name(lab))
}

#' `bf() + e2`.
#'
#' @noRd
plus_bf <- function(e1, e2) {
  if (inherits(e2, "frmtmb_mvformula") || inherits(e2, "frmtmb_formula")) {
    return(mvbf(e1, e2))
  }
  if ((inherits(e2, "frmtmb_lf") || inherits(e2, "frmtmb_nlf")) &&
      !is.null(e2[["resp"]]) && !bform_resp_matches(e1, e2[["resp"]])) {
    fn <- if (inherits(e2, "frmtmb_lf")) "lf()" else "nlf()"
    frm_stop(fn, " names resp = '", e2[["resp"]], "', but the bf() it is ",
             "added to models '", bform_resp_label(e1), "'", call. = FALSE)
  }
  if (inherits(e2, "frmtmb_lf")) {
    for (nm in names(e2$pforms)) {
      if (nm %in% c(names(e1$pforms), names(e1$pfix),
                    names(e1$nlforms))) {
        frm_stop("lf() sets '", nm, "', which the bf() it is added to ",
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
        frm_stop("nlf() sets '", nm, "', which the bf() it is added to ",
                 "already sets", call. = FALSE)
      }
      e1$nlforms[[nm]] <- e2$nlforms[[nm]]
    }
    for (nm in names(e2$pforms)) {
      if (set(nm)) {
        frm_stop("The linear parameter formulas passed to nlf() set '", nm,
                 "', which the bf() it is added to already sets",
                 call. = FALSE)
      }
      e1$pforms[[nm]] <- e2$pforms[[nm]]
    }
    return(e1)
  }
  if (inherits(e2, "frmtmb_rescor")) {
    frm_stop("set_rescor() applies to multivariate formulas; combine ",
             "responses with mvbf() or `bf() + bf()` first", call. = FALSE)
  }
  if (inherits(e2, "frmtmb_family") || inherits(e2, "family") ||
      is.function(e2)) {
    e1$family <- as_frmtmb_family(e2)
    return(e1)
  }
  frm_stop("Cannot add an object of class ", paste(class(e2), collapse = "/"),
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
#'   (gaussian only). `mvbf()` defaults to `FALSE`, `set_rescor()` to
#'   `TRUE`. It is brms's spelling on both, and brms's default on
#'   `set_rescor()`; brms's `mvbf()` defaults to `NULL` and decides
#'   later, which frmtmb settles at formula-assembly time instead.
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
  check_flag(rescor, "rescor")
  forms <- list(...)
  flat <- list()
  for (f in forms) {
    if (inherits(f, "frmtmb_mvformula")) {
      flat <- c(flat, f$forms)
      if (isTRUE(f$rescor)) rescor <- TRUE
    } else if (inherits(f, "frmtmb_formula")) {
      flat <- c(flat, list(f))
    } else {
      frm_stop("mvbf() takes bf() formulas", call. = FALSE)
    }
  }
  if (length(flat) < 2) {
    frm_stop("mvbf() needs at least two responses", call. = FALSE)
  }
  structure(list(forms = flat, rescor = isTRUE(rescor)),
            class = c("frmtmb_mvformula", "frmtmb_bform"))
}

#' @rdname mvbf
#' @param rescor_value The spelling `set_rescor()` shipped with, still
#'   accepted as an alias of `rescor`. It existed only because this Rd
#'   page documents two functions and could not carry two `rescor`
#'   entries; brms spells the argument `rescor`, and so does this
#'   function now. Give one spelling or the other, not both.
#' @export
set_rescor <- function(rescor = arg_unset(),
                       rescor_value = arg_unset()) {
  # both spellings default to the "not supplied" marker so that either
  # one alone is a setting and both together are a refusal; the value
  # the function acts on when neither is given is still TRUE
  rescor <- dual_arg(rescor, rescor_value, "rescor", "rescor_value",
                     "set_rescor()", default = TRUE)
  check_flag(rescor, "rescor",
             what = "(`rescor_value` is the same setting under its old name)")
  structure(list(rescor = isTRUE(rescor)), class = "frmtmb_rescor")
}

#' `mvbf() + e2`.
#'
#' @noRd
plus_mvbf <- function(e1, e2) {
  if (inherits(e2, "frmtmb_lf") || inherits(e2, "frmtmb_nlf")) {
    fn <- if (inherits(e2, "frmtmb_lf")) "lf()" else "nlf()"
    resps <- vapply(e1$forms, bform_resp_label, "")
    if (is.null(e2[["resp"]])) {
      eg <- if (inherits(e2, "frmtmb_lf")) "lf(sigma ~ z" else
        "nlf(sigma ~ a * z"
      frm_stop(fn, " does not say which response it belongs to. Name the ",
               "response, e.g. + ", eg, ", resp = \"", resps[length(resps)],
               "\"), or add it directly after the bf() it modifies. The ",
               "responses are: ", paste(resps, collapse = ", "),
               call. = FALSE)
    }
    at <- which(vapply(e1$forms, bform_resp_matches, TRUE,
                       resp = e2[["resp"]]))
    if (length(at) != 1L) {
      frm_stop(fn, " names resp = '", e2[["resp"]], "', which is not one of ",
               "the responses: ", paste(resps, collapse = ", "),
               call. = FALSE)
    }
    e1$forms[[at]] <- plus_bf(e1$forms[[at]], e2)
    return(e1)
  }
  if (inherits(e2, "frmtmb_rescor")) {
    e1$rescor <- e2$rescor
    return(e1)
  }
  if (inherits(e2, "frmtmb_mvformula") || inherits(e2, "frmtmb_formula")) {
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
  frm_stop("Cannot add an object of class ", paste(class(e2), collapse = "/"),
           " to a 'frmtmb_mvformula'", call. = FALSE)
}

#' @export
print.frmtmb_mvformula <- function(x, ...) {
  frm_check_dots(...)
  for (f in x$forms) print(f)
  cat("rescor:", x$rescor, "\n")
  invisible(x)
}

#' @export
print.frmtmb_formula <- function(x, ...) {
  frm_check_dots(...)
  cat(deparse1(x$formula), if (isTRUE(x$nl)) " (nonlinear)" else "", "\n",
      sep = "")
  for (f in x$nlforms) cat(deparse1(f), " (nonlinear)\n", sep = "")
  for (f in x$pforms) cat(deparse1(f), "\n")
  for (nm in names(x$pfix)) cat(nm, "=", x$pfix[[nm]], "\n")
  if (!is.null(x$family)) {
    cat("Family:", x$family[["family"]], "\n")
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
