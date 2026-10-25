#' Cox proportional hazards with a flexible baseline
#'
#' `cox()` is the flexible-parametric proportional hazards model brms
#' fits: the baseline hazard is an M-spline in time,
#' `h0(t) = sum_j s_j M_j(t)`, and the cumulative baseline hazard is the
#' I-spline integral of the same basis over the same weights, which
#' makes it monotone by construction. The weights `s` form a simplex,
#' which is what identifies the baseline against the intercept. They are
#' estimated as `sbhaz_raw`, their `Kbhaz - 1` free coordinates in the
#' chart a `mo()` simplex uses (see "Monotonic effects" in
#' [frmtmb::frm()]); [cox_baseline()] returns the simplex itself. The
#' default basis is brms's: `df = 5` cubic M-splines with an intercept,
#' internal knots on response quantiles and boundary knots just outside
#' the observed range.
#'
#' The hazard is `h0(t) exp(eta)`, so a coefficient is a log hazard
#' ratio, exactly as in `survival::coxph()`. Right, left, and interval
#' censoring come through the ordinary `cens()` addition term: an event
#' contributes the density and a censored observation the survivor
#' function, which is what this family's log-density and log-CDF are.
#' Random effects are the point: `time | cens(c) ~ x + (1 | g)` is a
#' frailty model, and the Laplace approximation integrates the frailties
#' out.
#'
#' The baseline is semiparametric only in spirit: it has `df`
#' parameters, not one per event time, so coefficients agree with
#' `coxph()` closely rather than exactly. A survival response has no
#' mean, so `fitted()` and `frm_linpred(type = "response")` are refused;
#' `frm_linpred(type = "link")` gives the log hazard ratio. `simulate()`
#' is not available.
#'
#' Maximum likelihood often puts one or more baseline weights ON the
#' simplex boundary, at exactly zero: an unpenalized flexible baseline
#' does that. A weight of zero is a finite point of the chart, which the
#' fit reaches with a finite Hessian, and a zero weight the likelihood
#' rises from is a saddle the fit is moved off
#' (`frmtmb_control(mo_search = )`). brms keeps the weights interior
#' with its `dirichlet(1)` prior on `sbhaz`; frmtmb.sample's
#' `frm_sample()` samples the baseline under the same prior, and
#' `set_prior("dirichlet(...)", class = "sbhaz")` puts another one on
#' it. Lower `df` until the baseline is one the data supports, and read
#' [cox_baseline()] to see which weights are at zero.
#'
#' @section Where it lives:
#' `cox()` was in frmtmb itself before this package existed. brms has it
#' in its core, so `family = "cox"` and `brmsfamily("cox")` resolve to
#' this function once this package is loaded, and refuse by name, naming
#' this package, when it is not.
#'
#' The baseline simplex reaches frmtmb through
#' `frmtmb_family(simplex_pars = c(sbhaz = "sbhaz_raw"))`, and its basis
#' through `frmtmb_family(aterm_data = )`. brms's `bhaz()` addition term
#' reaches it through `frmtmb_register_aterm(option = )`. Those seams
#' are all the family needs from core.
#'
#' @section brms's bhaz() addition term:
#' brms 2.23.0 sets the baseline basis on the response side:
#' `time | cens(c) + bhaz(df = 6) ~ x` with `family = cox()`. This
#' package takes that spelling as well. `bhaz(df =, degree =,
#' intercept =, knots =, Boundary.knots =)` sets the argument of the
#' same name of `cox()`, so these two models are the same model, with
#' the same basis to rounding as brms's `standata()`:
#'
#' ```
#' frm(bf(time | cens(c) + bhaz(df = 6) ~ x), family = cox(), data = d)
#' frm(bf(time | cens(c) ~ x), family = cox(df = 6), data = d)
#' ```
#'
#' The arguments of `bhaz()` are evaluated in the formula's
#' environment. An argument that both `bhaz()` and `cox()` set must
#' have the same value in both, or the fit stops and names the two.
#' `df` and `knots` must agree too, in either place: `length(knots) +
#' degree + intercept` basis functions. brms ignores `df` when `knots`
#' is given; frmtmb refuses the pair.
#'
#' Some differences from brms are refusals. `bhaz(gr = g)` gives each
#' level of `g` a baseline of its own (a stratified Cox model), and
#' this family has one baseline: fit each stratum as a model of its
#' own, or let `g` shift the one baseline by adding it to the formula.
#' `Boundary.knots` must cover every observed time. brms takes boundary
#' knots inside the data and extrapolates the basis there, where the
#' cumulative baseline hazard can be negative, so that the survivor
#' function exceeds one. An argument name that brms drops without a
#' message, such as `bhaz(int = FALSE)` or `bhaz(Bound = c(0, 30))`,
#' is refused with the name that was probably meant. A shortened name
#' that brms reads as an argument (`bhaz(d = 6)` for `df`, `kn =` for
#' `knots`, `deg =` for `degree`) is read the same way here.
#' `warn.outside` is taken and has no effect, as in brms.
#'
#' brms 2.23.0's `cox()` has no `bhaz =` argument (brms before 2.22.0
#' took the basis that way). `cox(bhaz = )` stops and names the two
#' spellings above.
#'
#' @param link The link for the hazard ratio `mu`: `"log"` (the
#'   default), `"identity"`, `"softplus"` or `"squareplus"`, the set
#'   brms 2.23.0 allows. Quoted or not.
#' @param df The number of M-spline basis functions in the baseline
#'   hazard (brms's `bhaz(df = 5)` default). With `knots`, it is
#'   `length(knots) + degree + intercept`, and need not be given.
#' @param degree The spline degree of that basis (cubic by default).
#' @param intercept Keep the basis function that is non-zero at the
#'   lower boundary knot.
#' @param knots The internal knots of the basis. The default, `NULL`,
#'   puts `df - degree - intercept` of them on equally spaced quantiles
#'   of the observed times, as brms and `splines2::mSpline()` do. Given,
#'   they must be distinct and lie strictly inside the boundary knots.
#' @param Boundary.knots The two boundary knots. The default, `NULL`,
#'   is brms's: the range of the observed times widened by 1/50 of its
#'   width at each end, with the lower knot not below zero. Given, they
#'   must cover every observed time.
#' @param bhaz Not an argument of brms 2.23.0's `cox()`. Given, it is
#'   refused with the spellings to use instead.
#' @return A `frmtmb_family` object, ready for [frmtmb::frm()].
#' @seealso [cox_baseline()] for the fitted baseline weights,
#'   [royston_parmar()] for a flexible parametric model on the log
#'   cumulative hazard.
#' @examples
#' set.seed(2)
#' n <- 200
#' dd <- data.frame(x = rnorm(n))
#' dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
#' dd$out <- rbinom(n, 1, 0.3)        # 1 = right censored
#' cox_fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
#' fixef_by_dpar(cox_fit)$mu          # log hazard ratios
#' cox_baseline(cox_fit)              # the baseline hazard weights
#' @export
cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE,
                knots = NULL, Boundary.knots = NULL, bhaz = NULL) {
  if (!missing(bhaz)) {
    frm_stop("cox(bhaz = ) is the spelling of brms before 2.22.0, and ",
             "brms 2.23.0's cox() has no such argument. Set the basis ",
             "in cox() itself, cox(df = 6), or with brms's addition ",
             "term, time | bhaz(df = 6) ~ x", call. = FALSE)
  }
  link <- sv_link_value(substitute(link), link, cox_links, "log")
  link <- sv_check_link(link, cox_links, "cox")
  # only what the call writes, so that bhaz() can tell an argument set
  # twice from one left at its default
  args <- list(df = df, degree = degree, intercept = intercept,
               knots = knots, Boundary.knots = Boundary.knots)
  given <- c(!missing(df), !missing(degree), !missing(intercept),
             !is.null(knots), !is.null(Boundary.knots))
  # a df too small for its degree is refused at the fit, not here: a
  # bhaz(intercept = FALSE) or bhaz(degree =) in the formula can still
  # make the merged basis valid (cox_bhaz_option() checks it then)
  cox_build(link, args[given],
            stats::setNames(rep("cox", sum(given)), names(args)[given]),
            check_df = FALSE)
}

#' A `cox()` family from its link and the basis arguments that were
#' written, each with the place it was written in (`"cox"` or
#' `"bhaz"`), which the refusals name. The arguments that were not
#' written take `cox()`'s defaults. The family records what `cox()`
#' itself was given, so that `bhaz()` applied to it again gives the
#' same family (`frmtmb_register_aterm(option =)`). `check_df = FALSE`
#' leaves a `df` below `degree + intercept` to the frame
#' (`bhaz_spec()`).
#'
#' @noRd
cox_build <- function(link, args, where, check_df = TRUE) {
  lab <- function(nm) paste0(where[[nm]], "(", nm, " =)")
  for (nm in intersect(c("df", "degree"), names(args))) {
    # df and degree size the I-spline baseline basis, so a length-2
    # value used to build a basis of one shape and record another
    sv_check_count(args[[nm]], if (where[[nm]] == "cox") nm else lab(nm),
                   min = 1L)
  }
  if ("intercept" %in% names(args)) {
    sv_check_flag(args$intercept,
                  if (where[["intercept"]] == "cox") "intercept" else
                    lab("intercept"))
  }
  df <- args$df %||% 5
  degree <- args$degree %||% 3
  intercept <- args$intercept %||% TRUE
  knots <- args$knots
  boundary <- args$Boundary.knots
  if (!is.null(boundary)) {
    ok <- is.numeric(boundary) && length(boundary) == 2L &&
      all(is.finite(boundary)) && boundary[1L] < boundary[2L]
    if (!ok) {
      frm_stop(lab("Boundary.knots"), " takes two finite numbers, the ",
               "lower boundary knot below the upper one, not ",
               sv_arg_desc(boundary), call. = FALSE)
    }
    boundary <- as.numeric(boundary)
  }
  if (!is.null(knots)) {
    if (!is.numeric(knots) || !all(is.finite(knots)) ||
          anyDuplicated(knots)) {
      frm_stop(lab("knots"), " takes distinct finite numbers, the ",
               "internal knots of the baseline basis, not ",
               sv_arg_desc(knots), call. = FALSE)
    }
    knots <- sort(as.numeric(knots))
    if (!is.null(boundary) &&
          any(knots <= boundary[1L] | knots >= boundary[2L])) {
      frm_stop("The internal knots of the baseline basis must lie ",
               "strictly inside its boundary knots: ", lab("knots"),
               " gives ", cox_num_text(knots), ", and ",
               lab("Boundary.knots"), " gives ", cox_num_text(boundary),
               call. = FALSE)
    }
    n_basis <- length(knots) + as.integer(degree) + as.integer(intercept)
    if ("df" %in% names(args) && df != n_basis) {
      frm_stop(lab("df"), " is ", df, ", and ", lab("knots"), " gives ",
               length(knots), " internal knots, which with degree ",
               degree, if (intercept) " and an intercept" else
                 " and no intercept", " make ", n_basis, " basis ",
               "functions. Give df or knots, not both (brms ignores df ",
               "when knots is given)", call. = FALSE)
    }
    df <- n_basis
  } else if (check_df && df < degree + intercept) {
    from <- if (any(where[intersect(c("df", "degree", "intercept"),
                                    names(where))] == "bhaz")) {
      "bhaz()"
    } else {
      "cox()"
    }
    cox_stop_df(from, df, degree, intercept)
  }
  fam <- fam_cox(link, df = df, degree = degree, intercept = intercept,
                 knots = knots, boundary = boundary)
  fam[["cox_args"]] <- list(link = link,
                            args = args[names(where)[where == "cox"]])
  fam
}

#' brms's `bhaz()` addition term, as the option of `cox()` it is
#' (`frmtmb_register_aterm(option =)`): `bhaz(gr = NA, df = 5, ...)` in
#' brms 2.23.0 (`resp_bhaz()`), whose `...` reaches
#' `splines2::mSpline()`. The arguments are named as brms names them
#' (`cox_bhaz_args()`), so an unnamed first argument is `gr`.
#'
#' @noRd
cox_bhaz_option <- function(family, call, env) {
  if (!identical(family[["family"]], "cox") ||
        is.null(family[["cox_args"]])) {
    frm_stop("Addition term `bhaz()` sets the baseline hazard of a cox() ",
             "model, and family '", family[["family"]], "' has none. ",
             "brms takes bhaz() for cox() alone too", call. = FALSE)
  }
  a <- cox_bhaz_args(call)
  gr <- a[["gr"]]
  if (!is.null(gr) && !(is.logical(gr) && length(gr) == 1L && is.na(gr))) {
    if (is.atomic(gr)) {
      # brms reads bhaz(6) as ONE stratum named "6", with the default
      # df, which is almost never what was meant
      frm_stop("bhaz(", deparse1(gr), "): the first argument of bhaz() is ",
               "gr, the grouping factor of a stratified baseline, not df. ",
               "brms reads it as one stratum named ", deparse1(gr),
               " and keeps df = 5. Write bhaz(df = ", deparse1(gr), ")",
               call. = FALSE)
    }
    frm_stop("bhaz(gr = ", deparse1(gr), ") gives each level of ",
             deparse1(gr), " a baseline hazard of its own (a stratified ",
             "Cox model), and frmtmb's cox() has one baseline. Fit each ",
             "stratum as a model of its own, or add ", deparse1(gr),
             " to the formula, which shifts the one baseline by a hazard ",
             "ratio per level", call. = FALSE)
  }
  a[["gr"]] <- NULL
  # splines2 reads warn.outside only for its warning, and frmtmb's
  # refusal of data outside the boundary knots replaces that warning
  a[["warn.outside"]] <- NULL
  vals <- lapply(names(a), function(nm) {
    tryCatch(eval(a[[nm]], env), error = function(e) {
      frm_stop("bhaz(): cannot evaluate ", nm, " = ", deparse1(a[[nm]]),
               ": ", conditionMessage(e), call. = FALSE)
    })
  })
  names(vals) <- names(a)
  # a NULL argument is one not written, as cox() reads it
  vals <- vals[!vapply(vals, is.null, NA)]
  prior <- family[["cox_args"]]
  given <- prior$args
  for (nm in intersect(names(vals), names(given))) {
    if (!cox_same_arg(vals[[nm]], given[[nm]])) {
      frm_stop("bhaz(", nm, " = ", cox_arg_text(vals[[nm]]), ") and cox(",
               nm, " = ", cox_arg_text(given[[nm]]), ") set the baseline ",
               "basis twice and disagree. Set ", nm, " in one place",
               call. = FALSE)
    }
  }
  args <- c(given[setdiff(names(given), names(vals))], vals)
  where <- stats::setNames(
    ifelse(names(args) %in% names(vals), "bhaz", "cox"), names(args))
  out <- cox_build(prior$link, args, where)
  # bhaz() applied again starts from what cox() itself was given
  out[["cox_args"]] <- prior
  out
}

#' The arguments of a `bhaz()` call, named as brms 2.23.0 names them,
#' as unevaluated expressions. Not `match.call()`: R's partial matching
#' there takes `int =` for `intercept`, while brms never does
#' (`dev/bhaz-brms-partial.R`). brms matches in two steps:
#'
#' 1. `resp_bhaz(gr = NA, df = 5, ...)`: exact names, then a unique
#'    prefix of `gr` or `df` (`d = 6` is `df`), then the first two
#'    unnamed arguments as `gr` and `df`. A third unnamed argument is
#'    passed on and ignored.
#' 2. Everything else goes to `splines2::mSpline()`, after brms has
#'    set `intercept` and `Boundary.knots` by their exact names. So a
#'    prefix of `knots`, `degree` or `warn.outside` is taken (`kn`,
#'    `deg`), and a prefix of `intercept` or `Boundary.knots` (`int`,
#'    `inter`, `Bound`) is dropped without a word, or reaches
#'    `integral`.
#'
#' What brms takes is taken here. What brms drops or misreads is
#' refused, with the argument that was probably meant: a setting that
#' is dropped silently fits a different model from the one written.
#'
#' @noRd
cox_bhaz_args <- function(call) {
  a <- as.list(call)[-1L]
  nms <- names(a) %||% rep("", length(a))
  out <- list()
  put <- function(full, i) {
    if (!is.null(out[[full]])) {
      # brms stops on a full name given twice; a shortened one beside
      # the full name never reaches the argument, and brms drops it
      frm_stop("bhaz() sets ", full, " twice, as `", nms[[i]], " =` ",
               "and earlier in the same call. ",
               if (identical(nms[[i]], full)) {
                 "brms stops on this too"
               } else {
                 paste0("brms keeps `", full, " =` and drops `", nms[[i]],
                        " =` without a message")
               }, call. = FALSE)
    }
    out[[full]] <<- a[[i]]
  }
  taken <- c("gr", "df", "degree", "intercept", "knots", "Boundary.knots",
             "warn.outside")
  for (i in which(nzchar(nms) & nms %in% taken)) put(nms[[i]], i)
  for (i in which(nzchar(nms) & !nms %in% taken)) {
    nm <- nms[[i]]
    pre <- function(cands) cands[startsWith(cands, nm)]
    step1 <- pre(setdiff(c("gr", "df"), names(out)))
    if (length(step1) == 1L) {
      put(step1, i)
      next
    }
    step2 <- pre(c("knots", "degree", "periodic", "derivs", "integral",
                   "warn.outside"))
    meant <- pre(c("intercept", "Boundary.knots"))
    if (length(step2) == 1L && step2 %in% c("knots", "degree",
                                             "warn.outside")) {
      put(step2, i)
      next
    }
    if (length(meant)) {
      frm_stop("bhaz(", nm, " = ) is not an argument of bhaz(): write ",
               meant[[1L]], " =. brms ",
               if (length(step2)) {
                 paste0("passes `", nm, " =` to splines2's ", step2[[1L]])
               } else {
                 "drops it without a message"
               }, " and keeps ", meant[[1L]], " at its default, a ",
               "different model from the one written", call. = FALSE)
    }
    if (length(step2) == 1L || nm %in% c("periodic", "derivs",
                                         "integral")) {
      full <- if (length(step2) == 1L) step2 else nm
      frm_stop("bhaz(", nm, " = ) sets splines2's ", full, ", which ",
               "changes the baseline basis in brms. cox() takes df, ",
               "degree, intercept, knots and Boundary.knots only",
               call. = FALSE)
    }
    frm_stop("bhaz() takes brms's arguments for the baseline basis, df, ",
             "degree, intercept, knots and Boundary.knots",
             if (length(step2) > 1L) {
               paste0(", and `", nm, " =` could be any of ",
                      paste(step2, collapse = " or "), " (brms stops on ",
                      "it too)")
             } else {
               paste0(". `", nm, " =` is none of them; brms drops it ",
                      "without a message")
             }, call. = FALSE)
  }
  pos <- which(!nzchar(nms))
  free <- setdiff(c("gr", "df"), names(out))
  if (length(pos) > length(free)) {
    frm_stop("bhaz() takes at most two unnamed arguments, gr and df, as ",
             "brms does; brms ignores a third. Name the argument: ",
             "bhaz(knots = ...)", call. = FALSE)
  }
  for (j in seq_along(pos)) out[[free[[j]]]] <- a[[pos[[j]]]]
  out
}

#' Whether two values of one basis argument are the same: knots as a
#' set, numbers by value whatever their storage mode.
#'
#' @noRd
cox_same_arg <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    return(identical(sort(as.numeric(a)), sort(as.numeric(b))))
  }
  identical(a, b)
}

#' One basis argument as a refusal prints it.
#'
#' @noRd
cox_arg_text <- function(x) {
  if (is.atomic(x) && length(x) == 1L) return(format(x))
  deparse1(x)
}

#' Refuse a basis with fewer functions than its degree needs; `from`
#' is the place the refusal names (`"cox()"` or `"bhaz()"`).
#'
#' @noRd
cox_stop_df <- function(from, df, degree, intercept) {
  frm_stop(from, ": df must be at least degree", if (intercept) " + 1",
           " (", degree + intercept, ") with",
           if (intercept) " an" else "out an",
           " intercept in the baseline basis, and it is ", df,
           call. = FALSE)
}

#' Numbers as a refusal prints them: each to 7 significant digits, with
#' no padding to a common width.
#'
#' @noRd
cox_num_text <- function(x) {
  paste(vapply(x, function(v) format(v, digits = 7L), ""), collapse = ", ")
}

#' The links brms 2.23.0 allows for `cox()`'s mean, its default first:
#' `brms:::.family_cox()$links`. `test-cox-links.R` reads the installed
#' brms and fails if this has drifted.
#'
#' @noRd
cox_links <- c("log", "identity", "softplus", "squareplus")

#' The fitted baseline-hazard simplex of a `cox()` fit
#'
#' @param fit A `frmtmb_fit` with a [cox()] family.
#' @return The `Kbhaz` M-spline weights, summing to one.
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = rnorm(200))
#' dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
#' fit <- frm(bf(time ~ x), family = cox(), data = dd)
#' cox_baseline(fit)
#' @export
cox_baseline <- function(fit) {
  rspec <- single_response(fit, "cox_baseline()")
  fam <- rspec$family
  if (is.null(fam[["cox_sbhaz"]])) {
    frm_stop("cox_baseline() needs a cox() family fit", call. = FALSE)
  }
  s <- fam[["cox_sbhaz"]](fit$estimates[["sbhaz_raw"]])
  stats::setNames(s, paste0("s", seq_along(s)))
}


#' B-spline basis of a given order (degree + 1) over a full knot vector,
#' by the Cox-de Boor recursion. Written out rather than taken from
#' `splines::splineDesign()` so the baseline hazard needs no dependency
#' outside the ones the package already carries; it agrees with
#' `splines2` to machine precision, which the tests assert.
#'
#' @noRd
bspline_basis <- function(x, knots, order) {
  n <- length(x)
  nb1 <- length(knots) - 1L
  # order 1: the indicator of [t_j, t_{j+1}); the last non-degenerate
  # interval closes on the right so the upper boundary is covered
  B <- matrix(0, n, nb1)
  last <- max(which(knots < knots[length(knots)]))
  for (j in seq_len(nb1)) {
    B[, j] <- if (j == last) {
      as.numeric(x >= knots[j] & x <= knots[j + 1L])
    } else {
      as.numeric(x >= knots[j] & x < knots[j + 1L])
    }
  }
  if (order == 1L) return(B[, seq_len(length(knots) - order), drop = FALSE])
  for (k in 2:order) {
    Bn <- matrix(0, n, length(knots) - k)
    for (j in seq_len(ncol(Bn))) {
      d1 <- knots[j + k - 1L] - knots[j]
      d2 <- knots[j + k] - knots[j + 1L]
      # a repeated knot gives a zero-width span; its term drops out
      t1 <- if (d1 > 0) (x - knots[j]) / d1 * B[, j] else 0
      t2 <- if (d2 > 0) (knots[j + k] - x) / d2 * B[, j + 1L] else 0
      Bn[, j] <- t1 + t2
    }
    B <- Bn
  }
  B
}

#' The knot placement brms uses for the Cox baseline: internal knots on
#' equally spaced response quantiles, boundary knots just outside the
#' observed range so every observation sits strictly inside the basis
#' (`splines2`'s default through brms's `bhaz_basis_matrix()`). Knots
#' and boundary knots the user gave (`cox(knots =, Boundary.knots =)`,
#' or brms's `bhaz()`) replace the ones placed here.
#'
#' @noRd
bhaz_spec <- function(y, df, degree, intercept, knots = NULL,
                      boundary = NULL) {
  n_int <- df - degree - as.integer(intercept)
  # cox(df =) alone reaches here unchecked (cox_build(check_df = FALSE))
  if (n_int < 0L) cox_stop_df("cox()", df, degree, intercept)
  rng <- range(y)
  d <- rng[2L] - rng[1L]
  if (is.null(boundary)) {
    boundary <- c(max(rng[1L] - d / 50, 0), rng[2L] + d / 50)
  } else if (rng[1L] < boundary[1L] || rng[2L] > boundary[2L]) {
    # brms extrapolates the basis past its boundary knots. With
    # Boundary.knots = c(0.5, 10) on times from 0.0044 to 14.42, 62 rows
    # lie below the lower knot and 1 above; the cumulative baseline
    # hazard at equal weights reaches -1.082 there, so S(t) > 1, and the
    # hazard is negative for some weights (no row is negative for every
    # weight). dev/review-bhaz-bnd.R, dev/bhaz-brms-probe2.R
    frm_stop("cox(): the boundary knots of the baseline basis, ",
             cox_num_text(boundary), ", must cover every observed time, ",
             "and the times run from ", cox_num_text(rng[1L]), " to ",
             cox_num_text(rng[2L]), ". Outside its ",
             "boundary knots the basis is an extrapolation, where the ",
             "cumulative baseline hazard can be negative. Widen ",
             "Boundary.knots, or leave them out for brms's default",
             call. = FALSE)
  }
  if (!is.null(knots)) {
    if (any(knots <= boundary[1L] | knots >= boundary[2L])) {
      frm_stop("cox(): the internal knots of the baseline basis, ",
               cox_num_text(knots), ", must lie strictly inside its ",
               "boundary knots, ", cox_num_text(boundary),
               ". The default boundary knots are the range of the ",
               "observed times widened by 1/50 of its width at each end",
               call. = FALSE)
    }
    return(list(internal = knots, boundary = boundary, degree = degree,
                intercept = intercept, df = df))
  }
  internal <- if (n_int > 0L) {
    p <- seq(0, 1, length.out = n_int + 2L)[-c(1L, n_int + 2L)]
    unname(stats::quantile(y, probs = p))
  } else {
    numeric(0)
  }
  list(internal = internal, boundary = boundary, degree = degree,
       intercept = intercept, df = df)
}

#' The full knot vector of a baseline-hazard spec: boundary knots
#' repeated to the spline's order at each end.
#'
#' @noRd
bhaz_knots <- function(sp, extra = 0L) {
  k <- sp$degree + 1L + extra
  c(rep(sp$boundary[1L], k), sp$internal, rep(sp$boundary[2L], k))
}

#' The M-spline basis of the baseline hazard: the B-spline basis scaled
#' so each function integrates to one over its support, which is what
#' makes a simplex of coefficients a density-like hazard.
#'
#' @noRd
mspline_design <- function(x, sp) {
  k <- sp$degree + 1L
  tt <- bhaz_knots(sp)
  N <- bspline_basis(x, tt, k)
  M <- N
  for (j in seq_len(ncol(N))) M[, j] <- N[, j] * k / (tt[j + k] - tt[j])
  if (!sp$intercept) M <- M[, -1L, drop = FALSE]
  M
}

#' The I-spline basis: `I_j(x)` is the integral of `M_j` from the lower
#' boundary knot, so it is the CUMULATIVE baseline hazard's basis and it
#' is monotone by construction. It equals the reverse cumulative sum of
#' the order-(degree + 2) B-spline basis on the once-more-repeated
#' boundary knots (Ramsay 1988); the tests check it against both
#' `splines2::iSpline()` and a numeric integral of `mspline_design()`.
#'
#' @noRd
ispline_design <- function(x, sp) {
  N2 <- bspline_basis(x, bhaz_knots(sp, extra = 1L), sp$degree + 2L)
  R <- t(apply(N2, 1L, function(v) rev(cumsum(rev(v)))))
  if (length(x) == 1L) R <- matrix(R, 1L, ncol(N2))
  I <- R[, -1L, drop = FALSE]
  if (!sp$intercept) I <- I[, -1L, drop = FALSE]
  I
}

#' Cox proportional hazards with a flexible (M-spline) baseline hazard,
#' brms's construction. The baseline hazard is `Zbhaz %*% s` and the
#' cumulative baseline hazard `Zcbhaz %*% s` over the SAME simplex `s`,
#' with `Zcbhaz` the I-spline integral of the M-spline `Zbhaz`. The
#' simplex is what identifies the baseline against the intercept, and it
#' rides in the parameter vector as `K - 1` free coordinates of the
#' chart `mo_simplex()` reads a `mo()` simplex in.
#'
#' @noRd
fam_cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE,
                    knots = NULL, boundary = NULL) {
  df <- as.integer(df)
  degree <- as.integer(degree)
  intercept <- isTRUE(intercept)
  # The simplex from its K-1 free coordinates, in the chart a mo()
  # simplex uses: a weight of 0 is at a finite point there, where the
  # softmax this used put it at an infinite coordinate and fits stopped
  # on the plateau toward it (weights of 1.1e-8, 8.3e-9 and 5.9e-8 on
  # seeds 2, 5 and 6 of dev/optima-rev-cox1.R; dev/moparity-findings.md,
  # item 4). `chart` is the response's `sbhaz_chart` addition term, which
  # frmtmb.sample's sampler sets to "softmax" on its own copy of the
  # frame (frmtmb::mo_chart_frame(), for simplex_pars below).
  sbhaz <- function(raw, chart = NULL) mo_simplex(raw, chart)
  fam <- frmtmb_family(
    "cox",
    # bhaz() is brms's spelling of this family's basis arguments
    # (cox_bhaz_option()); it carries no data
    accepts_aterms = c("weights", "cens", "trunc", "bhaz"),
    dpars = "mu",
    links = list(mu = link),
    lpdf = function(y, dpars, aterms, extra) {
      s <- sbhaz(extra$sbhaz_raw, aterms[["sbhaz_chart"]])
      bhaz <- as.vector(aterms[["Zbhaz"]] %*% s)
      cbhaz <- as.vector(aterms[["Zcbhaz"]] %*% s)
      # log h(t) + log S(t), with h(t) = h0(t) * mu and mu = exp(eta)
      log(bhaz) + log(dpars[["mu"]]) - cbhaz * dpars[["mu"]]
    },
    lcdf = function(q, dpars, aterms, extra) {
      # F(t) = 1 - exp(-H0(t) * mu). Every bound a censored or truncated
      # row asks about is plain data, so its I-spline row block was
      # built once at frame assembly and is looked up by value here.
      cbhaz <- as.vector(cox_cbhaz_design(q, aterms) %*%
                           sbhaz(extra$sbhaz_raw, aterms[["sbhaz_chart"]]))
      1 - exp(-cbhaz * dpars[["mu"]])
    },
    lccdf = function(q, dpars, aterms, extra) {
      # log S = -H0(t) * mu, the cumulative hazard itself
      cbhaz <- as.vector(cox_cbhaz_design(q, aterms) %*%
                           sbhaz(extra$sbhaz_raw, aterms[["sbhaz_chart"]]))
      -cbhaz * dpars[["mu"]]
    },
    valid_y = function(y, aterms) {
      if (any(y <= 0)) {
        frm_stop("cox: the response is a survival time and must be ",
                 "strictly positive", call. = FALSE)
      }
      if (length(unique(y)) < df) {
        frm_stop("cox: fewer distinct event times (", length(unique(y)),
                 ") than baseline basis functions (", df,
                 "); lower df in cox(df = )", call. = FALSE)
      }
    },
    init_dpars = list(
      # the baseline integrates to one over the observed window, so the
      # hazard ratio starts at the crude event rate
      mu = function(y, aterms) {
        ev <- if (is.null(aterms[["cens"]])) 1 else mean(aterms[["cens"]] == 0)
        max(ev, 0.05) / mean(y)
      }
    ),
    type = "continuous",
    post = list(
      # mu is the hazard RATIO, not a mean, and reporting it as one is
      # exactly the mistake this family invites. brms refuses the same
      # question (posterior_epred has no cox method).
      mean_fn = function(dpars, aterms) {
        frm_stop("cox: a survival time has no mean on the response scale ",
                 "here - the model fits a hazard, and the mean survival ",
                 "time would be an integral over the baseline that the ",
                 "censored rows do not identify. frm_linpred(type = \"link\") ",
                 "gives the log hazard ratio and cox_baseline() the fitted ",
                 "baseline weights", call. = FALSE)
      }
    ),
    extra_pars = function(y, aterms) {
      list(sbhaz_raw = rep(0, ncol(aterms[["Zbhaz"]]) - 1L))
    },
    # a deliberate omission, not a gap: every entry point repeats this
    # reason after its own refusal
    sim_refusal = paste0(
      "Drawing a survival time means inverting the cumulative baseline ",
      "hazard, and cox() carries no quantile function for it: the ",
      "I-spline baseline is identified only on the observed time ",
      "window, so a draw beyond the last event time has no defined ",
      "distribution"
    ),
    # the baseline bases are data, not parameters: they are built once
    # from the observed times and ride with the response's addition
    # terms
    aterm_data = function(y, aterms) {
      sp <- bhaz_spec(y, df, degree, intercept, knots, boundary)
      out <- list(Zbhaz = mspline_design(y, sp), bhaz_spec = sp,
                  cox_cb = list(cox_cb_entry(y, sp)))
      for (nm in c("cens_y2", "trunc_lb", "trunc_ub")) {
        if (!is.null(aterms[[nm]])) {
          out$cox_cb[[length(out$cox_cb) + 1L]] <-
            cox_cb_entry(aterms[[nm]], sp)
        }
      }
      out$Zcbhaz <- out$cox_cb[[1L]]$Z
      out
    },
    # the baseline weights are a simplex, which core searches, charts
    # and priors as brms's class sbhaz
    simplex_pars = c(sbhaz = "sbhaz_raw")
  )
  fam[["cox_sbhaz"]] <- sbhaz
  fam
}

#' One cached (bound, I-spline design) pair. A bound outside the
#' boundary knots is clamped first: below the lower knot the cumulative
#' baseline hazard is zero and above the upper one it is its total, and
#' the raw B-spline recursion returns an all-zero row for both, which
#' would read as "no risk accumulated" at an upper bound.
#'
#' @noRd
cox_cb_entry <- function(q, sp) {
  qq <- pmin(pmax(as.numeric(q), sp$boundary[1L]), sp$boundary[2L])
  list(q = as.numeric(q), Z = ispline_design(qq, sp))
}

#' The cumulative-baseline-hazard design for a bound the Cox likelihood
#' asks about. Frame assembly cached one per bound; anything else (a
#' post-fit query) is built on the spot.
#'
#' @noRd
cox_cbhaz_design <- function(q, aterms) {
  qn <- as.numeric(q)
  for (e in aterms[["cox_cb"]] %||% list()) {
    if (identical(e$q, qn)) return(e$Z)
  }
  cox_cb_entry(qn, aterms[["bhaz_spec"]])$Z
}

