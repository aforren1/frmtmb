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
#' through `frmtmb_family(aterm_data = )`. Those two seams are all the
#' family needs from core.
#'
#' @param link The link for the hazard ratio `mu`: `"log"` (the
#'   default), `"identity"`, `"softplus"` or `"squareplus"`, the set
#'   brms 2.23.0 allows. Quoted or not.
#' @param df The number of M-spline basis functions in the baseline
#'   hazard (brms's `bhaz(df = 5)` default).
#' @param degree The spline degree of that basis (cubic by default).
#' @param intercept Keep the basis function that is non-zero at the
#'   lower boundary knot.
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
cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE) {
  link <- sv_link_value(substitute(link), link, cox_links, "log")
  link <- sv_check_link(link, cox_links, "cox")
  # df and degree size the I-spline baseline basis, so a length-2 value
  # used to build a basis of one shape and record another
  sv_check_count(df, "df", min = 1L)
  sv_check_count(degree, "degree", min = 1L)
  sv_check_flag(intercept, "intercept")
  fam_cox(link, df = df, degree = degree, intercept = intercept)
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
#' (`splines2`'s default through brms's `bhaz_basis_matrix()`).
#'
#' @noRd
bhaz_spec <- function(y, df, degree, intercept) {
  n_int <- df - degree - as.integer(intercept)
  if (n_int < 0L) {
    frm_stop("cox(): df must be at least degree + 1 (", degree + 1L,
             ") with an intercept in the baseline basis", call. = FALSE)
  }
  rng <- range(y)
  d <- rng[2L] - rng[1L]
  boundary <- c(max(rng[1L] - d / 50, 0), rng[2L] + d / 50)
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
fam_cox <- function(link = "log", df = 5, degree = 3, intercept = TRUE) {
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
    accepts_aterms = c("weights", "cens", "trunc"),
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
      sp <- bhaz_spec(y, df, degree, intercept)
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

