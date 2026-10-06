# Standard errors the outer Hessian cannot give: found at fit time,
# recovered where the Hessian still carries them, and said once.
#
# sdreport() inverts the outer Hessian with solve(), which refuses a
# matrix whose reciprocal condition number is below machine epsilon and
# then fills EVERY standard error with NaN. On `ls ~ mo(income) * age`
# (brms_monotonic's own data code, dev/nanse-mo-sweep.R) that happened
# on 72 of 200 data sets, 55 of them with optimizer code 0 and no
# warning: the interaction's simplex sits with a weight at 0, its
# softmax coordinate has run to -20 or beyond, and the likelihood is
# flat along that one coordinate to 1e-20. The rest of the matrix is
# well conditioned. rcond() separated the two outcomes exactly: at most
# 1.54e-16 on every all-NaN fit, at least 2.68e-16 on every finite one.

#' The outer Hessian exactly as `autoscale_sdreport()` builds it, at the
#' point it builds it, with the objective's state left as it was found.
#'
#' `optimHess()` evaluates the objective at points around the optimum,
#' and TMB keeps the best point it has evaluated, which sdreport() then
#' reads; so the state is saved and restored, as nl_flat_message() does.
#'
#' Returns `list(H, at, Hq, E)`: `H` in natural units, `at` the outer
#' parameter vector it was taken at, `Hq` the autoscaled Hessian when
#' `par_units` are in play (what `autoscale_sdreport()` inverts), and
#' `E` the noise of `H` (fd_hessian()).
#'
#' @noRd
fit_outer_hessian <- function(fit) {
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  u <- fit$par_units
  if (is.null(u) || all(u == 1)) {
    at <- sdr_outer_point(obj)
    h <- fd_hessian(at, obj$gr)
    return(list(H = h$H, at = at, Hq = NULL, E = h$E))
  }
  q0 <- fit$opt$par / u
  h <- fd_hessian(q0, function(q) obj$gr(q * u) * u)
  list(H = h$H / outer(u, u), at = fit$opt$par, Hq = h$H,
       E = h$E / outer(u, u))
}

#' stats::optimHess() with its default controls, written out so that the
#' unsymmetrized differences are kept.
#'
#' The arithmetic is optimHess()'s C code step for step: the same
#' gradient calls in the same order, the base point moved by `+ eps`,
#' `- 2 * eps` and `+ eps` in place, column `i` as
#' `(g(+) - g(-)) / (2 * eps)`, then `0.5 * (H[i, j] + H[j, i])`. So `H`
#' is the matrix sdreport() would build, bit for bit
#' (test-se-check.R pins it with identical()), and sdreport() can be
#' handed it.
#'
#' `E` is half the asymmetry of the differences before they are
#' symmetrized: the Hessian's own noise, which on a Laplace objective
#' comes from the inner solve. Whether to build it at fit time at all is
#' se_check_at_fit()'s decision, made before the first gradient.
#'
#' @noRd
fd_hessian <- function(par, gr) {
  n <- length(par)
  eps <- 0.001
  dpar <- par
  H <- matrix(0, n, n)
  for (i in seq_len(n)) {
    dpar[i] <- dpar[i] + eps
    df1 <- as.numeric(gr(dpar))
    dpar[i] <- dpar[i] - 2 * eps
    df2 <- as.numeric(gr(dpar))
    H[, i] <- (1 * (df1 - df2)) / (2 * eps * 1 * 1)
    dpar[i] <- dpar[i] + eps
  }
  E <- abs(H - t(H)) / 2
  H <- 0.5 * (H + t(H))
  # optimHess() names both margins after the parameters, and sdreport()
  # names its covariance from the matrix it is handed
  dimnames(H) <- list(names(par), names(par))
  list(H = H, E = E)
}

#' The parts of a TMB objective's environment that evaluating it moves:
#' the best point seen, which sdreport() reads, and the last point,
#' which seeds the next inner solve.
#'
#' @noRd
obj_state_save <- function(obj) {
  env <- obj$env
  keep <- intersect(c("last.par", "last.par.best", "value.best",
                      "last.par.ok", "last.par1", "last.par2"),
                    ls(env, all.names = TRUE))
  mget(keep, envir = env)
}

#' @noRd
obj_state_restore <- function(obj, saved) {
  for (nm in names(saved)) assign(nm, saved[[nm]], envir = obj$env)
  invisible(NULL)
}

#' The outer point sdreport() reads when it is not given one.
#'
#' @noRd
sdr_outer_point <- function(obj) {
  par <- obj$env$last.par.best
  r <- obj$env$random
  if (length(r)) par[-r] else par
}

#' A covariance whose every variance is a positive number, which is when
#' every standard error is finite. A Hessian that is not positive
#' definite can still give one; that case keeps sdreport()'s own
#' verdict (`pdHess`) and is not repaired here.
#'
#' @noRd
cov_usable <- function(V) {
  !inherits(V, "try-error") && is.matrix(V) && all(is.finite(V)) &&
    all(diag(V) > 0)
}

#' The covariance of the outer parameters from their Hessian, and which
#' parameters have no standard error, with the reason for each.
#'
#' Three tiers, each tried only when the one before has failed, so a fit
#' whose plain inverse is usable gets exactly what sdreport() gives:
#'
#' 1. `solve(H)`, sdreport()'s own inverse.
#' 2. The same inverse computed on `H` scaled to unit diagonal. It is
#'    the same matrix, so this can only differ from tier 1 where tier 1
#'    refused. For a fit without random effects `H` is first replaced by
#'    the exact AD Hessian `obj$he()`, because optimHess()'s absolute
#'    step of 1e-3 is wrong for a coefficient of 2e-5 on a covariate
#'    spanning 1e5 (its diagonal came out 1.3e96 against an exact 2.4e14)
#'    and leaves the domain of a parameter estimated at 5.9e-5 (a NaN
#'    row).
#' 3. What tier 2 cannot invert is a property of the likelihood, and
#'    se_tier3() decides which parameters it costs.
#'
#' Returns `list(V, shown, lost, null, tier)`: `V` the covariance kept
#' for propagation, `shown` the one vcov(), summary() and confint() read
#' (lost rows and columns NaN), `lost` the named reasons, `null` a basis
#' of the directions `V` leaves out (natural units, one column each),
#' which a prediction must not move along (jc_nonest()).
#'
#' @noRd
cov_from_hessian <- function(fit, H, E = NULL) {
  np <- NROW(H)
  none <- stats::setNames(character(0), character(0))
  if (!np) return(list(V = H, shown = H, lost = none, tier = 1L))
  V <- try(solve(H), silent = TRUE)
  if (cov_usable(V)) return(list(V = V, shown = V, lost = none, tier = 1L))
  p <- fit$opt$par
  saved <- obj_state_save(fit$obj)
  on.exit(obj_state_restore(fit$obj, saved), add = TRUE)
  # the noise in H, which decides what tier 3 may read off it: an AD
  # Hessian is exact to rounding, a finite-difference one is not
  exact <- FALSE
  if (!length(fit$obj$env$random)) {
    Hx <- tryCatch(fit$obj$he(p), error = function(e) NULL)
    if (is.matrix(Hx) && all(dim(Hx) == np) && all(is.finite(Hx))) {
      E <- abs(Hx - t(Hx)) / 2
      H <- (Hx + t(Hx)) / 2
      exact <- TRUE
    }
  }
  D <- sqrt(abs(diag(H)))
  if (all(is.finite(H)) && isTRUE(all(D > 0))) {
    S <- H / outer(D, D)
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    # solve() inverts a scaled matrix with a condition number of 1e12,
    # which on `a + log(c0)` gave a standard error of 6810 on an
    # estimate of 0.99 along a ridge the likelihood is flat on. This
    # also asks for positive definiteness, which tier 1 did not: tier 2
    # is reached only when tier 1's standard errors were not finite.
    Vs <- if (min(ev) > se_flat_tol * max(abs(ev))) {
      try(solve(S), silent = TRUE)
    }
    if (!is.null(Vs) && cov_usable(Vs)) {
      V <- Vs / outer(D, D)
      return(list(V = V, shown = V, lost = none, tier = 2L))
    }
  }
  if (is.null(E)) {
    E <- tryCatch(fd_hessian_noise(fit), error = function(e) NULL)
    if (is.null(E) || !all(dim(E) == np)) E <- matrix(0, np, np)
    E[!is.finite(E)] <- 0
  }
  se_tier3(fit, H, E, p, exact)
}

#' Tier 3 of cov_from_hessian(): which parameters have no standard
#' error, why, and the covariance of the rest.
#'
#' In order, each step only over the parameters the earlier ones left:
#'
#' - `"bound"`: a bound holds the parameter (grad_bound_active()). The
#'   fit is a constrained optimum and the others are conditional on it.
#' - `"nonfinite"`: its Hessian row is not finite.
#' - `"flat"`, empty row: its row is no larger than the Hessian's own
#'   noise (`se_noise_mult` times the row's noise; on a finite-difference
#'   Hessian at least `se_empty_row`). A row of noise scaled to unit
#'   diagonal has a diagonal of +1 or -1 and would show up as curvature of
#'   either sign.
#' - Then the unit-diagonal matrix of the rest is decomposed. When its
#'   smallest eigenvalue is above `se_flat_tol` of the largest it is
#'   inverted as it is. Otherwise every direction whose eigenvalue is
#'   not above its threshold (`se_flat_tol` of the largest, or
#'   `se_dir_mult` times that direction's own noise `||Es v||`) is
#'   removed: flat when its eigenvalue is within the threshold of zero,
#'   negative below that.
#' - A flat direction is one the data do not determine. A
#'   parameter loses its standard error to it when its unit vector has a
#'   projection onto the flat subspace above `tau`, the larger of
#'   sqrt(`se_flat_tol`) and the bound that the matrix's noise puts on an
#'   eigenvector's error (`se_noise_mult` times the noise norm over the
#'   eigenvalue gap). That has no dilution: a flat direction shared by
#'   120 coefficients loads 0.065 on each, far above tau, and every one
#'   of them is named.
#' - A negative direction is named by its dominant parameters
#'   (loading at least `se_dominant` of the largest). It is `"concave"`
#'   when a step along it raises the log-likelihood by more than
#'   `grad_tol` (se_line_probe()), and `"flat"` otherwise. A parameter
#'   with a small loading on it keeps the curvature of the directions
#'   that are a maximum: on frmtmb.sample's `gr(g, by = f)` fixture a
#'   slope SD at exp(-6.6) and its correlation load 0.70 each on an
#'   eigenvalue of -1.18 and `x` loads 0.10; `x` keeps an SE close to
#'   lme4's 0.1009 (lane nanse review, B1).
#'
#' The covariance kept for propagation is the pseudo-inverse over the
#' kept directions, zero in the rows of a removed parameter.
#'
#' @noRd
se_tier3 <- function(fit, H, E, p, exact = FALSE) {
  np <- NROW(H)
  nm <- outer_par_names(fit)
  why <- character(np)
  g <- tryCatch(drop(fit$obj$gr(p)), error = function(e) rep(NA_real_, np))
  why[grad_bound_active(p, g, fit_outer_box(fit)) %in% TRUE] <- "bound"
  # one parameter whose row is not finite also spoils a cell of every
  # other row, so rows go one at a time, the worst first
  repeat {
    f <- which(why == "")
    nbad <- rowSums(!is.finite(H[f, f, drop = FALSE]))
    if (!length(f) || !any(nbad > 0)) break
    why[f[which.max(nbad)]] <- "nonfinite"
  }
  f <- which(why == "")
  if (length(f)) {
    u <- fit$par_units %||% rep(1, np)
    Hu <- abs(H * outer(u, u))[f, f, drop = FALSE]
    Eu <- (E * outer(u, u))[f, f, drop = FALSE]
    # empty: within the row's own noise, which on an AD Hessian is
    # rounding, so that exact zeros are what it empties. A
    # finite-difference Hessian also errs by more than its asymmetry
    # shows (truncation, the inner solve), so there the measured floor
    # se_empty_row applies as well. Never relative to the largest entry:
    # a coefficient estimated at 5.9e-5 has a diagonal of 6.9e11 and
    # would empty sigma's row of 400
    floor_i <- se_noise_mult * apply(Eu, 1L, max)
    if (!exact) floor_i <- pmax(floor_i, se_empty_row)
    why[f[apply(Hu, 1L, max) <= floor_i]] <- "flat"
  }
  # a zero diagonal beside a row that is not empty: some combination
  # with another parameter curves downward
  why[why == "" & (diag(H) == 0) %in% TRUE] <- "concave"
  V <- matrix(0, np, np)
  # the removed rows are directions a prediction must not move along
  out_rows <- which(why %in% c("flat", "nonfinite", "concave"))
  null <- diag(1, np)[, out_rows, drop = FALSE]
  free <- which(why == "")
  ev <- numeric(0)
  if (length(free)) {
    Hf <- H[free, free, drop = FALSE]
    Df <- sqrt(abs(diag(Hf)))
    e <- eigen(Hf / outer(Df, Df), symmetric = TRUE)
    ev <- e$values
    big <- max(abs(ev))
    Es <- E[free, free, drop = FALSE] / outer(Df, Df)
    enorm <- sqrt(sum(Es^2))
    # each direction against its own noise, ||Es v_k||, which bounds how
    # far the noise can move that eigenvalue. The norm of the whole
    # block's noise did not: on (1 + x | g2) with no g2 variation and a
    # 40-level factor it was 0.033, above the identified intercept-
    # against-contrasts eigenvalue of 0.022 whose own noise is 0.0012,
    # and the intercept and 39 contrasts lost SEs lme4 reports (RB4)
    nk <- sqrt(colSums((Es %*% e$vectors)^2))
    thr <- pmax(se_flat_tol * big, se_dir_mult * nk)
    # once the empty, bound and non-finite rows are out, a positive
    # definite remainder is inverted as it is
    keep <- if (min(ev) > se_flat_tol * big) rep(TRUE, length(ev)) else {
      ev > thr
    }
    rem <- which(!keep)
    gap <- if (any(keep) && length(rem)) {
      min(ev[keep]) - max(ev[rem])
    } else big
    tau <- max(sqrt(se_flat_tol), se_noise_mult * enorm / max(gap, 1e-300))
    neg <- rem[ev[rem] < -thr[rem]]
    flat_k <- setdiff(rem, neg)
    conc <- vapply(neg, function(k) {
      se_line_probe(fit, p, free, e$vectors[, k] / Df, ev[k])
    }, NA)
    dominant <- function(k) {
      a <- abs(e$vectors[, k])
      a >= se_dominant * max(a)
    }
    mark <- character(length(free))
    # confirmed ascents first, so a parameter on one is called concave
    for (k in neg[order(!conc)]) {
      on <- dominant(k) & mark == ""
      mark[on] <- if (conc[match(k, neg)]) "concave" else "flat"
    }
    if (length(flat_k)) {
      proj <- sqrt(rowSums(e$vectors[, flat_k, drop = FALSE]^2))
      on <- proj > tau | Reduce(`|`, lapply(flat_k, dominant))
      mark[on & mark == ""] <- "flat"
    }
    why[free[mark != ""]] <- mark[mark != ""]
    if (length(rem)) {
      # a removed direction stands for the parameters it named; the
      # small loading of a parameter that keeps its standard error is
      # not part of it, or every prediction through that parameter would
      # lose its band (mo() seed 12: zeta1_2 at 0.022 on zeta2_2's
      # concave direction took every conditional_effects() band)
      Vr <- e$vectors[, rem, drop = FALSE]
      Vr[mark == "", ] <- 0
      Vr <- Vr[, colSums(Vr^2) > 0, drop = FALSE]
      Nk <- matrix(0, np, NCOL(Vr))
      Nk[free, ] <- Vr / Df
      null <- cbind(null, Nk)
    }
    if (any(keep)) {
      U <- e$vectors[, keep, drop = FALSE]
      Vf <- U %*% (t(U) / ev[keep])
      V[free, free] <- Vf / outer(Df, Df)
    }
  }
  lost <- why != ""
  shown <- V
  shown[lost, ] <- NaN
  shown[, lost] <- NaN
  # `ev` is kept for diagnosis: the spectrum the verdicts were read from
  list(V = V, shown = shown, lost = stats::setNames(why[lost], nm[lost]),
       null = null, tier = 3L, ev = ev)
}

#' The noise of the finite-difference outer Hessian of a fit with random
#' effects, when the caller did not keep it: rebuilt, only when tier 3
#' is reached, at the cost of the Hessian itself.
#'
#' @noRd
fd_hessian_noise <- function(fit) fit_outer_hessian(fit)$E

#' Is a direction of negative curvature a real ascent? One step each
#' way along it, sized so that the quadratic model predicts a gain of
#' twice `grad_tol` in log-likelihood; the direction counts as concave
#' when the better step really gains more than `grad_tol`. On a
#' finite-difference Hessian a near-zero standard deviation's row turns
#' into an eigenvalue of -0.16 to -10.6 that the objective does not
#' follow (lane nanse review, m1: at most 1.2e-8 along it).
#'
#' @noRd
se_line_probe <- function(fit, p, free, dir_free, lambda) {
  tol <- fit$control$grad_tol %||% 1e-3
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  d <- numeric(length(p))
  d[free] <- dir_free * sqrt(4 * tol / abs(lambda))
  f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
  fu <- tryCatch(obj$fn(p + d), error = function(e) NA_real_)
  fd <- tryCatch(obj$fn(p - d), error = function(e) NA_real_)
  gain <- f0 - min(fu, fd, na.rm = TRUE)
  is.finite(gain) && gain > tol
}

#' A direction of the unit-diagonal Hessian with an eigenvalue at most
#' this fraction of the largest is one the likelihood does not curve
#' along. Only read after both inverses have failed. Exact ridges give
#' 1e-16 to 2e-16 there and a mo() simplex weight at 0 gives an exactly
#' zero diagonal or a decoupled row, so the value is not delicate; it
#' is the threshold nl_flat_message() uses for the same question. Its
#' square root is the smallest projection onto a flat subspace that
#' takes a parameter's standard error.
#'
#' @noRd
se_flat_tol <- 1e-9

#' A quantity this many times the Hessian's own noise is real; below it,
#' it is noise. The noise is the asymmetry of the Hessian before it is
#' symmetrized (fd_hessian_noise(), or the AD Hessian's rounding).
#'
#' @noRd
se_noise_mult <- 10

#' A direction of the unit-diagonal Hessian is noise when its eigenvalue
#' is within this many times its own noise `||Es v||` of zero. That
#' noise bounds the eigenvalue's error (Weyl, restricted to the
#' direction), so a small multiple is enough; the block's whole noise
#' norm was not a bound on any one direction and took identified ones
#' (lane nanse review, RB4).
#'
#' @noRd
se_dir_mult <- 3

#' A finite-difference Hessian row whose largest entry is at most this,
#' in the optimizer's units, is empty whatever its asymmetry says.
#' Measured on the 47 fits of the eight test suites where tier 3 ran in
#' lane nanse's first round (dev/nanse-rowmax-sum.R): the 71 rows it
#' emptied had largest entries of at most 5.6e-7, and the smallest row
#' it kept was 1.7e-6, a nearly saturated simplex coordinate. Not
#' applied to an exact AD Hessian, whose rows are emptied only at its
#' rounding.
#'
#' @noRd
se_empty_row <- 1e-6

#' On a removed direction, the parameters whose loading is at least this
#' fraction of the largest are its own and lose their standard error.
#' Flat directions also apply the projection rule of se_tier3(), which a
#' direction shared by many parameters cannot dilute.
#'
#' @noRd
se_dominant <- 0.5

#' Put the covariance of `cov_from_hessian()` into an sdreport, and
#' record which standard errors are lost and the directions they span.
#'
#' The fit-time check's own analysis wins when it found something
#' (`fit$cache$se_analysis`): on a model without random effects it was
#' read off the exact Hessian, and sdreport()'s finite-difference
#' Hessian can invert a matrix the exact one shows to be singular, which
#' would print finite standard errors under a warning that says they do
#' not exist. Otherwise sdreport()'s own inverse is kept whenever it is
#' usable.
#'
#' Random-effect standard errors (`diag.cov.random`) are left as
#' sdreport() computed them: they are built from the same failed
#' inverse inside sdreport(), and redoing that is outside this repair.
#'
#' @noRd
sdr_rescue <- function(fit, sdr, H = NULL) {
  V0 <- sdr$cov.fixed
  if (!length(V0)) return(sdr)
  an <- if (is.environment(fit$cache)) fit$cache$se_analysis
  if (is.null(an)) {
    if (cov_usable(V0)) return(sdr)
    # A fit that did not converge keeps sdreport()'s own answer: its
    # curvature is not that of an optimum, and its convergence warning
    # has already said what is wrong with it
    if (!isTRUE(fit$opt$convergence == 0) ||
          identical(fit$cache$se_explained, "convergence")) {
      return(sdr)
    }
    hc <- if (is.environment(fit$cache)) fit$cache$hessian_fixed
    E <- if (!is.null(H) && identical(H, hc$H)) hc$E
    if (is.null(H)) {
      h <- tryCatch(fit_outer_hessian(fit), error = function(e) NULL)
      H <- h$H
      E <- h$E
    }
    if (is.null(H)) return(sdr)
    an <- cov_from_hessian(fit, H, E)
    if (an$tier == 1L) return(sdr)
  }
  if (!all(dim(an$V) == dim(V0))) return(sdr)
  dn <- dimnames(V0)
  # sdreport()'s own inverse, which diagnose() reads its smallest
  # eigenvalue from
  sdr$cov_fixed_raw <- V0
  sdr$cov.fixed <- an$shown
  dimnames(sdr$cov.fixed) <- dn
  sdr$cov_fixed_prop <- an$V
  dimnames(sdr$cov_fixed_prop) <- dn
  sdr$se_lost <- se_relabel(fit, an$lost)
  sdr$se_null <- an$null
  sdr
}

#' The fit-time standard-error check.
#'
#' Runs on the user's own call, after every other fit-time verdict. An
#' optimizer that did not report convergence or a gradient that says the
#' fit stopped short is not a point whose curvature means anything, so
#' such a fit is left to its convergence warning; an earlier check that
#' named some parameters explains those parameters only
#' (se_explained_pars()).
#'
#' Cost. On a model without random effects the check reads the exact AD
#' Hessian, and the finite-difference one only when that shows a problem.
#' That is one gradient sweep per coefficient, so its share of the fit
#' grows with the number of coefficients (1.19 to 1.24 at 301). With
#' random effects it builds the finite-difference Hessian (fd_hessian(),
#' two gradients per outer parameter), which a later summary() reuses,
#' when se_check_at_fit() says so: a decision on counted work, never on
#' the clock, so the same call warns at the same place on any machine
#' (lane nanse review, RB1: a wall-clock budget deferred the same fit 1
#' time in 20 idle and 3 in 20 under load). Otherwise the check waits
#' for the first standard-error use, which builds the same Hessian
#' anyway and warns then (se_deferred_report()). The costs measured
#' against `check_se = "ignore"` are in the findings (dev/nanse-cost2.R,
#' dev/nanse-rev-cost.R). `frmtmb_control(check_se = "ignore")` skips
#' it.
#'
#' @noRd
se_check <- function(fit, control) {
  act <- control$check_se %||% "warning"
  if (identical(act, "ignore") || !length(fit$opt$par) ||
        !is.environment(fit$cache)) {
    return(invisible(NULL))
  }
  if (identical(fit$cache$se_explained, "convergence")) {
    # the convergence warning said it; vcov() must not say it again in
    # older words
    fit$cache$warned_nonfinite_cov <- TRUE
    return(invisible(NULL))
  }
  lost <- if (!is.null(fit$cache$sdr)) {
    fit$cache$sdr$se_lost
  } else {
    an <- tryCatch(se_fit_analysis(fit), error = function(e) e)
    if (inherits(an, "error")) {
      re_check_act(act, paste0("The standard errors could not be ",
                               "computed: building the Hessian at the ",
                               "optimum failed (", conditionMessage(an),
                               ")"))
      return(invisible(NULL))
    }
    if (is.null(an)) {
      fit$cache$se_deferred <- act
      return(invisible(NULL))
    }
    an$lost
  }
  se_report(fit, lost, act)
}

#' The fit-time analysis behind se_check(), or NULL when
#' se_check_at_fit() leaves it to the first standard-error use. Caches
#' what a later sdreport() reuses: the finite-difference Hessian
#' (`hessian_fixed`) and, when it found something, the analysis itself
#' (`se_analysis`).
#'
#' @noRd
se_fit_analysis <- function(fit) {
  obj <- fit$obj
  p <- fit$opt$par
  if (!length(obj$env$random)) {
    saved <- obj_state_save(obj)
    Hx <- tryCatch(obj$he(p), error = function(e) NULL)
    obj_state_restore(obj, saved)
    if (is.matrix(Hx) && all(dim(Hx) == length(p)) && all(is.finite(Hx))) {
      Hs <- (Hx + t(Hx)) / 2
      if (cov_usable(try(solve(Hs), silent = TRUE))) {
        return(list(lost = character(0)))
      }
      an <- cov_from_hessian(fit, Hs, abs(Hx - t(Hx)) / 2)
      if (an$tier > 1L) {
        an$lost <- se_relabel(fit, an$lost)
        fit$cache$se_analysis <- an
      }
      return(an)
    }
  }
  if (!se_check_at_fit(fit)) return(NULL)
  h <- fit_outer_hessian(fit)
  fit$cache$hessian_fixed <- h
  an <- cov_from_hessian(fit, h$H, h$E)
  if (an$tier > 1L) {
    an$lost <- se_relabel(fit, an$lost)
    fit$cache$se_analysis <- an
  }
  an
}

#' Does the fit-time check build the finite-difference Hessian of a fit
#' with random effects now, or leave it to the first standard-error use?
#'
#' Counted work, so the answer is the same on every run and every
#' machine: the build costs two gradients per outer parameter, and it
#' runs when that is at most `se_check_share` of the objective and
#' gradient evaluations the optimizer made (`fit$opt$evals`, counted by
#' optimize_obj()), or when the fit has at most `se_check_np_free` outer
#' parameters. A fit without a count (a path that did not go through
#' optimize_obj()) has only the second rule.
#'
#' @noRd
se_check_at_fit <- function(fit) {
  np <- length(fit$opt$par)
  if (np <= se_check_np_free) return(TRUE)
  ev <- fit$opt$evals
  is.numeric(ev) && length(ev) == 1L && is.finite(ev) &&
    2 * np <= se_check_share * ev
}

#' A gradient costs more than an objective value (on a Laplace
#' objective about two to four times), so a quarter of the counted
#' evaluations keeps the build under about a third of the fit's own
#' optimization work; see the findings for the measured costs.
#'
#' @noRd
se_check_share <- 0.25

#' Twenty gradients: what the smallest models always pay, so that
#' ordinary small fits are checked where they are made.
#'
#' @noRd
se_check_np_free <- 10L

#' Warn about the lost standard errors that no earlier check explained.
#'
#' @noRd
se_report <- function(fit, lost, act) {
  if (!length(lost)) return(invisible(NULL))
  # vcov() would otherwise say it again
  fit$cache$warned_nonfinite_cov <- TRUE
  ex <- se_explained_pars(fit)
  # a family's own verdict is about identification, so it explains every
  # lost parameter but one a bound holds, which it does not mention
  lost <- lost[!names(lost) %in% ex$pars & !(ex$family & lost != "bound")]
  if (!length(lost)) return(invisible(NULL))
  re_check_act(act, se_lost_message(fit, lost))
  invisible(NULL)
}

#' The deferred half of se_check(): the first standard-error use of a fit
#' whose check waited (se_check_at_fit()) reports what it finds, once.
#'
#' Once means once to the CALLER. Several readers of the covariance run
#' it inside their own suppressWarnings() (fixef() reads vcov() that
#' way), and a report raised there was muffled and then marked as given,
#' so the fit lost its warning for good (lane nanse review, RB1). Inside
#' such a call the report stays pending for the next standard-error use;
#' the readers that suppress call se_flush_deferred() first, so that
#' their own caller gets it.
#'
#' @noRd
se_deferred_report <- function(fit, sdr) {
  act <- if (is.environment(fit$cache)) fit$cache$se_deferred
  if (is.null(act) || se_muffled_inside()) return(invisible(NULL))
  fit$cache$se_deferred <- NULL
  se_report(fit, sdr$se_lost, act)
}

#' Is this call running inside a suppressWarnings() that frmtmb or one of
#' its extensions wrote? A user's own suppressWarnings() is the user
#' declining the warning, which counts as delivered.
#'
#' @noRd
se_muffled_inside <- function() {
  sw <- base::suppressWarnings
  par <- sys.parents()
  for (i in seq_len(sys.nframe())) {
    if (!identical(sys.function(i), sw) || par[i] < 1L) next
    # the function that called it, not the environment it ran in: a
    # testthat test environment also has the namespace above it
    te <- environment(sys.function(par[i]))
    if (is.environment(te) && isNamespace(te) &&
          startsWith(getNamespaceName(te), "frmtmb")) {
      return(TRUE)
    }
  }
  FALSE
}

#' Deliver a pending deferred report before a covariance reader that
#' runs under suppressWarnings(), so that the report reaches its caller.
#' Errors are left to the reader, which handles them its own way.
#'
#' @noRd
se_flush_deferred <- function(fit) {
  if (is.environment(fit$cache) && !is.null(fit$cache$se_deferred)) {
    tryCatch(sdr_of(fit), error = function(e) NULL)
  }
  invisible(NULL)
}

#' What other checks have already explained: `list(family, pars)`.
#'
#' Precedence is per parameter. An earlier warning explains the
#' parameters it is about and no others:
#'
#' - the structural checks that run before the fit name theirs, and
#'   record them as the frame attribute `se_explained` only when they
#'   warned, not under "ignore": a grouping factor with one level
#'   explains that term's standard deviation, an observation-level
#'   effect its standard deviation and the residual sd (`pars`);
#' - a family's fit-end check (fit_end_checks()) cannot name outer
#'   parameters, and what it reports is identification, so it explains
#'   every lost parameter except one a bound holds (`family`);
#' - the convergence verdict explains the whole fit, and se_check()
#'   returns before asking.
#'
#' @noRd
se_explained_pars <- function(fit) {
  fam <- !is.null(fit$cache$se_explained)
  ent <- attr(fit$frame, "se_explained")
  if (!length(ent)) return(list(family = fam, pars = character(0)))
  om <- outer_par_map(fit)
  th <- om$names[om$comp == "theta"]
  bd <- names(fit$frame[["par_template"]][["betad"]])
  out <- character(0)
  for (e in ent) {
    out <- c(out, th[e$theta_idx])
    if (identical(e$kind, "olre")) {
      lp <- fit$frame[["linpreds"]][[linpred_key(e$resp, "sigma")]]
      if (!is.null(lp) && identical(lp[["par"]], "betad")) {
        out <- c(out, bd[lp[["idx"]]])
      }
    }
  }
  list(family = fam, pars = out[!is.na(out)])
}

#' A flat mean coefficient of a binomial-type fit that has run far out
#' on the link is separation, which says what happened better than
#' "flat" does (lane nanse review, m7).
#'
#' @noRd
se_relabel <- function(fit, lost) {
  if (!length(lost) || fit$REML || isTRUE(fit$control$profile)) {
    return(lost)
  }
  om <- outer_par_map(fit)
  pb <- which(om$comp == "beta")
  sep <- character(0)
  for (lp in fit$frame[["linpreds"]]) {
    fam <- fit$spec$responses[[lp[["resp"]]]]$family
    if (!identical(lp[["par"]], "beta") ||
          !fam[["family"]] %in% separation_families ||
          !lp[["dpar"]] %in% (fam[["primary_dpars"]] %||% "mu")) next
    pos <- pb[lp[["idx"]]]
    pos <- pos[!is.na(pos)]
    sep <- c(sep, om$names[pos][abs(fit$opt$par[pos]) > 10])
  }
  hit <- names(lost) %in% sep & lost == "flat"
  lost[hit] <- "separation"
  lost
}

#' One clause per reason, naming the parameters. A mo() simplex
#' coordinate also names its term, since `zeta2_1` alone says nothing.
#'
#' @noRd
se_lost_clauses <- function(fit, lost) {
  mo_lab <- character(0)
  for (lp in fit$frame[["linpreds"]]) {
    for (mi in lp[["mo"]] %||% list()) mo_lab[mi$zeta] <- mi$label
  }
  label <- function(n) {
    z <- sub("_[0-9]+$", "", n)
    if (z %in% names(mo_lab)) {
      paste0(n, " (simplex of ", mo_lab[[z]], ")")
    } else n
  }
  text <- function(r, one) {
    it <- if (one) "it" else "them"
    switch(r,
      bound = paste0("a bound holds ", it, ", so the fit is a constrained ",
                     "optimum there"),
      separation = paste0("the data separate the outcomes, so ",
                          if (one) "the estimate runs" else
                            "the estimates run",
                          " off toward infinity and the optimizer stopped ",
                          "where its tolerances did"),
      flat = paste0("the likelihood is flat along ", it, " at the ",
                    "estimates, so the data do not determine ", it,
                    " there (a standard deviation or a mo() simplex weight ",
                    "at 0 does this, and so do parameters that enter only ",
                    "through a combination)",
                    # flat_par_note()'s remedy, for the model it is about
                    if (fit_has_nlpars(fit)) {
                      paste0(". A nonlinear term that has left its own ",
                             "support does this too; give it a starting ",
                             "value that puts it back (see par_template())")
                    }),
      concave = paste0("a step along ", it, " raises the log-likelihood, ",
                       "so the estimates are not a maximum in that ",
                       "direction and other starting values may reach a ",
                       "higher one"),
      nonfinite = paste0("the Hessian is not finite in ",
                         if (one) "its row" else "their rows"))
  }
  out <- character(0)
  for (r in c("bound", "separation", "flat", "concave", "nonfinite")) {
    n <- names(lost)[lost == r]
    if (!length(n)) next
    out <- c(out, paste0(paste(vapply(n, label, ""), collapse = ", "),
                         ": ", text(r, length(n) == 1L)))
  }
  out
}

#' The warning.
#'
#' Under REML or `profile = TRUE` the coefficients' standard errors come
#' from the joint precision, which this repair does not reach, so the
#' warning does not promise them (lane nanse review, m5).
#'
#' @noRd
se_lost_message <- function(fit, lost) {
  np <- length(fit$opt$par)
  paste0("Standard errors are not available for ", length(lost), " of ",
         np, " parameters. ", paste(se_lost_clauses(fit, lost),
                                    collapse = ". "),
         if (needs_jp(fit)) {
           paste0(". The other outer parameters keep theirs; the ",
                  "coefficients' standard errors come from the joint ",
                  "precision, which this does not repair")
         } else {
           paste0(". The other standard errors are kept, and a ",
                  "prediction that moves along these directions gets ",
                  "none")
         },
         "; summary() and diagnose() list them. See the ",
         "'Convergence problems' section of vignette('diagnostics')")
}

#' Rows of a prediction gradient that move along a direction the fit
#' does not determine.
#'
#' `jc$null` (sdr_rescue(), through get_joint_cov()) spans the directions
#' the propagation covariance leaves out. A prediction whose gradient
#' `G` (rows, over the joint-covariance positions `pos`) has a component
#' along one of them has no finite standard error, and the covariance
#' would otherwise give it a finite and too small one (lane nanse
#' review, B2). The component is a cosine in the optimizer's units, so
#' it does not depend on how a coefficient is scaled; a prediction from a
#' saturated simplex weight has an exactly zero component and keeps its
#' band. This is predict()'s `alias_null` estimability test, for the
#' directions the Hessian lost rather than the design.
#'
#' @noRd
jc_nonest <- function(jc, G, pos) {
  N <- jc$null
  n <- NROW(G)
  if (is.null(N) || !NCOL(N) || !n || !length(pos)) return(rep(FALSE, n))
  u <- jc$units %||% rep(1, NROW(N))
  Nq <- N / u
  Gq <- as.matrix(G) * rep(u[pos], each = n)
  num <- abs(Gq %*% Nq[pos, , drop = FALSE])
  den <- outer(sqrt(rowSums(Gq^2)), sqrt(colSums(Nq^2)))
  rowSums(num > se_pred_tol * den, na.rm = TRUE) > 0
}

#' The cosine above which a prediction counts as moving along a lost
#' direction. The directions come from an exact AD Hessian on the only
#' models that use them (no random effects), so they are accurate to
#' rounding and the test can be tight.
#'
#' @noRd
se_pred_tol <- 1e-6

#' The one warning a call gives when some of its predictions (or
#' hypotheses) lost their standard error to jc_nonest().
#'
#' @noRd
se_pred_warn <- function(bad, what = "predictions") {
  if (!any(bad)) return(invisible(NULL))
  frm_warning(sum(bad), " of ", length(bad), " ", what, " move along a ",
              "direction the fit does not determine (a parameter without ",
              "a standard error), so their standard errors are NaN; ",
              "summary() lists those parameters", call. = FALSE,
              class = "frmtmb_se_lost_prediction")
}
