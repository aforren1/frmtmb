# Standard errors the outer Hessian cannot give: found at fit time,
# recovered where the Hessian still carries them, and said once.
#
# sdreport() inverts the outer Hessian with solve(), which refuses a
# matrix whose reciprocal condition number is below machine epsilon and
# then fills EVERY standard error with NaN. On `ls ~ mo(income) * age`
# (brms_monotonic's own data code, dev/nanse-mo-sweep.R) that happened
# on 72 of 200 data sets, 55 of them with optimizer code 0 and no
# warning: the interaction's simplex sat with a weight at 0, its
# softmax coordinate had run to -20 or beyond, and the likelihood was
# flat along that one coordinate to 1e-20 (frmtmb 0.68.1; lane optima
# moved the simplex to a chart that reaches a face at a finite
# coordinate, dev/optima-findings.md). The rest of the matrix is
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
#' Tiers 1 and 2 are kept only when se_tier3_could_act() finds nothing
#' tier 3 would remove. Without that, whether a parameter had a standard
#' error depended on the sign of noise: a log sd at exp(-29) on
#' `(1 | Subject/a)` has a finite-difference Hessian row of +7.1e-12 with
#' the reference BLAS and -7.1e-12 with OpenBLAS, so solve() gave it a
#' standard error of 375,150 on one platform and NaN on the other. Over
#' the eight test suites 117 fits accepted such a row with no warning,
#' 5 of 104 of them changing verdict with the BLAS
#' (dev/reviews/2026-10-06-cifix.md, item 4); and a nonlinear ridge with
#' a random effect kept a prediction standard error of 822,571. When
#' tier 3 then removes nothing, the tier 1 inverse is still the one
#' returned, so a healthy fit keeps sdreport()'s covariance bit for bit.
#'
#' `exact` says `H` is the AD Hessian (rounding is its only noise).
#'
#' Returns `list(V, shown, lost, null, tier)`: `V` the covariance kept
#' for propagation, `shown` the one vcov(), summary() and confint() read
#' (lost rows and columns NaN), `lost` the named reasons, `null` a basis
#' of the directions `V` leaves out (natural units, one column each),
#' which a prediction must not move along (jc_nonest()).
#'
#' @noRd
cov_from_hessian <- function(fit, H, E = NULL, exact = FALSE) {
  np <- NROW(H)
  none <- stats::setNames(character(0), character(0))
  if (!np) return(list(V = H, shown = H, lost = none, tier = 1L))
  V <- try(solve(H), silent = TRUE)
  tier1 <- if (cov_usable(V)) {
    list(V = V, shown = V, lost = none, tier = 1L)
  }
  p <- fit$opt$par
  held <- se_bound_held(fit, p)
  edge <- se_edge_sd(fit, H)
  if (!is.null(tier1) &&
        !se_tier3_could_act(fit, H, E, exact, held, edge)) {
    return(tier1)
  }
  saved <- obj_state_save(fit$obj)
  on.exit(obj_state_restore(fit$obj, saved), add = TRUE)
  # the noise in H, which decides what tier 3 may read off it: an AD
  # Hessian is exact to rounding, a finite-difference one is not
  replaced <- FALSE
  if (!exact && !length(fit$obj$env$random)) {
    Hx <- tryCatch(fit$obj$he(p), error = function(e) NULL)
    if (is.matrix(Hx) && all(dim(Hx) == np) && all(is.finite(Hx))) {
      E <- abs(Hx - t(Hx)) / 2
      H <- (Hx + t(Hx)) / 2
      exact <- TRUE
      replaced <- TRUE
    }
  }
  D <- sqrt(abs(diag(H)))
  # tier 2 only where tier 1 refused; a tier 1 inverse that tier 3 might
  # act on goes straight to tier 3, which is what decides it
  if (is.null(tier1) && all(is.finite(H)) && isTRUE(all(D > 0))) {
    S <- H / outer(D, D)
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    # solve() inverts a scaled matrix with a condition number of 1e12,
    # which on `a + log(c0)` gave a standard error of 6810 on an
    # estimate of 0.99 along a ridge the likelihood is flat on. This
    # also asks for positive definiteness, which tier 1 did not: tier 2
    # is reached only when tier 1's standard errors were not finite.
    Vs <- if (min(ev) > se_flat_tol * max(abs(ev)) &&
                !se_tier3_could_act(fit, H, E, exact, held, edge)) {
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
  an <- se_tier3(fit, H, E, p, exact, held, edge)
  if (!length(an$lost) && !is.null(tier1) && !replaced) tier1 else an
}

#' Could se_tier3() remove a parameter from this Hessian? Asked before
#' tier 1 or tier 2 is accepted, with tier 3's own two tests, so that the
#' verdict does not depend on which tier ran:
#'
#' - a row no larger than the Hessian's noise (tier 3's empty-row rule,
#'   the same floors, in the optimizer's units);
#' - a unit-diagonal eigenvalue at or below `se_flat_tol` of the largest.
#'   Tier 3 inverts a remainder whose smallest eigenvalue is above that
#'   as it is, so nothing else can make it remove a direction.
#'
#' Two more of tier 3's reasons send it there too, so that they do not
#' depend on an unrelated eigenvalue: a parameter a bound holds (`held`;
#' a raw polynomial's top coefficient under an active `ub` was reported
#' at degree 5 and not at degree 3, dev/setier-rev-bound.R), and a log
#' standard deviation at its edge (`edge`, se_edge_sd()). A `TRUE` only
#' sends the Hessian to tier 3; tier 3 decides.
#'
#' @noRd
se_tier3_could_act <- function(fit, H, E, exact, held = logical(0),
                               edge = logical(0)) {
  np <- NROW(H)
  if (any(held) || any(edge)) return(TRUE)
  if (!all(is.finite(H))) return(TRUE)
  if (is.null(E) || !all(dim(E) == np)) E <- matrix(0, np, np)
  E[!is.finite(E)] <- 0
  u <- fit$par_units %||% rep(1, np)
  floor_i <- se_noise_mult * apply(E * outer(u, u), 1L, max)
  if (!exact) floor_i <- pmax(floor_i, se_empty_row)
  if (any(apply(abs(H * outer(u, u)), 1L, max) <= floor_i)) return(TRUE)
  D <- sqrt(abs(diag(H)))
  if (!all(D > 0)) return(TRUE)
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  min(ev) <= se_flat_tol * max(abs(ev))
}

#' The parameters a bound holds (grad_bound_active()), read only when the
#' box has a finite bound, since a gradient costs a sweep.
#'
#' @noRd
se_bound_held <- function(fit, p) {
  np <- length(p)
  box <- fit_outer_box(fit)
  if (!any(is.finite(c(box$lower, box$upper)))) return(logical(np))
  g <- tryCatch(drop(fit$obj$gr(p)), error = function(e) rep(NA_real_, np))
  grad_bound_active(p, g, box) %in% TRUE
}

#' Log standard deviations at their edge: the likelihood does not change
#' when one moves further toward zero (se_at_edge()). Asked of every log
#' sd whose curvature could allow it, because the empty-row floor is
#' absolute and the row of an sd at its boundary grows with n: on 100
#' groups of 10 an lme4-singular sd had a row of 4.99e-6, a standard
#' error of 1515 and no report (dev/setier-rev-smallsd.R, seed 10). Near
#' zero the log-likelihood goes as c - a exp(2 theta), so a move of 2
#' down changes it by about a quarter of the curvature; a curvature above
#' `se_edge_screen` (in the optimizer's units) rules the edge out without
#' an evaluation.
#'
#' @noRd
se_edge_sd <- function(fit, H) {
  np <- NROW(H)
  out <- logical(np)
  # Under quadrature the boundary is not asked: the upward check
  # (se_sd_gain_up()) cannot run on that objective, and without it a fit
  # stuck near sd = 0 is called a boundary (4 of 4 rescaled-response fits
  # 0.0007 to 0.825 below lme4, 20 of 20 GLMMs started at exp(-12),
  # setier review, final check RQ1); such an sd keeps base's "flat"
  # verdict and its warning
  if (isTRUE(fit$quadrature)) return(out)
  th <- fit$estimates[["theta"]]
  if (!length(th)) return(out)
  om <- outer_par_map(fit)
  thpos <- which(om$comp == "theta")
  if (length(thpos) != length(th) || length(om$names) != np) return(out)
  sd_i <- tryCatch(log_sd_theta_index(fit), error = function(e) integer(0))
  u <- fit$par_units %||% rep(1, np)
  for (k in unique(sd_i)) {
    j <- thpos[k]
    hd <- abs(H[j, j]) * u[j]^2
    if (!is.finite(hd) || hd > se_edge_screen) next
    out[j] <- se_at_edge(fit, om$names[j], -1) &&
      se_sd_gain_up(fit, k) <= se_edge_tol
  }
  # the one correlation of a block whose sd is at zero is undefined
  # there, however much curvature is left in it: on the gr(g, by = f)
  # fixture the slope sd at exp(-6.6) left its correlation a finite
  # standard error once the sd was out. A wider block's correlations
  # cannot be told apart here and keep their own verdicts.
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    ti <- bk[["theta_idx"]]
    ti <- ti[!is.na(ti) & ti >= 1L & ti <= length(th)]
    other <- setdiff(ti, sd_i)
    if (length(other) == 1L && any(out[thpos[intersect(ti, sd_i)]])) {
      out[thpos[other]] <- TRUE
    }
  }
  out
}

#' The curvature above which se_edge_sd() does not ask: 25 times the
#' largest curvature a move of 2 toward zero within `grad_tol` (1e-3)
#' allows.
#'
#' @noRd
se_edge_screen <- 0.1

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
se_tier3 <- function(fit, H, E, p, exact = FALSE, held = NULL,
                     edge = logical(0)) {
  np <- NROW(H)
  nm <- outer_par_names(fit)
  why <- character(np)
  if (is.null(held)) held <- se_bound_held(fit, p)
  why[held %in% TRUE] <- "bound"
  # a log sd the likelihood ignores below its estimate is at its edge
  # whatever the size of its row, which grows with n (se_edge_sd())
  if (length(edge) == np) why[why == "" & edge] <- "flat"
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
    # A small positive eigenvalue can be real curvature: a raw polynomial
    # of degree 5 on [1, 2] has one at 1e-11 of the largest, and lm()
    # and lme4 report every standard error of it. With (1 | g) and degree
    # 6 the finite-difference eigenvalue is 2.1e-13 of the largest
    # against an asymmetry noise of 1.9e-12, and sdreport()'s inverse
    # still matched lme4's standard errors to 0.4 percent
    # (dev/setier-collin.R, dev/setier-deg6.R), so the noise bound
    # cannot decide it. A step along the direction that the quadratic
    # model says loses 2 * grad_tol of log-likelihood does: a ridge loses
    # nothing (se_curvature_real()). Below 100 machine epsilons of the
    # largest the eigenvalue itself has no digits left to trust.
    try_k <- which(!keep & ev > 100 * .Machine$double.eps * big)
    for (k in try_k) {
      keep[k] <- se_curvature_real(fit, p, free, e$vectors[, k] / Df, ev[k])
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
  # both steps can leave the domain (c0 + exp(a)^k seeds 4 and 8), and
  # min() of nothing warns "no non-missing arguments"
  best <- c(fu, fd)[is.finite(c(fu, fd))]
  gain <- if (length(best)) f0 - min(best) else NA_real_
  is.finite(gain) && gain > tol
}

#' Is a small positive eigenvalue real curvature? Symmetric second
#' differences along its direction, `c(t) = f(p + t) + f(p - t) - 2 f(p)`,
#' at a full step sized so that the quadratic model predicts a loss of
#' twice `grad_tol` each way (`c` = 4 `grad_tol`) and at half of it. The
#' curvature is real when `c(full)` is at least 2 `grad_tol` and
#' `c(full) / c(half)` is between 3 and 5.5 (a quadratic gives 4).
#'
#' Symmetric, because the fit need not be stationary along a weak
#' direction: on a poisson pair at a correlation of 1 - 5e-11 the
#' gradient's linear term moved the two full steps by +0.072 and -0.068
#' against a curvature loss of 0.002, and a rule on each side's own loss
#' took the standard errors of 20 of 20 identified fits, which glm()
#' reports (setier review, re-check RB1; dev/setier-rev2-window.R). The
#' second difference cancels that term. The ratio keeps out a straight
#' step off a curved ridge, which loses as the fourth power of the step
#' (ratio 16): on c0 + exp(a)^k with a ~ 1 + (1 | g), seed 77, the ridge
#' direction lost 5.6e-3 at a ratio of 15.5, while the objective moved
#' at most 5.2e-11 along the ridge itself (review B1). Along a straight
#' ridge the likelihood does not move, and a step that leaves the
#' domain is not a confirmation.
#'
#' @noRd
se_curvature_real <- function(fit, p, free, dir_free, lambda) {
  tol <- fit$control$grad_tol %||% 1e-3
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  d <- numeric(length(p))
  d[free] <- dir_free * sqrt(4 * tol / lambda)
  f <- function(x) tryCatch(obj$fn(x), error = function(e) NA_real_)
  f0 <- f(p)
  full <- f(p + d) + f(p - d) - 2 * f0
  half <- f(p + d / 2) + f(p - d / 2) - 2 * f0
  r <- full / half
  all(is.finite(c(full, half, r))) && full >= 2 * tol &&
    r > se_quad_lo && r < se_quad_hi
}

#' The band of full-step over half-step second difference that
#' se_curvature_real() takes as quadratic (4 exactly); the fourth power
#' of a curved ridge gives 16.
#'
#' @noRd
se_quad_lo <- 3

#' @noRd
se_quad_hi <- 5.5

#' A direction of the unit-diagonal Hessian with an eigenvalue at most
#' this fraction of the largest is one the likelihood does not curve
#' along. se_tier3_could_act() reads it before either inverse is kept.
#' Exact ridges give 1e-16 to 2e-16 there and a mo() simplex coordinate
#' the likelihood does not read gives an exactly zero diagonal or a
#' decoupled row (a softmax weight run to 0 did too, before lane
#' optima), so the value is not delicate; it is the threshold
#' nl_flat_message() uses for the same question. Its
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

#' The coefficients the separation warning named (check_convergence())
#' lose their standard error on every path, as the warning says: NaN in
#' what vcov() shows, held at their estimates in what predictions read,
#' reason "separation". Without this, a fit that stopped at nlminb's
#' limit kept sdreport()'s values, 6.2e132 for x with the reference BLAS
#' (setier review, m3), and a probit fit at code 0 whose Hessian was not
#' finite kept them NaN with no reason.
#'
#' @noRd
sdr_sep_lost <- function(fit, sdr, named) {
  V0 <- sdr$cov.fixed
  nm <- outer_par_names(fit)
  j <- match(named, nm)
  j <- j[!is.na(j)]
  if (!length(j) || length(nm) != NROW(V0)) return(sdr)
  lost <- sdr$se_lost %||% stats::setNames(character(0), character(0))
  j <- setdiff(j, match(names(lost), nm))
  if (!length(j)) return(sdr)
  dn <- dimnames(V0)
  prop <- sdr$cov_fixed_prop %||% V0
  prop[j, ] <- 0
  prop[, j] <- 0
  shown <- V0
  shown[j, ] <- NaN
  shown[, j] <- NaN
  sdr$cov_fixed_raw <- sdr$cov_fixed_raw %||% V0
  sdr$cov.fixed <- shown
  sdr$cov_fixed_prop <- prop
  dimnames(sdr$cov.fixed) <- dimnames(sdr$cov_fixed_prop) <- dn
  sdr$se_lost <- c(lost, stats::setNames(rep("separation", length(j)),
                                         nm[j]))
  sdr$se_null <- cbind(sdr$se_null %||% matrix(0, length(nm), 0L),
                       diag(1, length(nm))[, j, drop = FALSE])
  sdr
}

#' Put the covariance of `cov_from_hessian()` into an sdreport, and
#' record which standard errors are lost and the directions they span.
#'
#' The fit-time check's own analysis wins when it found something
#' (`fit$cache$se_analysis`): on a model without random effects it was
#' read off the exact Hessian, and sdreport()'s finite-difference
#' Hessian can invert a matrix the exact one shows to be singular, which
#' would print finite standard errors under a warning that says they do
#' not exist. Otherwise sdreport()'s own inverse is kept whenever
#' cov_from_hessian() would keep it: usable, and on a fit with random
#' effects with nothing tier 3 would remove.
#'
#' Random-effect standard errors (`diag.cov.random`) are left as
#' sdreport() computed them: they are built from the same failed
#' inverse inside sdreport(), and redoing that is outside this repair.
#'
#' @noRd
sdr_rescue <- function(fit, sdr, H = NULL) {
  sdr <- sdr_rescue_hessian(fit, sdr, H)
  sep <- if (is.environment(fit$cache)) fit$cache$se_sep_named
  if (length(sep)) sdr <- sdr_sep_lost(fit, sdr, sep)
  sdr
}

#' @noRd
sdr_rescue_hessian <- function(fit, sdr, H = NULL) {
  V0 <- sdr$cov.fixed
  if (!length(V0)) return(sdr)
  an <- if (is.environment(fit$cache)) fit$cache$se_analysis
  # an analysis asked for on a fit that did not converge (by
  # check_convergence() or nl_flat_message()) repairs it only when it
  # made the fit a boundary stop (se_boundary_stop()); otherwise the fit
  # keeps sdreport()'s own answer, as below
  if (!is.null(an) && !isTRUE(fit$opt$convergence == 0) &&
        !isTRUE(fit$cache$se_boundary_stop)) {
    an <- NULL
  }
  if (is.null(an)) {
    # Without random effects the fit-time check read the exact Hessian
    # and found nothing, and the finite-difference one sdreport() got
    # must not overrule it. With random effects an inverse that is
    # usable can still rest on a row of noise (cov_from_hessian()), so
    # it is looked at unless there is no Hessian to look at.
    hc <- if (is.environment(fit$cache)) fit$cache$hessian_fixed
    if (cov_usable(V0) &&
          (!length(fit$obj$env$random) || is.null(H) ||
             (isTRUE(hc$clean) && identical(H, hc$H)))) {
      return(sdr)
    }
    # A fit that did not converge keeps sdreport()'s own answer: its
    # curvature is not that of an optimum, and its convergence warning
    # has already said what is wrong with it
    if (!isTRUE(fit$opt$convergence == 0) ||
          identical(fit$cache$se_explained, "convergence")) {
      return(sdr)
    }
    same <- !is.null(H) && identical(H, hc$H)
    E <- if (same) hc$E
    if (is.null(H)) {
      h <- tryCatch(fit_outer_hessian(fit), error = function(e) NULL)
      H <- h$H
      E <- h$E
    }
    if (is.null(H)) return(sdr)
    an <- cov_from_hessian(fit, H, E)
    if (an$tier == 1L) {
      # the joint-precision sdreport of the same fit need not ask again
      if (same) fit$cache$hessian_fixed$clean <- TRUE
      return(sdr)
    }
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
    # under "stop" the check does not wait: a deferred check stopped at
    # the first summary() instead of at frm() (setier review, m10)
    an <- tryCatch(se_fit_analysis(fit, force = identical(act, "stop")),
                   error = function(e) e)
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
#' (`se_analysis`). `force` builds the Hessian even where
#' se_check_at_fit() would wait. The result is kept (`se_fit_an`), so
#' that check_convergence() and nl_flat_message(), which ask before
#' se_check() does, do not make it build twice.
#'
#' @noRd
se_fit_analysis <- function(fit, force = FALSE) {
  if (is.environment(fit$cache) && !is.null(fit$cache$se_fit_an)) {
    return(fit$cache$se_fit_an)
  }
  an <- se_fit_analysis0(fit, force)
  if (!is.null(an) && is.environment(fit$cache)) fit$cache$se_fit_an <- an
  an
}

#' @noRd
se_fit_analysis0 <- function(fit, force) {
  obj <- fit$obj
  p <- fit$opt$par
  if (!length(obj$env$random)) {
    saved <- obj_state_save(obj)
    Hx <- tryCatch(obj$he(p), error = function(e) NULL)
    obj_state_restore(obj, saved)
    if (is.matrix(Hx) && all(dim(Hx) == length(p)) && all(is.finite(Hx))) {
      Hs <- (Hx + t(Hx)) / 2
      an <- cov_from_hessian(fit, Hs, abs(Hx - t(Hx)) / 2, exact = TRUE)
      if (an$tier > 1L) {
        an$lost <- se_relabel(fit, an$lost)
        fit$cache$se_analysis <- an
      }
      return(an)
    }
  }
  if (!force && !se_check_at_fit(fit)) return(NULL)
  h <- fit_outer_hessian(fit)
  an <- cov_from_hessian(fit, h$H, h$E)
  # sdr_rescue() then keeps sdreport()'s inverse of this same Hessian
  # without asking again (tier 3's line probes evaluate the objective)
  h$clean <- an$tier == 1L
  fit$cache$hessian_fixed <- h
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
  bnd <- lost == "boundary"
  # a smooth, gp() or hsgp() block at zero is a term reduced to its
  # unpenalized part, which mgcv reports in silence and which most GAM
  # fits end with; only a block over grouping levels is lme4's singular
  # fit, the same split check_re_structure() makes
  said <- bnd & se_par_grouped(fit, names(lost))
  if (any(said)) se_boundary_act(act, se_boundary_message(fit, lost[said]))
  if (any(!bnd)) re_check_act(act, se_lost_message(fit, lost[!bnd]))
  invisible(NULL)
}

#' The term labels of the random-effect blocks with a parameter at the
#' edge of its parameter space (`se_lost` "boundary").
#'
#' @noRd
se_boundary_blocks <- function(fit) {
  lost <- tryCatch(sdr_of(fit)$se_lost, error = function(e) NULL)
  bn <- names(lost)[lost == "boundary"]
  if (!length(bn)) return(character(0))
  om <- outer_par_map(fit)
  thn <- om$names[om$comp == "theta"]
  out <- character(0)
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    if (any(thn[bk[["theta_idx"]]] %in% bn)) {
      out <- c(out, bk[["term_label"]])
    }
  }
  out
}

#' Does each named outer parameter belong to a random-effect block over
#' grouping levels (rather than a smooth, gp() or hsgp() block)?
#'
#' @noRd
se_par_grouped <- function(fit, nms) {
  om <- outer_par_map(fit)
  thn <- om$names[om$comp == "theta"]
  out <- logical(length(nms))
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    if (is.null(bk[["levels"]])) next
    out[nms %in% thn[bk[["theta_idx"]]]] <- TRUE
  }
  out
}

#' A variance component at its boundary is reported at the level lme4
#' reports a singular fit: a message, not a warning. It is an ordinary
#' outcome of a correct fit (lme4 tells it apart from a convergence
#' problem for that reason): on 600 simulated random-intercept and
#' binomial fits it fired on 185, every one a fit lme4 calls singular,
#' so as a warning it would reach a third of such fits
#' (dev/setier-findings.md). The `check_se` setting still decides:
#' "ignore" is silent and "stop" stops.
#'
#' @noRd
se_boundary_act <- function(what, msg) {
  switch(what %||% "warning",
         ignore = invisible(NULL),
         stop = frm_stop(msg, call. = FALSE),
         frm_message(msg, class = "frmtmb_boundary_fit"))
}

#' The covariance parameters at the edge of their parameter space: a
#' parameter of a random-effect block that the standard-error check found
#' flat, and where the likelihood does not change as it moves further
#' toward its own end (se_at_edge()), toward zero for a log standard
#' deviation and toward the nearer infinity for the others (a
#' correlation or mixing parameter at its limit, a range run off). The
#' other flat parameters of a block whose standard deviation is at zero
#' join them, since a correlation or a range is undefined there. `lost`
#' is the named reasons of se_tier3().
#'
#' Asked of the likelihood and not of the estimate's size, because the
#' size is in the response's units: a cut at an sd of 1e-4, which
#' diagnose() reads, called a binomial group sd of 1.35e-4 that lme4
#' reports as singular an unidentified parameter (dev/setier-singular.R,
#' bin15 seed 57). A parameter that trades off against another (an
#' observation-level sd against sigma) is flat too, but moving it alone
#' costs likelihood, so it keeps its warning.
#'
#' Each name's kind ("zero", "limit" or "undefined") is kept in
#' `fit$cache$se_boundary_kind` for the message.
#'
#' @noRd
se_boundary_names <- function(fit, lost) {
  # Under quadrature the boundary is not asked: the upward check
  # (se_sd_gain_up()) cannot run on that objective, and without it a fit
  # stuck near sd = 0 is called a boundary (4 of 4 rescaled-response fits
  # 0.0007 to 0.825 below lme4, 20 of 20 GLMMs started at exp(-12),
  # setier review, final check RQ1); such an sd keeps base's "flat"
  # verdict and its warning
  if (isTRUE(fit$quadrature)) return(character(0))
  th <- fit$estimates[["theta"]]
  flat <- names(lost)[lost == "flat"]
  if (!length(th) || !length(flat)) return(character(0))
  om <- outer_par_map(fit)
  thn <- om$names[om$comp == "theta"]
  if (length(thn) != length(th)) return(character(0))
  sd_i <- tryCatch(log_sd_theta_index(fit), error = function(e) integer(0))
  kind <- character(0)
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    ti <- bk[["theta_idx"]]
    ti <- ti[!is.na(ti) & ti >= 1L & ti <= length(th)]
    cand <- ti[thn[ti] %in% flat]
    if (!length(cand)) next
    dir <- ifelse(cand %in% sd_i, -1, sign(th[cand]))
    edge <- cand[vapply(seq_along(cand), function(i) {
      dir[i] != 0 && se_at_edge(fit, thn[cand[i]], dir[i])
    }, NA)]
    # an sd the optimizer left near zero short of an interior maximum is
    # flat toward zero too; a larger value that gains tells the two apart
    short <- edge[edge %in% sd_i][vapply(edge[edge %in% sd_i], function(j) {
      se_sd_gain_up(fit, j) > se_edge_tol
    }, NA)]
    edge <- setdiff(edge, short)
    k <- stats::setNames(ifelse(edge %in% sd_i, "zero", "limit"),
                         thn[edge])
    if (any(edge %in% sd_i)) {
      # a correlation or a range of a block whose variance is zero is
      # not at a limit, even where moving it changes nothing: it is
      # undefined there
      und <- setdiff(cand, edge[edge %in% sd_i])
      k <- k[edge %in% sd_i]
      k <- c(k, stats::setNames(rep("undefined", length(und)), thn[und]))
    } else if (is.null(bk[["levels"]]) && length(ti) >= 2L &&
                 length(cand) == length(ti) &&
                 se_block_ridge_edge(fit, thn[ti])) {
      # a gp() or hsgp() term none of whose hyperparameters the data
      # determine, flat along a joint move (an hsgp length scale at
      # exp(-11.9) where only sd^2 * length scale matters,
      # dev/setier-gpprobe.R), which no single coordinate's probe sees;
      # the curve keeps its uncertainty given them, as mgcv's does. A
      # one-hyperparameter block (an s() term) must pass the probe
      # itself: s(x) + s(x2) with x2 = x has both sds at log 0.42, flat
      # only jointly, and a move of either by -2 costs 0.99
      # (dev/setier-rev-smooth2.R)
      k <- stats::setNames(rep("limit", length(cand)), thn[cand])
    }
    kind <- c(kind, k, stats::setNames(rep("short", length(short)),
                                       thn[short]))
  }
  kind <- kind[!duplicated(names(kind))]
  if (is.environment(fit$cache)) {
    fit$cache$se_boundary_kind <- c(fit$cache$se_boundary_kind, kind)
  }
  names(kind)[kind != "short"]
}

#' How much log-likelihood a larger value of a log standard deviation
#' gains: the best of the sd at 0.01, 0.03, 0.1, 0.3 and 1 times the
#' scale of its term (the residual sd where the response has one, else
#' 1 on the link scale), where that is above the estimate, the rest
#' held; 0 when none gains. At a boundary every interior value is worse,
#' so the callers refuse the boundary on any gain above `se_edge_tol`.
#'
#' A fit can stop far short of an interior maximum with an sd near zero,
#' where the log-likelihood is flat toward zero, so the edge test alone
#' calls it a boundary: on y ~ x + (1 | g) with the response in units of
#' 1e-3, seed 36, nlminb stopped with "false convergence (8)" at an sd
#' of 1.7e-12 times sigma, 0.825 below lme4's log-likelihood at 0.33
#' (setier review, re-check RB2; dev/setier-rev2-trap.R). With the
#' response in units of 1e3, seed 23, the fit stopped 7e-4 below lme4,
#' whose sd is 0.053 times sigma; held at the fit's other values the
#' log-likelihood gains 4.3e-4 at 0.03 sigma and loses at 0.1 sigma,
#' which is why the grid starts low and the threshold is not `grad_tol`
#' (dev/setier-trapdbg.R).
#'
#' @noRd
se_sd_gain_up <- function(fit, k) {
  th <- fit$estimates[["theta"]]
  om <- outer_par_map(fit)
  thpos <- which(om$comp == "theta")
  p <- fit$opt$par
  if (length(thpos) != length(th) || length(om$names) != length(p)) {
    return(0)
  }
  scale <- 1
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    if (!k %in% bk[["theta_idx"]]) next
    key <- bk[["components"]][[1L]]$lp_key
    lp <- fit$frame[["linpreds"]][[key]]
    dp <- tryCatch(eval_dpars(fit, b = NULL)[[lp[["resp"]]]],
                   error = function(e) NULL)
    sg <- dp[["sigma"]]
    if (is.numeric(sg) && length(sg) && all(is.finite(sg))) {
      scale <- stats::median(sg)
    }
    break
  }
  up <- log(c(0.01, 0.03, 0.1, 0.3, 1) * scale)
  up <- up[up > p[thpos[k]]]
  # a quadrature objective is not accurate above the sd it was fitted
  # at: on test-quadrature-defects.R's nested beta fit, an sd at
  # exp(-11.2) moved up by 2 lowered it by 0.009 and moved to 0.01 by
  # 17, where the Laplace objective rose by 0.0015 (dev/setier-quadup.R),
  # so the question is not asked there
  if (isTRUE(fit$quadrature)) return(0)
  if (!length(up)) return(0)
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
  if (!is.finite(f0)) return(0)
  gain <- vapply(up, function(v) {
    q <- p
    q[thpos[k]] <- v
    f0 - tryCatch(obj$fn(q), error = function(e) NA_real_)
  }, 0)
  gain <- gain[is.finite(gain)]
  if (length(gain)) max(0, gain) else 0
}

#' Is a block of hyperparameters flat along its own flattest direction?
#' The block's rows of the outer Hessian the check kept
#' (`hessian_fixed`), scaled to unit diagonal; a move along the smallest
#' eigenvector, scaled so that its largest component is `se_edge_step`,
#' leaves the likelihood within `grad_tol` on at least one side.
#'
#' @noRd
se_block_ridge_edge <- function(fit, names) {
  H <- if (is.environment(fit$cache)) fit$cache$hessian_fixed$H
  nm <- outer_par_names(fit)
  j <- match(names, nm)
  p <- fit$opt$par
  if (is.null(H) || anyNA(j) || NROW(H) != length(p)) return(FALSE)
  Hb <- H[j, j, drop = FALSE]
  D <- sqrt(abs(diag(Hb)))
  if (!all(is.finite(Hb)) || !all(D > 0)) return(FALSE)
  v <- eigen(Hb / outer(D, D), symmetric = TRUE)$vectors
  v <- v[, ncol(v)] / D
  v <- v * se_edge_step / max(abs(v))
  tol <- fit$control$grad_tol %||% 1e-3
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  f <- function(x) tryCatch(obj$fn(x), error = function(e) NA_real_)
  f0 <- f(p)
  ch <- vapply(c(1, -1), function(s) {
    q <- p
    q[j] <- q[j] + s * v
    abs(f(q) - f0)
  }, 0)
  is.finite(f0) && any(is.finite(ch) & ch <= tol)
}

#' Is the estimate at the end of this parameter's range, as far as the
#' likelihood can tell? One objective evaluation with the parameter
#' moved by `se_edge_step` in direction `dir` (a factor of 0.14 on a
#' standard deviation), the rest held: at an end, moving further changes
#' the log-likelihood by no more than `se_edge_tol`.
#'
#' Not within `grad_tol`: an interior estimate with little curvature
#' also changes by less than that. Asked of every log sd (se_edge_sd()),
#' a `grad_tol` test called 5 of 331 fits lme4 does not call singular
#' boundary fits (sds of 0.014 to 0.074, dev/setier-singular.R), whose
#' log-likelihood falls by 1e-5 or more; at a boundary it moves by
#' about 1e-8 (ri20 seed 30, dev/setier-rev-miss.R). Two-sided, so that
#' a parameter whose log-likelihood still rises toward its end (a
#' correlation beside an sd at exp(-6.6), rising 9.3e-5 over a move of 2,
#' dev/setier-grby.R) is not called settled there.
#'
#' @noRd
se_at_edge <- function(fit, name, dir) {
  p <- fit$opt$par
  j <- match(name, outer_par_names(fit))
  if (is.na(j) || length(p) != length(outer_par_names(fit))) return(FALSE)
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  q <- p
  q[j] <- q[j] + dir * se_edge_step
  f0 <- tryCatch(obj$fn(p), error = function(e) NA_real_)
  f1 <- tryCatch(obj$fn(q), error = function(e) NA_real_)
  is.finite(f0) && is.finite(f1) && abs(f1 - f0) <= se_edge_tol
}

#' The rounding allowed in se_at_edge(): a Laplace objective's inner
#' solve is accurate to far below this, and a quadrature one moved by 2
#' changed by 3e-9 (dev/setier-quadprobe.R).
#'
#' @noRd
se_edge_tol <- 1e-6

#' How far se_at_edge() moves a parameter on its internal scale.
#'
#' A step of 2 on the log scale tells a standard deviation the data
#' ignore (its row is noise, so the change is below 1e-6) from one that
#' trades off against another parameter (an observation-level sd moved
#' alone costs whole units of log-likelihood). A step of 20 left the
#' range where a quadrature objective is accurate: on the nested beta
#' fit of test-quadrature-defects.R an sd at exp(-11.2) moved by -5
#' changed it by 1.8, by -20 gave NaN, and by -2 changed it by 3e-9
#' (dev/setier-quadprobe.R).
#'
#' @noRd
se_edge_step <- 2

#' The message for variance components at their boundary.
#'
#' @noRd
se_boundary_message <- function(fit, lost) {
  nms <- names(lost)
  kind <- fit$cache$se_boundary_kind[nms]
  kind[is.na(kind)] <- "limit"
  lab <- se_boundary_labels(fit, nms)
  part <- function(k, one, more) {
    n <- lab[kind == k]
    if (!length(n)) return(NULL)
    paste(paste(n, collapse = ", "),
          if (length(n) == 1L) paste("is", one) else paste("are", more))
  }
  parts <- c(part("zero", "at zero", "at zero"),
             part("limit", "at the limit of its range",
                  "at the limits of their ranges"),
             part("undefined",
                  "undefined there, since a standard deviation of its term
                   is zero",
                  "undefined there, since a standard deviation of their
                   terms is zero"))
  parts <- gsub("\\s+", " ", parts)
  one <- length(lost) == 1L
  paste0("Boundary (singular) fit: ", paste(parts, collapse = "; "),
         ". The likelihood does not change as ",
         if (one) "it moves" else "they move", " further that way, so ",
         if (one) "it has" else "they have",
         " no standard error (NaN in summary(), vcov() and confint()); ",
         "the other standard errors are kept. lme4 reports such a fit as ",
         "singular; diagnose() lists the components at the boundary")
}

#' brms's name for a covariance parameter where it has one of its own
#' (`sd_g__Intercept`, `cor_g__Intercept__x`, par_alias_index()), and
#' otherwise `theta_k (sd of <term> <column>)` or `theta_k (in <term>)`.
#'
#' @noRd
se_boundary_labels <- function(fit, nms) {
  out <- flat_par_display(fit, nms)
  om <- outer_par_map(fit)
  thn <- om$names[om$comp == "theta"]
  for (bk in fit$frame[["re_blocks"]] %||% list()) {
    hit <- out %in% thn[bk[["theta_idx"]]]
    out[hit] <- paste0(out[hit], " (in ", bk[["term_label"]], ")")
  }
  ai <- tryCatch(par_alias_index(fit), error = function(e) integer(0))
  pos <- match(nms, om$names)
  for (i in seq_along(nms)) {
    a <- names(ai)[ai == pos[i]]
    if (length(a) && !is.na(pos[i])) out[i] <- a[1L]
  }
  out
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
#' - the separation warning of check_convergence() explains the
#'   coefficients it names (`pars`);
#' - the convergence verdict explains the whole fit, and se_check()
#'   returns before asking.
#'
#' @noRd
se_explained_pars <- function(fit) {
  fam <- !is.null(fit$cache$se_explained)
  # the coefficients check_convergence() named in its separation warning
  sep <- fit$cache$se_explained_sep %||% character(0)
  ent <- attr(fit$frame, "se_explained")
  if (!length(ent)) return(list(family = fam, pars = sep))
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
  list(family = fam, pars = c(out[!is.na(out)], sep))
}

#' A flat mean coefficient of a binomial-type fit that has run far out
#' on the link, in a predictor separation_check() found separated, is
#' separation, which says what happened better than "flat" does (lane
#' nanse review, m7). A flat variance component at
#' zero is a boundary fit (se_boundary_names()).
#'
#' @noRd
se_relabel <- function(fit, lost) {
  if (!length(lost)) return(lost)
  lost[names(lost) %in% se_boundary_names(fit, lost) &
         lost == "flat"] <- "boundary"
  kind <- fit$cache$se_boundary_kind
  lost[names(lost) %in% names(kind)[kind == "short"] &
         lost == "flat"] <- "short"
  if (fit$REML || isTRUE(fit$control$profile)) return(lost)
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
  # only where check_convergence()'s certificate proved separation in
  # that predictor: a flat coefficient above 10 on unseparated data (a
  # bernoulli polynomial) is not one (setier review, re-check m-c)
  sep <- intersect(sep, fit$cache$se_explained_sep %||% character(0))
  hit <- names(lost) %in% sep & lost == "flat"
  lost[hit] <- "separation"
  lost
}

#' Separation in a binomial-type mean, proved rather than guessed.
#'
#' The maximum likelihood of a binomial-type mean does not exist when
#' some combination `d` of the coefficients moves every observation it
#' moves toward that observation's own outcome (complete separation when
#' it moves all of them, quasi-complete otherwise; Albert and Anderson
#' 1984). The likelihood then rises without end along `d`, the optimizer
#' stops wherever its limits do, and that stop is platform-dependent: on
#' `yb ~ z + x` with `yb = z > 0` nlminb stopped at code 0 after 1939
#' evaluations with the reference BLAS and at code 9 ("function
#' evaluation limit reached") after 1997 with OpenBLAS, which then said
#' nothing about separation (dev/reviews/2026-10-06-cifix.md, m5).
#'
#' The fit itself supplies `d`. The observations whose fitted
#' probability is still away from both ends hold the estimate in place,
#' so the runaway direction lies in the null space of their rows of the
#' design; the projection of the estimate onto that null space is the
#' candidate. It is a certificate, checked exactly as the definition
#' reads: no observation it moves against its outcome (an observation
#' with both outcomes among its trials must not move), at least one
#' moved. A design without such a direction is never named, whatever the
#' size of its estimates, so a large but finite effect is not a false
#' alarm. brglm2's detect_separation() decides the same question with a
#' linear program; this asks it only along the direction the fit took.
#'
#' Returns `NULL` or `list(message, pars)`, `pars` the outer names of the
#' coefficients `d` loads on.
#'
#' @noRd
separation_check <- function(fit) {
  if (!length(fit$opt$par) || fit$REML || isTRUE(fit$control$profile)) {
    return(NULL)
  }
  om <- outer_par_map(fit)
  pb <- which(om$comp == "beta")
  hits <- list()
  for (lp in fit$frame[["linpreds"]]) {
    rspec <- fit$spec$responses[[lp[["resp"]]]]
    fam <- rspec$family
    if (!identical(lp[["par"]], "beta") || !is.null(lp[["nl_body"]]) ||
          !fam[["family"]] %in% c("binomial", "bernoulli", "beta_binomial") ||
          !lp[["dpar"]] %in% (fam[["primary_dpars"]] %||% "mu") ||
          !(lp[["link"]]$name %||% "") %in% sep_links ||
          is.null(lp[["X"]]) || !ncol(lp[["X"]])) {
      next
    }
    h <- tryCatch(sep_certificate(fit, lp, rspec), error = function(e) NULL)
    if (is.null(h)) next
    # every coefficient of this predictor: the separated observations
    # carry no curvature, so a standard error any of them loses is lost
    # to the separation, named or not (on complete separation the
    # Hessian of the whole predictor vanishes)
    pos <- pb[lp[["idx"]]]
    h$pars <- om$names[pos[!is.na(pos)]]
    h$named <- om$names[pos[h$cols][!is.na(pos[h$cols])]]
    # a few names and the count: a 2050-level factor with 339 separated
    # levels listed every one of them
    cn <- colnames(lp[["X"]])[h$cols]
    h$label <- paste0(coef_block_key(fit, lp), ": ",
                      paste(utils::head(cn, sep_name_max), collapse = ", "),
                      if (length(cn) > sep_name_max) {
                        paste0(" and ", length(cn) - sep_name_max,
                               " more (", length(cn), " coefficients)")
                      })
    hits[[length(hits) + 1L]] <- h
  }
  if (!length(hits)) return(NULL)
  complete <- all(vapply(hits, `[[`, NA, "complete"))
  moved <- sum(vapply(hits, `[[`, 0L, "moved"))
  msg <- paste0(
    "The data separate the outcomes (", if (complete) "complete" else
      "quasi-complete", " separation): along a combination of ",
    paste(vapply(hits, `[[`, "", "label"), collapse = "; "), " every ",
    "observation that moves (", moved, ") moves toward its own outcome, ",
    "so the likelihood rises without end and the maximum likelihood ",
    "estimates do not exist. The reported values are where the optimizer ",
    "stopped",
    if (!isTRUE(fit$opt$convergence == 0)) {
      paste0(" (", fit$opt$message, ")")
    },
    ", and these coefficients have no standard error. Remove or merge ",
    "the separating predictor, or set a prior on it (set_prior()) to ",
    "keep it finite")
  list(message = msg, pars = unique(unlist(lapply(hits, `[[`, "pars"))),
       named = unique(unlist(lapply(hits, `[[`, "named"))))
}

#' How many coefficients the separation warning names before it gives
#' the count.
#'
#' @noRd
sep_name_max <- 5L

#' Links whose inverse reaches 0 and 1 only in the limit, so that
#' separation sends the linear predictor to infinity.
#'
#' @noRd
sep_links <- c("logit", "probit", "cloglog", "cauchit", "loglog")

#' The certificate of separation_check() for one linear predictor, or
#' `NULL`: `list(cols, moved, complete)`.
#'
#' The design stays as the frame holds it, sparse or dense: an earlier
#' version made it dense and took a singular value decomposition of it,
#' 54 to 68 s and 2.3 GB on a 500-level factor at n = 1e5 with
#' `sparse_x = TRUE`, where the fit itself took 96 s and 297 MB
#' (dev/setier-rev-sepcost.R). The candidates are, in order:
#'
#' - the estimate's own direction, which needs no decomposition and is
#'   what a complete separation that stopped short shows (a probit fit
#'   whose underflow stopped it at a slope of 19.1 kept nine
#'   observations above a tail of 1e-3, so the holding rows had full
#'   rank; dev/setier-rev-probit.R);
#' - the null space of the holding rows, from the eigenvectors of their
#'   p x p cross product at about zero, scaled to unit columns. Squaring
#'   the conditioning does no harm: every candidate is verified on every
#'   row, so a near-null direction can only yield a true certificate. It
#'   is skipped when no observation has left a tail of 1e-3 (no row
#'   outside the holding set); above `sep_p_max` columns single columns
#'   are tried instead (sep_column_cert()).
#'
#' @noRd
sep_certificate <- function(fit, lp, rspec) {
  rn <- rspec$resp_name %||% lp[["resp"]]
  y <- fit$frame[["y"]][[rn]]
  X <- patch_mo_cols(fit, lp, lp[["X"]])
  n <- nrow(X)
  if (!is.numeric(y) || is.matrix(y) || length(y) != n) return(NULL)
  av <- fit$frame[["aterm_values"]][[rn]] %||% list()
  size <- rep_len(av[["trials"]] %||% 1, n)
  w <- rep_len(av[["weights"]] %||% 1, n)
  use <- is.finite(y) & size > 0 & w > 0
  if (!any(use)) return(NULL)
  beta <- fit$estimates[["beta"]][lp[["idx"]]]
  xb <- as.numeric(X %*% beta)
  eta <- xb
  b <- fit$estimates[["b"]]
  if (!is.null(lp[["Z"]]) && length(b)) {
    eta <- eta + as.numeric(lp[["Z"]] %*% expand_b(fit$frame, b,
                                                   fit$estimates[["theta"]]))
  }
  if (!is.null(lp[["offset"]])) eta <- eta + lp[["offset"]]
  mu <- lp[["link"]]$linkinv(eta)
  # +1 every trial a success, -1 every trial a failure, 0 both occur
  s <- ifelse(y >= size, 1, ifelse(y <= 0, -1, 0))[use]
  tail <- ifelse(s > 0, 1 - mu[use], ifelse(s < 0, mu[use], 1))
  Xu <- if (all(use)) X else X[use, , drop = FALSE]
  sc <- sqrt(as.numeric(Matrix::colSums(Xu^2)))
  sc[!is.finite(sc) | sc == 0] <- 1
  # the largest entry of the unit-column design, without forming it
  xs_max <- if (inherits(Xu, "sparseMatrix")) {
    Xd <- methods::as(Xu, "CsparseMatrix")
    j <- rep.int(seq_len(ncol(Xd)), diff(Xd@p))
    if (length(Xd@x)) max(abs(Xd@x) / sc[j]) else 0
  } else {
    max(vapply(seq_len(ncol(Xu)), function(k) max(abs(Xu[, k])) / sc[k], 0))
  }
  verify <- function(d_s, m = NULL) {
    if (!all(is.finite(d_s)) || !any(d_s != 0)) return(NULL)
    nrm <- sqrt(sum(d_s^2))
    d_s <- d_s / nrm
    m <- if (is.null(m)) as.numeric(Xu %*% (d_s / sc)) else m / nrm
    # a direction the whole design hardly moves (near-collinear
    # columns) moves observations by rounding, whose signs prove nothing
    if (max(abs(m)) <= 1e-6 * xs_max) return(NULL)
    tol <- 1e-8 * max(abs(m))
    mv <- abs(m) > tol
    if (!any(mv) || any(s[mv] == 0) || any(s[mv] * m[mv] < 0)) return(NULL)
    dd <- abs(d_s)
    list(cols = which(dd > 1e-6 * max(dd)), moved = sum(mv),
         complete = all(mv))
  }
  bs <- beta * sc
  out <- verify(bs, if (all(use)) xb else xb[use])
  if (!is.null(out)) return(out)
  if (all(s == 0 | tail > 1e-3)) return(NULL)
  if (ncol(Xu) > sep_p_max) return(sep_column_cert(Xu, s))
  # tighter thresholds keep more observations in the holding set; the
  # certificate is exact at any of them, so trying several only widens
  # what is found
  for (tau in c(1e-3, 1e-6, 1e-9)) {
    hold <- s == 0 | tail > tau
    if (all(hold)) next
    N <- sep_null_space(Xu[hold, , drop = FALSE], sc)
    if (!NCOL(N)) next
    cand <- cbind(N %*% crossprod(N, bs), N, -N)
    for (k in seq_len(NCOL(cand))) {
      out <- verify(cand[, k])
      if (!is.null(out)) return(out)
    }
  }
  NULL
}

#' Single columns as certificates, for a design too wide to decompose:
#' a column whose nonzero rows all move toward their own outcome under
#' +e_j or -e_j (a level of a factor whose rows are all 0, or all 1) is
#' separation by the definition, checked in O(nnz). On a 2050-level
#' factor with `sparse_x = TRUE` and 39 levels of all 0 nothing was named
#' (setier review, re-check m-b). Every certified column is named.
#'
#' @noRd
sep_column_cert <- function(X, s) {
  if (inherits(X, "sparseMatrix")) {
    Xd <- methods::as(X, "CsparseMatrix")
    j <- rep.int(seq_len(ncol(Xd)), diff(Xd@p))
    v <- Xd@x * s[Xd@i + 1L]
    mixed <- s[Xd@i + 1L] == 0
  } else {
    nz <- which(X != 0)
    j <- (nz - 1L) %/% nrow(X) + 1L
    i <- (nz - 1L) %% nrow(X) + 1L
    v <- X[nz] * s[i]
    mixed <- s[i] == 0
  }
  keep <- v != 0 | mixed
  j <- j[keep]
  v <- v[keep]
  mixed <- mixed[keep]
  if (!length(j)) return(NULL)
  bad <- tapply(mixed, j, any)
  pos <- tapply(v > 0, j, all) & !bad
  neg <- tapply(v < 0, j, all) & !bad
  cols <- as.integer(names(bad))[pos | neg]
  if (!length(cols)) return(NULL)
  moved <- sum(tabulate(j, ncol(X))[cols])
  list(cols = cols, moved = moved, complete = FALSE)
}

#' The most columns sep_certificate() decomposes: its cross product is
#' p x p and dense, and its eigen decomposition grows as p^3.
#'
#' @noRd
sep_p_max <- 2000L

#' An orthonormal basis of the null space of `A`'s columns scaled by
#' `sc` (`A diag(1 / sc) v = 0`), from the eigenvectors of the p x p
#' cross product, which a sparse `A` gives without being made dense.
#' The cut, 1e-10 of the largest eigenvalue (a singular value of 1e-5
#' of the largest), errs toward too many directions, which the caller's
#' verification on every row makes harmless.
#'
#' @noRd
sep_null_space <- function(A, sc) {
  p <- ncol(A)
  if (!nrow(A)) return(diag(1, p))
  G <- as.matrix(Matrix::crossprod(A)) / outer(sc, sc)
  e <- eigen((G + t(G)) / 2, symmetric = TRUE)
  keep <- e$values <= 1e-10 * max(e$values, 0)
  e$vectors[, keep, drop = FALSE]
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
      short = paste0("the fit stopped short of the maximum: ",
                     if (one) "this standard deviation is" else
                       "these standard deviations are",
                     " near zero, where the likelihood is flat, but a ",
                     "larger value raises the log-likelihood, so ",
                     if (one) "it is" else "they are",
                     " not at a boundary. Refit from a larger starting ",
                     "value (start =; see par_template())"),
      boundary = paste0("at the edge of the parameter space (a variance ",
                        "at zero, or a correlation, mixing or range ",
                        "parameter at its limit), where the likelihood no ",
                        "longer changes: a boundary, or singular, fit"),
      separation = paste0("the data separate the outcomes, so ",
                          if (one) "the estimate runs" else
                            "the estimates run",
                          " off toward infinity and the optimizer stopped ",
                          "where its tolerances did"),
      flat = paste0("the likelihood is flat along ", it, " at the ",
                    "estimates, so the data do not determine ", it,
                    " there (a standard deviation at 0 does this, and so do ",
                    "a mo() simplex split by a category no row is in and ",
                    "parameters that enter only through a combination)",
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
  for (r in c("bound", "boundary", "separation", "short", "flat",
              "concave", "nonfinite")) {
    n <- names(lost)[lost == r]
    if (!length(n)) next
    lab <- if (r == "boundary") se_boundary_labels(fit, n) else {
      vapply(n, label, "")
    }
    out <- c(out, paste0(paste(lab, collapse = ", "), ": ",
                         text(r, length(n) == 1L)))
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
#' it does not depend on how a coefficient is scaled; a prediction that
#' does not read a lost simplex coordinate has an exactly zero component
#' and keeps its band. This is predict()'s `alias_null` estimability
#' test, for the directions the Hessian lost rather than the design.
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
#' direction. Without random effects the directions come from an exact
#' AD Hessian, so they are accurate to rounding and the test can be
#' tight. With random effects (joint_cov_repair()) they come from a
#' finite-difference Hessian, but se_tier3() zeroes every loading of a
#' parameter that keeps its standard error, so a prediction that loads
#' no lost parameter has a component of exactly zero and the threshold
#' does not decide it.
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

#' Is a fit the optimizer stopped short at a boundary, as lme4 calls a
#' singular fit? The largest gradient (`g`, in the optimizer's units) is
#' within `grad_tol`, and every standard error the check loses there is
#' a boundary one. The analysis is built even where se_check_at_fit()
#' would wait, since its answer decides which report the fit gets; under
#' `check_se = "ignore"` it is not asked, and the convergence warning
#' stays.
#'
#' @noRd
se_boundary_stop <- function(fit, g, control) {
  # (never under quadrature, where no boundary verdict is given:
  # se_boundary_names())
  if (isTRUE(fit$quadrature) ||
        !is.finite(g) || g > (control$grad_tol %||% 1e-3) ||
        identical(control$check_se, "ignore") ||
        !is.environment(fit$cache)) {
    return(FALSE)
  }
  an <- tryCatch(se_fit_analysis(fit, force = TRUE),
                 error = function(e) NULL)
  lost <- an$lost
  ok <- length(lost) > 0L && all(lost == "boundary")
  fit$cache$se_boundary_stop <- ok
  # an sdreport taken at the fit (se = TRUE) predates the analysis and
  # kept its own covariance, as for any fit that stopped short
  if (ok && !is.null(fit$cache$sdr)) {
    fit$cache$sdr <- tryCatch(autoscale_sdreport(fit),
                              error = function(e) fit$cache$sdr)
  }
  ok
}
