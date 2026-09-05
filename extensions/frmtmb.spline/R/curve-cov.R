# The joint covariance of a grid prediction, and why it is assembled
# here rather than asked for.
#
# A simultaneous band, a delta-method derivative and an implicit-function
# feature all need the same object: the covariance of the WHOLE grid
# prediction, Sigma = C V C', where C maps the coefficient vector to the
# curve and V is the joint covariance of the fixed and random
# coefficients. A smooth's wiggly part is a random-effect block even when
# the smooth is a population term, so V must carry the b block. vcov(fit,
# full = TRUE) does not carry it under either of its branches, and
# predict(se.fit = TRUE) forms Sigma internally and returns only
# sqrt(diag(Sigma)). See dev/spline-seam-proposal.md.
#
# So C is rebuilt and V is re-derived, and the result is CHECKED against
# the one exported route that overlaps: diag(C V C') must equal
# predict(se.fit = TRUE)^2. The check is not decoration. It is what
# licenses the rest of the package to use a covariance that no exported
# function returned.

#' One prediction on the link scale, as a plain numeric vector.
#'
#' @noRd
sp_predict_eta <- function(fit, newdata, dpar, resp, re.form) {
  as.numeric(stats::predict(fit, newdata = newdata, type = "link",
                            dpar = dpar, resp = resp, re.form = re.form))
}

#' Is the prediction linear in the coefficients?
#'
#' The unit difference below is a design column only for a predictor
#' that is linear in `(beta, b)`. Two probes settle it: a linear map
#' satisfies `f(c + 2v) - f(c) = 2 (f(c + v) - f(c))` for every `v`, and
#' a random `v` makes an accidental pass improbable.
#'
#' This probe is not what catches a reduced-rank block, and the
#' distinction is worth stating because it is easy to assume otherwise.
#' An `rr` block's loadings live in `theta`, not in `b`, so `eta` IS
#' linear in `b` at fixed `theta` and this probe passes. What the
#' perturbation cannot see is the derivative with respect to the
#' loadings, which core's own delta method carries separately as
#' `rr_jacobians()`. The covariance check at the end of
#' `sp_curve_parts()` is what catches that: measured, the assembled
#' standard errors come out 27 percent away from
#' `predict(se.fit = TRUE)`'s on an `rr` fit at `re.form = NULL`, and
#' the call refuses. Nothing here is redundant with anything else.
#'
#' @noRd
sp_check_linear <- function(fit, newdata, dpar, resp, re.form, eta0) {
  v <- list(beta = sp_probe(length(fit$estimates[["beta"]]), 9314L),
            betad = sp_probe(length(fit$estimates[["betad"]]), 7717L),
            b = sp_probe(length(fit$estimates[["b"]]), 5501L))
  bump <- function(mult) {
    f2 <- fit
    for (cp in names(v)) {
      if (length(v[[cp]])) {
        f2$estimates[[cp]] <- f2$estimates[[cp]] + mult * v[[cp]]
      }
    }
    sp_predict_eta(f2, newdata, dpar, resp, re.form) - eta0
  }
  d1 <- bump(1)
  d2 <- bump(2)
  max(abs(d2 - 2 * d1)) <= 1e-8 * max(abs(d1), 1)
}

#' A reproducible random probe of a given length.
#'
#' Seeded from a constant so that two calls on one fit agree and a user's
#' stream is left alone: the probe decides which coefficients are live,
#' and that decision must not move between calls.
#'
#' @noRd
sp_probe <- function(n, seed) {
  if (!n) return(numeric(0))
  old <- if (exists(".Random.seed", globalenv())) {
    get(".Random.seed", globalenv())
  }
  set.seed(seed)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()))
  stats::rnorm(n)
}

#' The grid design `C`, and the coefficient positions its columns sit at.
#'
#' `eta` is LINEAR in the coefficient vector, so the difference of two
#' predictions one unit apart is the design column itself, with no
#' truncation error and nothing to tune. That is the whole trick.
#'
#' Coefficients are probed in chunks first. A chunk whose random
#' perturbation leaves the prediction alone contains no contributing
#' coefficient, so the whole chunk is skipped without a call each: under
#' `re.form = NA` every grouping block is such a chunk, and a
#' per-subject model is mostly grouping blocks. Cancellation inside a
#' chunk would have to hit a random vector exactly, which is a measure
#' zero event rather than a case to guard.
#'
#' @noRd
sp_curve_design <- function(fit, newdata, dpar, resp, re.form, eta0,
                            chunk = 24L) {
  cols <- list()
  comp <- character(0)
  idx <- integer(0)
  ncall <- 0L
  scale0 <- max(abs(eta0), 1)
  # betad carries the coefficients of every dpar that is not the
  # location one, so a curve read off dpar = "sigma" or dpar = "gamma1"
  # lives there and not in beta. Its FIXED entries are excluded up
  # front: a dpar held at a constant moves the prediction when it is
  # perturbed and has no row in the joint precision to pair the column
  # with, so including it would build an A wider than V.
  bd_fixed <- fit$frame[["betad_fixed_idx"]] %||% integer(0)
  for (cp in c("beta", "betad", "b")) {
    cv <- fit$estimates[[cp]]
    if (!length(cv)) next
    live <- rep(TRUE, length(cv))
    if (cp == "betad" && length(bd_fixed)) live[bd_fixed] <- FALSE
    if (length(cv) > chunk) {
      pr <- sp_probe(length(cv), 4471L)
      grp <- split(seq_along(cv), ceiling(seq_along(cv) / chunk))
      for (g in grp) {
        f2 <- fit
        f2$estimates[[cp]][g] <- cv[g] + pr[g]
        d <- sp_predict_eta(f2, newdata, dpar, resp, re.form) - eta0
        ncall <- ncall + 1L
        if (max(abs(d)) <= 1e-10 * scale0) live[g] <- FALSE
      }
    }
    for (j in which(live)) {
      f2 <- fit
      f2$estimates[[cp]][j] <- cv[j] + 1
      cols[[length(cols) + 1L]] <-
        sp_predict_eta(f2, newdata, dpar, resp, re.form) - eta0
      ncall <- ncall + 1L
      comp <- c(comp, cp)
      idx <- c(idx, j)
    }
  }
  if (!length(cols)) {
    stop("frm_curve(): the prediction depends on no estimated ",
         "coefficient, so it has no covariance. Every term of the linear ",
         "predictor is either fixed or dropped by `re.form`",
         call. = FALSE)
  }
  C <- do.call(cbind, cols)
  keep <- colSums(abs(C)) > 0
  list(C = C[, keep, drop = FALSE], comp = comp[keep], idx = idx[keep],
       n_predict = ncall)
}

#' The joint covariance of `(beta, b)` at the positions `C`'s columns
#' name.
#'
#' Read from the fit's OWN cache wherever possible. `predict(se.fit =
#' TRUE)` runs `get_joint_cov()` (`R/predict.R:12`), which memoizes the
#' inverted joint precision in `fit$cache$Vjoint`; `sp_curve_parts()`
#' makes that call anyway, for the check, so by the time this function
#' runs the object it needs already exists and is free.
#'
#' Recomputing it instead, which this package did until the review
#' measured the cost, was wrong on three counts. It paid for a second
#' `sdreport()` per call. It then took a Schur complement over the
#' coefficients the curve does NOT touch, on a matrix densified out of
#' the sparse one TMB had produced, which was the whole scaling wall:
#' 114 s and 2.1 GB at 8000 random coefficients against 0.28 s for the
#' design rebuild. And it went round `autoscale_sdreport()`
#' (`R/autoscale.R:134`), so a fit with `par_units` set got a covariance
#' built on the unscaled Hessian.
#'
#' The cache route has none of those problems: it is core's own object,
#' autoscaling included, and subsetting it is one indexing operation.
#'
#' `fit$cache$Vjoint` is a READ INTO CORE'S INTERNALS and is named as
#' such in the package documentation, not only here. `get_joint_cov()`
#' is `@noRd` and `frmtmb-extension-api` documents neither `fit$cache`
#' nor the `list(V =, names =)` shape of the memo, so unlike `fit$obj`
#' and `fit$estimates` this reach has no precedent to point at. It is
#' made anyway because the alternative was worse, and its failure modes
#' are bounded: see `?frm_curve`, section "The one internal this reaches
#' into", which states them. The standing ask is
#' `dev/spline-seam-proposal.md` Part 1a, an exported accessor.
#'
#' The fallback below runs only when the cache is empty AFTER the check
#' call has run, which happens when `predict(se.fit = TRUE)` returned
#' without core ever forming a joint covariance. It keeps `Q` SPARSE and
#' uses `Matrix::solve()`, for the reason the review gives.
#'
#' `ref` is the result of that check call, and it is an argument rather
#' than a comment ON PURPOSE. The fallback is the OLD route: it goes
#' round `autoscale_sdreport()` (`R/autoscale.R:134`) and would hand an
#' autoscaled fit a covariance built on the unscaled Hessian. It is
#' unreachable today only because the check runs first and warms the
#' cache, and an ordering that load-bearing must not rest on a comment
#' that a later edit can move. Requiring the check's own return value
#' makes the ordering structural: this function cannot be called before
#' the call it depends on, because it needs that call's answer.
#'
#' @noRd
sp_joint_cov <- function(fit, des, ref) {
  if (missing(ref) || is.null(ref) || is.null(ref$se.fit)) {
    stop("frm_curve(): the joint covariance was asked for before the ",
         "predict(se.fit = TRUE) check that fills core's cache. That ",
         "ordering is not cosmetic: without the cache this falls back to ",
         "a fresh sdreport(), which goes round autoscale_sdreport() and ",
         "would return an unscaled covariance on an autoscaled fit. ",
         "Call the check first and pass its result", call. = FALSE)
  }
  jc <- fit$cache$Vjoint
  if (!is.null(jc) && !is.null(jc$V) && !is.null(jc$names)) {
    pos <- sp_coef_pos(fit, des, jc$names, "the cached joint covariance")
    return(as.matrix(jc$V[pos, pos, drop = FALSE]))
  }
  # The cache is empty even though the check has run, so core formed no
  # joint covariance for this fit. The fallback is correct for a fit
  # that carries no autoscaling and wrong for one that does, and there
  # is no exported route to the autoscaled joint precision, so the
  # autoscaled case refuses rather than guesses.
  pu <- fit$par_units
  if (!is.null(pu) && !all(pu == 1)) {
    stop("frm_curve(): this fit is autoscaled (par_units is not all 1) ",
         "and core formed no cached joint covariance for it, so the only ",
         "route left here is a fresh sdreport() on the UNSCALED ",
         "objective, which would not be this fit's covariance. frmtmb ",
         "exports no accessor for the autoscaled joint precision; see ",
         "dev/spline-seam-proposal.md Part 1a", call. = FALSE)
  }
  sdr <- RTMB::sdreport(fit$obj, getJointPrecision = TRUE)
  Q <- sdr$jointPrecision
  if (is.null(Q)) {
    stop("frm_curve(): this fit reports no joint precision matrix, so ",
         "the covariance of a curve cannot be assembled. A model with no ",
         "random-effect block has no penalized smooth in it either",
         call. = FALSE)
  }
  pos <- sp_coef_pos(fit, des, rownames(Q), "the joint precision")
  rest <- setdiff(seq_len(nrow(Q)), pos)
  A <- Q[pos, pos, drop = FALSE]
  S <- if (length(rest)) {
    A - Q[pos, rest, drop = FALSE] %*%
      Matrix::solve(Q[rest, rest, drop = FALSE], Q[rest, pos, drop = FALSE])
  } else {
    A
  }
  as.matrix(Matrix::solve(S))
}

#' Rows of the joint covariance the design's columns sit at.
#'
#' The joint precision numbers `betad` rows over the ESTIMATED entries
#' only, while the design indexes the full `betad` vector, so the two
#' differ by however many entries are held fixed. Every other component
#' is one to one.
#'
#' @noRd
sp_coef_pos <- function(fit, des, rn, what) {
  bd_keep <- setdiff(seq_along(fit$estimates[["betad"]]),
                     fit$frame[["betad_fixed_idx"]] %||% integer(0))
  pos <- integer(length(des$comp))
  for (cp in unique(des$comp)) {
    at <- which(rn == cp)
    hit <- des$comp == cp
    want <- if (cp == "betad") match(des$idx[hit], bd_keep) else des$idx[hit]
    if (anyNA(want) || !length(at) || max(want) > length(at)) {
      stop("frm_curve(): ", what, " carries ", length(at), " '", cp,
           "' rows and the grid design needs ", max(des$idx[hit]),
           ". The fit and its objective disagree about the parameter ",
           "vector; refit before reading a curve off it", call. = FALSE)
    }
    pos[hit] <- at[want]
  }
  pos
}

#' Everything the three exported functions share: the grid, the design,
#' the covariance, and the check that the covariance is the right one.
#'
#' @noRd
sp_curve_parts <- function(fit, newdata, dpar, resp, re.form, tol) {
  if (!inherits(fit, "frmtmb_fit")) {
    stop("`object` must be a frmtmb_fit from frm(), or a frmtmb_curve ",
         "from frm_curve(), not an object of class ", class(fit)[1L],
         call. = FALSE)
  }
  if (!is.data.frame(newdata) || !nrow(newdata)) {
    stop("`newdata` must be a data frame with at least one row: it is ",
         "the grid the curve is evaluated on", call. = FALSE)
  }
  eta0 <- sp_predict_eta(fit, newdata, dpar, resp, re.form)
  if (!sp_check_linear(fit, newdata, dpar, resp, re.form, eta0)) {
    stop("frm_curve(): this linear predictor is not linear in its own ",
         "coefficients, which a nonlinear (nl = TRUE) body makes true. ",
         "The curve covariance is assembled by unit perturbation, and ",
         "that construction holds only for a predictor whose value is a ",
         "linear function of beta and b", call. = FALSE)
  }
  des <- sp_curve_design(fit, newdata, dpar, resp, re.form, eta0)
  # this call comes FIRST because it is what fills fit$cache$Vjoint: the
  # check the package owes anyway is also how the covariance is paid for
  # once instead of twice
  ref <- stats::predict(fit, newdata = newdata, type = "link", dpar = dpar,
                        resp = resp, re.form = re.form, se.fit = TRUE)
  V <- sp_joint_cov(fit, des, ref)
  Sigma <- des$C %*% V %*% t(des$C)
  se <- sqrt(pmax(diag(Sigma), 0))
  se_ref <- as.numeric(ref$se.fit)
  rel <- max(abs(se / pmax(se_ref, .Machine$double.eps) - 1))
  if (!is.finite(rel) || rel > tol) {
    stop("frm_curve(): the assembled curve covariance disagrees with ",
         "predict(se.fit = TRUE) by ", format(rel, digits = 3),
         " relative, which is above the tolerance ", format(tol),
         ". frmtmb exports no route to the joint covariance of a grid ",
         "prediction, so this package rebuilds it and checks it against ",
         "the standard errors that ARE exported. A fit where the two ",
         "disagree is one this package must not report a band for",
         call. = FALSE)
  }
  list(eta = eta0, C = des$C, V = V, Sigma = Sigma, se = se,
       se_ref = se_ref, rel = rel, n_predict = des$n_predict,
       newdata = newdata, dpar = dpar, resp = resp, re.form = re.form,
       fit = fit)
}

#' The max-deviation simultaneous critical value of Ruppert, Wand and
#' Carroll (2003, ch. 6).
#'
#' Draw the curve's own deviation process, standardize it by `div`, take
#' the largest absolute value over the grid, and report the `level`
#' quantile. `div` is a separate argument rather than `sqrt(diag(S))`
#' because the critical value is comparable between two packages only
#' when the divisor is: gratia standardizes a smooth-only deviation by
#' the FULL predictor's standard error, and the critical value that
#' comes back is 8 percent smaller than the self-standardized one on the
#' same fit. The BAND is right either way, because the same divisor
#' calibrates it and scales it.
#'
#' @noRd
sp_sim_crit <- function(S, div, nsim, level, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- nrow(S)
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)
  z <- matrix(stats::rnorm(m * nsim), m, nsim)
  mx <- apply(abs((L %*% z) / div), 2L, max)
  crit <- unname(stats::quantile(mx, level, type = 8))
  # The standard error of a sample quantile, sqrt(p(1-p)/n) / f(q). A
  # critical value reported without it invites a comparison across
  # packages that simulation noise alone would fail.
  d <- stats::density(mx)
  f <- stats::approx(d$x, d$y, xout = crit)$y
  list(crit = crit,
       mcse = if (isTRUE(is.finite(f) && f > 0)) {
         sqrt(level * (1 - level) / nsim) / f
       } else NA_real_)
}
