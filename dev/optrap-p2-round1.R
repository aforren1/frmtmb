# Lane optrap, punch round 2: punch round 1's predictor_resid_sd() (the
# sd of the residuals; above the cap of n * p cells no least squares,
# so the raw spread) and newton_polish() (no Cholesky check; only
# nlminb's tolerance names gate it), patched into this process's frmtmb
# namespace, so the round-2 tests can be seen to fail on the round-1
# rules with everything else of the final build. The cap reads the
# same internal option the final build reads. Sourced after
# library(frmtmb); dev/optrap-testfile.R does so with OPTRAP_ROUND1=1.
local({
  ns <- asNamespace("frmtmb")
  resid1 <- function(lp, z, ok) {
    z <- z[ok]
    if (length(z) < 2L) return(NA_real_)
    mag <- max(abs(z))
    X <- if (is.null(lp[["constant"]])) lp[["X"]]
    off <- lp[["offset"]]
    if (is.numeric(off) && length(off) == length(ok) &&
          all(is.finite(off))) {
      z <- z - off[ok]
    }
    if (design_ok(X) && nrow(X) == length(ok) &&
          as.numeric(sum(ok)) * ncol(X) <=
            getOption("frmtmb.unit_ls_cells", 5e6)) {
      Xk <- as.matrix(X[ok, , drop = FALSE])
      z <- as.numeric(base::qr.resid(base::qr(Xk), z))
    }
    s <- stats::sd(z)
    if (!is.finite(s) || !(s > 1e-10 * mag)) NA_real_ else s
  }
  polish1 <- function(obj, opt, bounds, par_units = NULL, control = NULL) {
    p <- opt$par
    np <- length(p)
    tols <- c("rel.tol", "x.tol", "abs.tol", "xf.tol")
    if (!np || np > 100L || length(obj$env$random) ||
          !isTRUE(opt$convergence == 0) || !is.finite(opt$objective) ||
          any(names(control$optCtrl) %in% tols)) {
      return(opt)
    }
    lower <- rep_len(bounds$lower %||% -Inf, np)
    upper <- rep_len(bounds$upper %||% Inf, np)
    u <- rep_len(par_units %||% 1, np)
    v <- opt$objective
    saved <- obj_state_save(obj)
    moved <- FALSE
    steps <- 0L
    for (it in 1:5) {
      g <- tryCatch(drop(obj$gr(p)), error = function(e) NULL)
      H <- tryCatch(obj$he(p), error = function(e) NULL)
      if (is.null(g) || is.null(H) || length(g) != np ||
            !identical(dim(H), c(np, np)) || !all(is.finite(g)) ||
            !all(is.finite(H))) break
      gs <- g * u
      st <- tryCatch(-solve(H * outer(u, u), gs), error = function(e) NULL)
      if (is.null(st) || !all(is.finite(st))) break
      gain <- -sum(gs * st) / 2
      if (!is.finite(gain) || gain <= 1e-6) break
      q <- p + st * u
      if (any(q <= lower) || any(q >= upper)) break
      vq <- tryCatch(obj$fn(q), error = function(e) NA_real_)
      if (!is.finite(vq) || vq >= v) break
      p <- q
      v <- vq
      moved <- TRUE
      steps <- steps + 1L
    }
    if (!moved) {
      obj_state_restore(obj, saved)
      return(opt)
    }
    out <- opt
    out$par <- stats::setNames(p, names(opt$par))
    out$objective <- v
    out$polish <- c(steps = steps, gain = opt$objective - v)
    obj_settle(obj, opt, out, saved)
    out
  }
  environment(resid1) <- ns
  environment(polish1) <- ns
  utils::assignInNamespace("predictor_resid_sd", resid1, "frmtmb")
  utils::assignInNamespace("newton_polish", polish1, "frmtmb")
  cat("predictor_resid_sd, newton_polish: punch round 1's\n")
})
