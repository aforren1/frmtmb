# Item 8: a marginal tape built from the frame, as a warm start for frm(). Generic in
# the blocks through cs_entry()$vcov(): gaussian identity, one predictor, sigma with an
# intercept only (or se()); gp blocks and positionwise blocks; no autocor or mixture.
# Everything else returns NULL. This is the whole of the branch the review measures.
marg_applicable <- function(frame) {
  fam_ok <- TRUE
  lps <- frame$linpreds
  mu <- Filter(function(l) identical(l$dpar, "mu"), lps)
  length(mu) == 1L && !length(frame$autocor) && !length(frame$me) &&
    is.null(frame$par_template$miss) &&
    all(vapply(frame$re_blocks, function(k) {
      k$covstruct %in% c("gp", "us", "diag", "cs", "homdiag", "smooth", "hsgp")
    }, NA))
}
marg_start <- function(frame, ngrid = TRUE) {
  ns <- asNamespace("frmtmb")
  lps <- frame$linpreds
  mu <- Filter(function(l) identical(l$dpar, "mu"), lps)[[1L]]
  sg <- Filter(function(l) identical(l$dpar, "sigma"), lps)
  X <- as.matrix(mu$X); Z <- mu$Z; y <- frame$y; if (is.list(y)) y <- y[[1L]]; y <- as.numeric(y)
  at <- frame$aterm_values; if (length(at) && is.list(at[[1L]])) at <- at[[1L]]; se2 <- at[["se"]]
  se2 <- if (is.null(se2)) 0 else as.numeric(se2)^2
  blocks <- frame$re_blocks
  Vfn <- function(theta, sig2) {
    # V on numbers or on the tape, all blocks, dense in rows
    "[<-" <- RTMB::ADoverload("[<-")
    V <- 0
    for (bk in blocks) {
      if (bk$covstruct %in% c("smooth") || (bk$covstruct == "us" && bk$dim == 1L)) {
        # one scale on every coefficient: vcov() on a scalar block drops the tape's class
        V <- V + exp(2 * theta[bk$theta_idx[1L]]) * tcrossprod(as.matrix(Z[, bk$c_idx, drop = FALSE]))
        next
      }
      Cb <- ns$cs_entry(bk$covstruct)$vcov(theta[bk$theta_idx], bk)
      Cb <- RTMB::matrix(Cb, bk$dim, bk$dim)
      Zb <- as.matrix(Z[, bk$c_idx, drop = FALSE])
      if (bk$n_levels > 1L) {
        for (l in seq_len(bk$n_levels)) {
          ii <- (l - 1L) * bk$dim + seq_len(bk$dim)
          Zl <- Zb[, ii, drop = FALSE]
          V <- V + Zl %*% Cb %*% t(Zl)
        }
      } else {
        V <- V + Zb %*% Cb %*% t(Zb)
      }
    }
    D <- RTMB::diag(sig2 + se2, length(y))
    V + D
  }
  nll <- function(p) {
    th <- p$theta
    sig2 <- exp(2 * p$betad)
    V <- Vfn(th, sig2)
    e <- y - as.vector(X %*% p$beta)
    -sum(RTMB::dmvnorm(e, 0, V, log = TRUE))
  }
  list(nll = nll, Vfn = Vfn, X = X, y = y)
}
# Coarse scan on numbers over each gp block's rho and the noise, profile of beta and a
# common scale left to the optimizer; returns the start and the seconds.
marg_opt <- function(frame, keep = 2L) {
  ns <- asNamespace("frmtmb")
  ms <- marg_start(frame)
  tp <- frame$par_template
  gpb <- Filter(function(k) identical(k$covstruct, "gp"), frame$re_blocks)
  # the scale of y: residual sd of the fixed design, for the sd and sigma grids
  res <- ms$y - ms$X %*% qr.solve(ms$X, ms$y)
  sy <- sqrt(mean(res^2))
  dmax <- vapply(gpb, function(k) sqrt(max(Reduce(`+`, k$aux_D2))), 0)
  val <- function(th, ls) {
    V <- ms$Vfn(th, exp(2 * ls))
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(Inf)
    wy <- backsolve(R, ms$y, transpose = TRUE); wX <- backsolve(R, ms$X, transpose = TRUE)
    b <- qr.solve(wX, wy)
    e <- wy - wX %*% b
    sum(log(diag(R))) + 0.5 * sum(e^2) + 0.5 * length(ms$y) * log(2 * pi)
  }
  grid <- list()
  th0 <- tp$theta; th0[] <- th0
  for (bi in seq_along(blocks <- frame$re_blocks)) {
    k <- blocks[[bi]]
    if (!identical(k$covstruct, "gp")) th0[k$theta_idx[1]] <- log(sy) - 1
  }
  for (k in gpb) th0[k$theta_idx[1]] <- log(sy)
  res1 <- list()
  for (i in seq_len(max(1L, length(gpb)))) {
    # one sweep per block over rho; the other blocks at d / 4
  }
  for (k in gpb) th0[k$theta_idx[-1]] <- log(sqrt(max(Reduce(`+`, k$aux_D2))) / 4)
  cands <- list()
  for (ib in seq_along(gpb)) {
    k <- gpb[[ib]]
    d <- dmax[ib]
    for (rho in d / 4^(0:4)) for (nr in c(1, 0.1, 1e-2, 1e-3)) {
      th <- th0
      th[k$theta_idx[-1]] <- log(rho)
      cands[[length(cands) + 1L]] <- list(th = th, ls = log(sy * nr))
    }
  }
  vals <- vapply(cands, function(cc) val(cc$th, cc$ls), 0)
  ord <- order(vals)[seq_len(min(keep, length(vals)))]
  nll <- ms$nll
  tpl <- list(beta = tp$beta, betad = tp$betad, theta = tp$theta)
  best <- NULL
  old <- RTMB::TapeConfig(); on.exit(RTMB::TapeConfig(old), add = TRUE)
  RTMB::TapeConfig(matmul = "atomic")
  for (i in ord) {
    t2 <- tpl
    t2$theta <- cands[[i]]$th
    t2$betad[] <- cands[[i]]$ls
    t2$beta[] <- qr.solve(ms$X, ms$y)
    ob <- RTMB::MakeADFun(nll, t2, silent = TRUE)
    o <- nlminb(ob$par, ob$fn, ob$gr, control = list(iter.max = 300, eval.max = 400))
    if (is.null(best) || o$objective < best$o$objective) best <- list(o = o, ob = ob)
  }
  pl <- best$ob$env$parList(best$o$par)
  list(start = pl[c("beta", "betad", "theta")], ll = -best$o$objective, iter = best$o$iterations, evals = best$o$evaluations)
}
