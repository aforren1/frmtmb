# Shared by the review scripts of lane ncgp. No frmtmb code in the marginal:
# it is written from the model statement, y ~ N(mu, sd^2 Z (K + nug I) Z' +
# sigma^2 I), in position space so that replicates stay cheap.
REV_LANE <- "C:/Users/adf44/source/r/wt-ncgp-lib"
REV_BASE <- "C:/Users/adf44/source/r/rellib-r16"
REV_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
REV_WT <- "C:/Users/adf44/source/r/frmtmb-wt-ncgp"
REV_OUT <- file.path(REV_WT, "dev/review-ncgp-out")
dir.create(REV_OUT, showWarnings = FALSE)
rev_libs <- function(arm = c("lane", "base")) {
  arm <- match.arg(arm)
  .libPaths(c(if (arm == "lane") REV_LANE, REV_BASE, REV_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(REV_WT, "dev/stan-cache"),
             NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "\n")
  invisible(arm)
}

# Sufficient statistics of a gaussian y over one exact field.
cf_setup <- function(y, x, X = NULL, nug = 1e-12) {
  x <- as.matrix(x)
  key <- apply(x, 1, function(r) paste(format(r, digits = 17), collapse = "_"))
  uk <- unique(key)
  J <- match(key, uk)
  m <- length(uk)
  ux <- x[match(uk, key), , drop = FALSE]
  N <- tabulate(J, m)
  D2 <- 0
  for (j in seq_len(ncol(ux))) D2 <- D2 + outer(ux[, j], ux[, j], "-")^2
  if (is.null(X)) X <- matrix(1, length(y), 1)
  list(y = y, J = J, m = m, N = N, D2 = D2, n = length(y), nug = nug, X = X)
}

# Profile: mu and sd^2 out, given rho and r = sigma^2 / sd^2.
cf_prof <- function(s, lrho, lr) {
  r <- exp(lr)
  C <- exp(-s$D2 / (2 * exp(2 * lrho)))
  diag(C) <- diag(C) + s$nug
  sN <- sqrt(s$N)
  S <- C * outer(sN, sN)
  diag(S) <- diag(S) + r
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) return(list(ll = -Inf))
  ldS <- 2 * sum(log(diag(R)))
  ld <- (s$n - s$m) * lr + ldS
  gsum <- function(v) rowsum_fast(v, s$J, s$m)
  qfm <- function(a, b) {
    ua <- gsum(a) / sN
    ub <- gsum(b) / sN
    # within-position SS from centered rows: a'b - ua'ub cancels at small r
    wa <- a - (ua / sN)[s$J]
    wb <- b - (ub / sN)[s$J]
    sum(wa * wb) / r +
      sum(backsolve(R, ua, transpose = TRUE) * backsolve(R, ub, transpose = TRUE))
  }
  p <- ncol(s$X)
  XtX <- matrix(0, p, p)
  Xty <- numeric(p)
  for (i in seq_len(p)) {
    Xty[i] <- qfm(s$X[, i], s$y)
    for (j in i:p) XtX[i, j] <- XtX[j, i] <- qfm(s$X[, i], s$X[, j])
  }
  beta <- solve(XtX, Xty)
  e <- as.vector(s$y - s$X %*% beta)
  q <- qfm(e, e)
  sd2 <- q / s$n
  list(ll = -0.5 * s$n * log(2 * pi * sd2) - 0.5 * ld - 0.5 * s$n, sd = sqrt(sd2),
       sigma = sqrt(r * sd2), beta = beta, rho = exp(lrho), r = r)
}
rowsum_fast <- function(v, J, m) {
  out <- numeric(m)
  for (i in seq_along(v)) out[J[i]] <- out[J[i]] + v[i]
  out
}

# ll at a full point (mu, sigma, sd, rho), without the profile.
cf_ll_at <- function(s, beta, sigma, sd, rho) {
  r <- (sigma / sd)^2
  C <- exp(-s$D2 / (2 * rho^2))
  diag(C) <- diag(C) + s$nug
  sN <- sqrt(s$N)
  S <- C * outer(sN, sN)
  diag(S) <- diag(S) + r
  R <- chol(S)
  e <- as.vector(s$y - s$X %*% beta)
  ue <- rowsum_fast(e, s$J, s$m) / sN
  w <- e - (ue / sN)[s$J]
  q <- sum(w^2) / r + sum(backsolve(R, ue, transpose = TRUE)^2)
  ld <- (s$n - s$m) * log(r) + 2 * sum(log(diag(R)))
  -0.5 * s$n * log(2 * pi * sd^2) - 0.5 * ld - 0.5 * q / sd^2
}

# Global maximum: a grid over the two profile coordinates, then a polish from
# the best few cells and from any extra start.
cf_max <- function(s, extra = list(), ngrid = 36L, lr_lo = -46, lr_hi = 8) {
  d <- sqrt(s$D2[s$D2 > 0])
  lo <- log(min(d)) - 1
  hi <- log(max(d)) + 6
  g <- expand.grid(lrho = seq(lo, hi, length.out = ngrid),
                   lr = seq(lr_lo, lr_hi, length.out = ngrid))
  v <- mapply(function(a, b) cf_prof(s, a, b)$ll, g$lrho, g$lr)
  ord <- order(-v)
  # distinct cells only: the best cell of each local region
  starts <- list()
  for (i in ord) {
    if (!is.finite(v[i])) next
    p <- c(g$lrho[i], g$lr[i])
    if (all(vapply(starts, function(q) sum(abs(q - p)) > 4, NA))) {
      starts[[length(starts) + 1L]] <- p
    }
    if (length(starts) == 6L) break
  }
  starts <- c(starts, extra)
  best <- NULL
  f <- function(p) {
    z <- cf_prof(s, p[1], p[2])$ll
    if (is.finite(z)) -z else 1e300
  }
  for (p in starts) {
    o <- nlminb(p, f, control = list(iter.max = 300, eval.max = 400,
                                     rel.tol = 1e-13, x.tol = 1e-10))
    if (!all(is.finite(o$par))) next
    o2 <- optim(o$par, f, method = "Nelder-Mead",
                control = list(maxit = 800, reltol = 1e-14))
    o3 <- nlminb(o2$par, f, control = list(iter.max = 200, rel.tol = 1e-14))
    if (is.null(best) || o3$objective < best$objective) best <- o3
  }
  cf_prof(s, best$par[1], best$par[2])
}

# 256-bit reading of cf_ll_at, for the rows double cannot settle.
cf_ll_mpfr <- function(s, beta, sigma, sd, rho, prec = 256L) {
  stopifnot(requireNamespace("Rmpfr", quietly = TRUE))
  mp <- function(x) Rmpfr::mpfr(x, prec)
  r <- mp(sigma)^2 / mp(sd)^2
  m <- s$m
  sN <- sqrt(mp(s$N))
  a <- matrix(list(), m, m)
  # a column-major Cholesky on an mpfr array
  S <- Rmpfr::mpfrArray(0, prec, dim = c(m, m))
  rho2 <- 2 * mp(rho)^2
  for (j in seq_len(m)) for (i in j:m) {
    v <- exp(-mp(s$D2[i, j]) / rho2) * sN[i] * sN[j]
    if (i == j) v <- v + mp(s$nug) * mp(s$N[i]) + r
    S[i, j] <- v
  }
  L <- S
  for (j in seq_len(m)) {
    if (j > 1L) {
      for (k in seq_len(j - 1L)) {
        L[j:m, j] <- L[j:m, j] - L[j:m, k] * L[j, k]
      }
    }
    L[j, j] <- sqrt(L[j, j])
    if (j < m) L[(j + 1L):m, j] <- L[(j + 1L):m, j] / L[j, j]
  }
  e <- mp(s$y) - mp(as.vector(s$X %*% beta))
  ue <- mp(rep(0, m))
  for (i in seq_along(s$y)) ue[s$J[i]] <- ue[s$J[i]] + e[i]
  ue <- ue / sN
  w <- ue
  for (i in seq_len(m)) {
    if (i > 1L) w[i] <- w[i] - sum(L[i, seq_len(i - 1L)] * w[seq_len(i - 1L)])
    w[i] <- w[i] / L[i, i]
  }
  wr <- e - (ue / sN)[s$J]
  q <- sum(wr^2) / r + sum(w^2)
  ld <- (s$n - m) * log(r)
  for (i in seq_len(m)) ld <- ld + 2 * log(L[i, i])
  ll <- -mp(s$n) / 2 * log(2 * Rmpfr::Const("pi", prec) * mp(sd)^2) - ld / 2 -
    q / (2 * mp(sd)^2)
  as.numeric(ll)
}

# What frm() read on a y ~ gp(x) fit, warnings and messages kept.
rev_fit <- function(fml, data, ...) {
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frmtmb::frm(fml, data = data, ...),
             error = function(e) e),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warn = w)
}
rev_gp_read <- function(fit) {
  bks <- Filter(function(k) identical(k$covstruct, "gp"), fit$frame$re_blocks)
  bk <- bks[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  list(ll = as.numeric(logLik(fit)), sd = exp(th[1L]), rho = exp(th[-1L]),
       sigma = tryCatch(exp(fit$estimates$betad[[1L]]), error = function(e) NA),
       beta = fit$estimates$beta, conv = fit$opt$convergence)
}
