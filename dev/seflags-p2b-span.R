# Lane seflags, punch round 2b: where does the remaining few percent of
# error in the coefficient SEs of an "illcond" fit come from? Degree 8,
# no random effects, seeds 6 and 7 (SE / exact 1.02 and 1.04 with the
# direction put back). The exact AD Hessian, its unit-diagonal
# eigendecomposition, and the covariance with the span of the k
# smallest directions re-measured by se_span_remeasure(), k = 1 to 4,
# against the exact SEs (QR of the design at the fit's sigma).
#   Rscript dev/seflags-p2b-span.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
for (s in c(6, 7, 3)) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:8, collapse = " + ")
  f <- suppressWarnings(frm(stats::as.formula(paste("y ~", rhs)), data = d))
  p <- f$opt$par
  on <- ns$outer_par_names(f)
  H <- f$obj$he(p)
  H <- (H + t(H)) / 2
  free <- seq_along(p)
  Df <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(Df, Df), symmetric = TRUE)
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  sg <- exp(p[grepl("^sigma", on)])
  q <- qr(X / sg, tol = 1e-30)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  ex <- numeric(ncol(X))
  ex[q$pivot] <- sqrt(rowSums(Ri^2))
  j <- match(colnames(X), on)
  n <- length(e$values)
  cat(sprintf("seed %d | smallest eigenvalues %s\n", s,
              paste(signif(rev(e$values)[1:5], 3), collapse = " ")))
  for (k in 1:4) {
    S <- (n - k + 1):n
    sp <- ns$se_span_remeasure(f, p, free, e, Df, S)
    if (is.null(sp)) {
      cat("  k", k, "not confirmed\n")
      next
    }
    ev <- e$values
    U <- e$vectors
    U[, S] <- sp$vectors
    ev[S] <- sp$values
    V <- U %*% (t(U) / ev) / outer(Df, Df)
    cat(sprintf("  span of %d re-measured: values %s confirmed %s | SE/exact %s\n",
                k, paste(signif(sp$values, 4), collapse = " "),
                paste(sp$confirmed, collapse = " "),
                paste(range(signif(sqrt(diag(V))[j] / ex, 4)),
                      collapse = " to ")))
  }
}
