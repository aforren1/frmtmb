# Reviewer of lane seflags, final check of B4: is the re-measured
# curvature of an "illcond" direction accurate? If so the direction can
# go back into the covariance, V + d diag(1 / lambda) d', instead of
# taking SEs away. Coefficient SEs from that covariance against the
# exact ones (QR of the whitened design at the fit's sigma and sd), for
# the illcond fits of dev/seflags-rev2-illpred.R.
#   Rscript dev/seflags-rev3-illfull.R <lane|op3>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
for (cs in list(c(8, 3, 1), c(8, 8, 1), c(8, 10, 1), c(8, 4, 0),
                c(9, 3, 1), c(9, 4, 0))) {
  deg <- cs[1]; s <- cs[2]; re <- cs[3] == 1
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  f <- suppressWarnings(frm(stats::as.formula(paste("y ~", rhs,
                                                     if (re) "+ (1 | g)")),
                            data = d))
  sdr <- suppressWarnings(ns$sdr_of(f))
  il <- sdr$se_null_flat$ill
  on <- ns$outer_par_names(f)
  if (is.null(il)) {
    cat(sprintf("degree %d seed %d re %s | code %d | lost %s | no illcond record\n",
                deg, s, re, f$opt$convergence, substr(lost_of(f), 1, 60)))
    next
  }
  Vf <- il$V + il$d %*% (t(il$d) / il$lam)
  cn <- c("(Intercept)", paste0("x", 1:deg))
  cn <- intersect(cn, on[il$rows])
  j <- match(cn, on[il$rows])
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)[, cn, drop = FALSE]
  sg <- exp(f$opt$par[grepl("^sigma", on)])
  tg <- if (re) exp(f$opt$par[on == "theta_1"]) else 0
  Vv <- sg^2 * diag(nrow(X)) + if (re) tg^2 * tcrossprod(stats::model.matrix(~ 0 + g, d)) else 0
  q <- qr(backsolve(chol(Vv), X, transpose = TRUE), tol = 1e-30)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  ex <- numeric(ncol(X)); ex[q$pivot] <- sqrt(rowSums(Ri^2))
  cat(sprintf(paste0("degree %d seed %d re %s | code %d | %d illcond dirs,",
                     " lambda %s\n  coef SE with the directions back / exact: %s\n"),
              deg, s, re, f$opt$convergence, length(il$lam),
              paste(signif(il$lam, 3), collapse = " "),
              paste(format(sqrt(diag(Vf))[j] / ex, digits = 4), collapse = " ")))
}
