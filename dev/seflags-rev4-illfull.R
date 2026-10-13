# Reviewer of lane seflags, check of punch 2b: on the "illcond" fits of
# dev/seflags-rev2-illpred.R, the coefficient SEs as fixef() and vcov()
# report them, against the exact ones (QR of the whitened design at the
# fit's sigma and sd), with the fit's warnings.
#   Rscript dev/seflags-rev4-illfull.R <lane|op4>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-common.R")
for (cs in list(c(8, 3, 1), c(8, 8, 1), c(8, 10, 1), c(8, 4, 0),
                c(8, 4, 1), c(9, 3, 1), c(9, 4, 0))) {
  deg <- cs[1]; s <- cs[2]; re <- cs[3] == 1
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  r <- cap(frm(stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)")),
               data = d))
  f <- r$v
  on <- ns$outer_par_names(f)
  cn <- intersect(c("(Intercept)", paste0("x", 1:deg)), on)
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)[, cn,
                                                                  drop = FALSE]
  sg <- exp(f$opt$par[grepl("^sigma", on)])
  tg <- if (re) exp(f$opt$par[on == "theta_1"]) else 0
  Vv <- sg^2 * diag(nrow(X)) +
    if (re) tg^2 * tcrossprod(stats::model.matrix(~ 0 + g, d)) else 0
  q <- qr(backsolve(chol(Vv), X, transpose = TRUE), tol = 1e-30)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  ex <- numeric(ncol(X)); ex[q$pivot] <- sqrt(rowSums(Ri^2))
  fe <- cap(fixef(f))
  vc <- sqrt(diag(stats::vcov(f)))
  cat(sprintf(paste0("degree %d seed %d re %s | code %d | lost %s | fit warnings %d%s\n",
                     "  fixef Est.Error / exact: %s (warn %d)\n",
                     "  sqrt(diag(vcov())) / exact: %s\n"),
              deg, s, re, f$opt$convergence, lost_of(f), length(r$w),
              if (length(r$w)) paste0(" (", substr(r$w[1], 1, 60), ")") else "",
              paste(format(fe$v[, "Est.Error"] / ex, digits = 4), collapse = " "),
              length(fe$w),
              paste(format(vc[seq_along(ex)] / ex, digits = 4), collapse = " ")))
}
