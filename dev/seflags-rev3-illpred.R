# Reviewer of lane seflags, re-check: on an "illcond" fit the
# coefficients have no SE but fitted(newdata) keeps one. Is it right?
# The prediction variance x0' (X' V^-1 X)^-1 x0 at the fit's own sigma
# and sd, computed accurately by QR of the whitened design (accurate to
# the square root of the condition number), against fitted()'s SE, for
# degree-8 and degree-9 raw polynomials, rows in and beyond the data.
# Also which coefficients keep an SE.
#   Rscript dev/seflags-rev3-illpred.R <lane|base|op2>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
qr_cov_rows <- function(X, Z, s, t, X0) {
  Vv <- s^2 * diag(nrow(X)) + if (is.null(Z)) 0 else t^2 * Z %*% t(Z)
  L <- chol(Vv)
  Xw <- backsolve(L, X, transpose = TRUE)
  q <- qr(Xw, tol = 1e-30)
  R <- qr.R(q)
  A <- backsolve(R, t(X0[, q$pivot, drop = FALSE]), transpose = TRUE)
  list(pred = sqrt(colSums(A^2)),
       coef = {
         Ri <- backsolve(R, diag(ncol(X)))
         v <- numeric(ncol(X)); v[q$pivot] <- sqrt(rowSums(Ri^2)); v
       })
}
for (cs in list(c(8, 3, 1), c(8, 8, 1), c(8, 10, 1), c(8, 4, 0),
                c(8, 4, 1), c(9, 3, 1), c(9, 4, 0))) {
  deg <- cs[1]; s <- cs[2]; re <- cs[3] == 1
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  on <- ns$outer_par_names(f)
  cn <- c("(Intercept)", paste0("x", 1:deg))
  nd <- data.frame(x = c(1.0, 1.5, 2.0, 2.5))
  for (k in 1:deg) nd[[paste0("x", k)]] <- nd$x^k
  if (re) nd$g <- d$g[1]
  p <- cap(fitted(f, newdata = nd, re_formula = if (re) NA else NULL))
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  X0 <- stats::model.matrix(stats::as.formula(paste("~", rhs)), nd)
  sg <- exp(f$opt$par[grepl("^sigma", on)])
  tg <- if (re) exp(f$opt$par[on == "theta_1"]) else 0
  Z <- if (re) stats::model.matrix(~ 0 + g, d) else NULL
  ex <- qr_cov_rows(X, Z, sg, tg, X0)
  se <- se_outer(f)[cn]
  cat(sprintf(paste0("degree %d seed %d re %s | code %d | lost %s\n",
                     "  fitted se %s (warn %d)\n  accurate  %s\n",
                     "  coef se/accurate %s\n"),
              deg, s, re, f$opt$convergence, substr(lost_of(f), 1, 90),
              paste(signif(p$v[, 2], 4), collapse = " "), length(p$w),
              paste(signif(ex$pred, 4), collapse = " "),
              paste(format(se / ex$coef, digits = 3), collapse = " ")))
}
