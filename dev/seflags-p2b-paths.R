# Lane seflags, punch round 2b: every route a prediction's SE takes on
# an "illcond" fit, against the exact SE (QR of the whitened design at
# the fit's sigma and sd). Degree-8 raw polynomial, seed 8 with (1 | g)
# and seed 4 without; rows x = 1, 1.5, 2, 2.5 at re_formula = NA.
#   Rscript dev/seflags-p2b-paths.R <lane|m9|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
ratio <- function(a, b) paste(format(a / b, digits = 4), collapse = " ")
for (cs in list(c(8, 1), c(4, 0))) {
  s <- cs[1]; re <- cs[2] == 1
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:8, collapse = " + ")
  r <- cap(frm(stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)")),
               data = d))
  f <- r$v
  on <- ns$outer_par_names(f)
  nd <- data.frame(x = c(1, 1.5, 2, 2.5), g = d$g[1])
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  X0 <- stats::model.matrix(stats::as.formula(paste("~", rhs)), nd)
  sg <- exp(f$opt$par[grepl("^sigma", on)])
  Vv <- sg^2 * diag(nrow(X))
  if (re) {
    tg <- exp(f$opt$par[on == "theta_1"])
    Vv <- Vv + tg^2 * tcrossprod(stats::model.matrix(~ 0 + g, d))
  }
  q <- qr(backsolve(chol(Vv), X, transpose = TRUE), tol = 1e-30)
  A <- backsolve(qr.R(q), t(X0[, q$pivot, drop = FALSE]), transpose = TRUE)
  ex <- sqrt(colSums(A^2))
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  exc <- numeric(ncol(X)); exc[q$pivot] <- sqrt(rowSums(Ri^2))
  rf <- if (re) NA else NULL
  cat(sprintf("\n== seed %d re %s | code %d | lost %s | fit warnings %d: %s\n",
              s, re, f$opt$convergence, lost_of(f), length(r$w),
              substr(paste(r$w, collapse = " | "), 1, 120)))
  cat("  exact rows   ", format(ex, digits = 4), "\n")
  cat("  coef SE       ", ratio(se_outer(f)[colnames(X)], exc), "\n")
  a <- cap(fitted(f, newdata = nd, re_formula = rf))
  cat("  fitted()      ", ratio(a$v[, "Est.Error"], ex), "warn", length(a$w), "\n")
  a <- cap(frm_linpred(f, newdata = nd, re_formula = rf, se.fit = TRUE))
  cat("  frm_linpred() ", ratio(a$v$se.fit, ex), "warn", length(a$w), "\n")
  a <- cap(frm_lp_basis(f, newdata = nd, re_formula = rf, extra_cov = TRUE))
  lb <- a$v
  v1 <- rowSums((lb$A %*% lb$V) * lb$A) + lb$extra_var
  v2 <- diag(as.matrix(lb$A %*% lb$V %*% t(lb$A) + lb$extra_cov))
  cat("  lp_basis var  ", ratio(sqrt(v1), ex), "| full cov",
      ratio(sqrt(v2), ex), "\n")
  jc <- cap(frm_joint_cov(f))$v
  cat("  joint cov diag", ratio(sqrt(diag(jc$V))[match(colnames(X), on)], exc),
      "\n")
  h <- paste0("Intercept + ", paste0(1.5^(1:8), " * x", 1:8, collapse = " + "),
              " = 0")
  a <- cap(hypothesis(f, h))
  cat("  hypothesis x=1.5", ratio(a$v$hypothesis$Est.Error, ex[2]), "warn",
      length(a$w), "\n")
  lin <- function(g) {
    matrix(suppressWarnings(frm_linpred(g, newdata = nd, re_formula = rf)),
           ncol = 1L)
  }
  a <- cap(ns$fit_fd_se(f, lin))
  cat("  fit_fd_se()   ", ratio(as.vector(a$v), ex), "warn", length(a$w), "\n")
  set.seed(1)
  dr <- ns$predict_par_drawer(f, TRUE, 2000L)
  E <- vapply(seq_len(2000L), function(s) as.vector(lin(dr(s))), numeric(4))
  cat("  draws sd      ", ratio(apply(E, 1, stats::sd), ex),
      "(2000 draws, MC 1.6 percent)\n")
}
