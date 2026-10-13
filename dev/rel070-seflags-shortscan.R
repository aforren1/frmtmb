# rel070 seflags, item 4 widened: raw polynomials of degree 5 to 8 with
# (1 | g), the design of rel070-seflags-interplay.R, seeds 1 to 20. Per
# fit: code, the shortfall against lmer ML on poly(x, deg), the number of
# warnings and whether one names stopping short, the lost count, and the
# largest |coef SE / exact - 1| (exact: QR of the whitened design at the
# fit's own sigma and sd), as TSV on stdout.
#   Rscript rel070-seflags-shortscan.R <lib> [<lib2>]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
qr_coef <- function(X, Z, s, t) {
  L <- chol(s^2 * diag(nrow(X)) + t^2 * Z %*% t(Z))
  q <- qr(backsolve(L, X, transpose = TRUE), tol = 1e-30)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  v <- numeric(ncol(X)); v[q$pivot] <- sqrt(rowSums(Ri^2)); v
}
cat("lib\tdeg\tseed\tcode\tshort\tn_warn\tshort_named\tlarge_grad\tn_lost\tunit\tmax_se_dev\n")
for (deg in 5:8) for (s in 1:20) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  r <- cap(frm(stats::as.formula(paste("y ~", rhs, "+ (1 | g)")), data = d))
  f <- r$v
  m <- suppressWarnings(suppressMessages(lme4::lmer(
    stats::as.formula(paste("y ~ poly(x,", deg, ") + (1 | g)")), data = d,
    REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 1e5)))))
  on <- ns$outer_par_names(f)
  l <- suppressWarnings(ns$sdr_of(f)$se_lost)
  u <- if (exists("predictor_units", ns)) ns$predictor_units(f$frame) else NULL
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  Z <- stats::model.matrix(~ 0 + g, d)
  ex <- qr_coef(X, Z, exp(f$opt$par[grepl("^sigma", on)]),
                exp(f$opt$par[on == "theta_1"]))
  v <- suppressWarnings(vcov(f, full = TRUE))
  se <- stats::setNames(suppressWarnings(sqrt(diag(as.matrix(v)))), on)
  cn <- c("(Intercept)", paste0("x", 1:deg))
  dev <- suppressWarnings(max(abs(se[cn] / ex - 1)))
  cat(basename(a[1]), deg, s, f$opt$convergence,
      sprintf("%.4f", as.numeric(logLik(m)) - as.numeric(logLik(f))),
      length(r$w), any(grepl("stopped short", r$w)),
      any(grepl("^Large maximum absolute gradient", r$w)), length(l),
      if (length(u)) u[[1]] else 1, signif(dev, 4), sep = "\t")
  cat("\n")
}
