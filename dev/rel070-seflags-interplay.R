# rel070 seflags, item 4: the degree-7 and degree-8 raw polynomial with
# (1 | g) (the design of dev/seflags-rev4-illpred.R and
# dev/seflags-rev3-deg7exact.R) on one build. Per fit: code, logLik
# against lmer ML (bobyqa), the lost set, optrap's units, every warning
# in full, coefficient SEs and fitted(newdata, re_formula = NA) SEs over
# the exact values (QR of the whitened design at the fit's own sigma and
# sd), and the sigma and sd SEs over the poly(x, deg) fit's.
#   Rscript rel070-seflags-interplay.R <lib> [<lib2>]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "| R", R.home(), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
qr_cov <- function(X, Z, s, t, X0) {
  Vv <- s^2 * diag(nrow(X)) + t^2 * Z %*% t(Z)
  L <- chol(Vv)
  Xw <- backsolve(L, X, transpose = TRUE)
  q <- qr(Xw, tol = 1e-30)
  R <- qr.R(q)
  A <- backsolve(R, t(X0[, q$pivot, drop = FALSE]), transpose = TRUE)
  Ri <- backsolve(R, diag(ncol(X)))
  v <- numeric(ncol(X)); v[q$pivot] <- sqrt(rowSums(Ri^2))
  list(pred = sqrt(colSums(A^2)), coef = v)
}
se_outer <- function(fit) {
  v <- suppressWarnings(vcov(fit, full = TRUE))
  stats::setNames(sqrt(diag(as.matrix(v))), ns$outer_par_names(fit))
}
f4 <- function(x) paste(formatC(x, digits = 4, format = "f"), collapse = " ")
for (deg in 7:8) for (s in c(3, 4, 7, 8, 10)) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
    stats::rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  rhs <- paste0("x", 1:deg, collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  g <- suppressWarnings(frm(stats::as.formula(
    paste0("y ~ poly(x, ", deg, ") + (1 | g)")), data = d))
  m <- suppressWarnings(suppressMessages(lme4::lmer(
    stats::as.formula(paste("y ~ poly(x,", deg, ") + (1 | g)")), data = d,
    REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 1e5)))))
  on <- ns$outer_par_names(f)
  l <- suppressWarnings(ns$sdr_of(f)$se_lost)
  u <- if (exists("predictor_units", ns)) ns$predictor_units(f$frame) else NULL
  cat(sprintf(paste0("\n== degree %d seed %d | code %d | logLik %.4f, lmer ",
                     "ML %.4f, poly() fit %.4f (short %.4f) | units %s | ",
                     "lost %s\n"),
              deg, s, f$opt$convergence, as.numeric(logLik(f)),
              as.numeric(logLik(m)), as.numeric(logLik(g)),
              as.numeric(logLik(m)) - as.numeric(logLik(f)),
              if (length(u)) paste(names(u), u) else "-",
              if (length(l)) paste(names(l), l, sep = ":", collapse = " ")
              else "-"))
  for (w in r$w) cat("  W:", gsub("\n", " ", w), "\n")
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
  nd <- data.frame(x = c(1.0, 1.5, 2.0, 2.5), g = d$g[1])
  for (k in 1:deg) nd[[paste0("x", k)]] <- nd$x^k
  X0 <- stats::model.matrix(stats::as.formula(paste("~", rhs)), nd)
  Z <- stats::model.matrix(~ 0 + g, d)
  sg <- exp(f$opt$par[grepl("^sigma", on)])
  tg <- exp(f$opt$par[on == "theta_1"])
  ex <- qr_cov(X, Z, sg, tg, X0)
  se <- se_outer(f)
  cn <- c("(Intercept)", paste0("x", 1:deg))
  p <- cap(fitted(f, newdata = nd, re_formula = NA))
  cat("  coef SE / exact (QR):   ", f4(se[cn] / ex$coef), "\n")
  cat("  fitted SE / exact (QR): ", f4(p$v[, "Est.Error"] / ex$pred),
      "| call warnings", length(p$w), "\n")
  sg2 <- se_outer(g)
  og <- ns$outer_par_names(g)
  cat("  sigma, sd SE / poly() fit's:",
      f4(c(se[grepl("^sigma", on)] / sg2[grepl("^sigma", og)],
           se[on == "theta_1"] / sg2[og == "theta_1"])),
      "| sigma", signif(sg, 5), "sd", signif(tg, 5), "\n")
}
