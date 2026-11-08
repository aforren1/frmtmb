# Review of lane mfxse, attack 1 (repaired covariance): predictions() SE
# against fitted()'s Est.Error where the standard-error check lost or
# repaired a direction, so get_joint_cov() reads cov_fixed_prop and `null`.
#   Rscript dev/review-mfxse-lost.R [lane|before]
# Data seeds 9401 (gp, no signal), 9402 (separation beside (1 | g)),
# 9403 (sd of (1 | g) at 0), 9404 (correlation at the bound).
a <- commandArgs(TRUE)
arm <- if (length(a)) a[1] else "lane"
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(switch(arm, lane = c(lane, base, user), before = c(base, user)))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "arm:", arm, "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
show <- function(lab, fit, nd) {
  sdr <- ns$sdr_of(fit)
  ws <- character(0)
  p <- tryCatch(withCallingHandlers(predictions(fit, newdata = nd),
    warning = function(w) { ws <<- c(ws, substr(conditionMessage(w), 1, 70)); invokeRestart("muffleWarning") }),
    error = function(e) conditionMessage(e))
  ft <- quiet(fitted(fit, newdata = nd))[, "Est.Error"]
  Vg <- tryCatch(quiet(get_vcov(fit)), error = function(e) NULL)
  ev <- if (!is.null(Vg) && all(is.finite(Vg))) min(eigen((Vg + t(Vg)) / 2, TRUE, only.values = TRUE)$values) else NA
  cat(sprintf("\n%s: conv %d, lost: %s\n", lab, fit$opt$convergence,
              paste(names(sdr$se_lost), sdr$se_lost, sep = "=", collapse = ", ")))
  if (is.character(p)) { cat("  predictions ERROR:", p, "\n"); return() }
  cat("  mfx    ", sprintf("%.4g", p$std.error), "\n")
  cat("  fitted ", sprintf("%.4g", ft), "\n")
  cat("  get_vcov min eigenvalue", sprintf("%.3g", ev), " non-finite entries",
      if (is.null(Vg)) NA else sum(!is.finite(Vg)), "\n")
  if (length(ws)) cat("  warnings:", paste(unique(ws), collapse = " || "), "\n")
}
set.seed(9401)
n <- 80
d1 <- data.frame(x = runif(n, 0, 6))
d1$y <- rbinom(n, 1, 0.5)
f1 <- quiet(frm(y ~ gp(x), family = bernoulli(), data = d1))
show("bernoulli ~ gp(x), no signal", f1, data.frame(x = c(1, 3, 5)))
set.seed(9402)
J <- 10; m <- 20
d2 <- data.frame(g = factor(rep(1:J, each = m)), x = rnorm(J * m))
d2$s <- as.integer(d2$x > 0.8)
d2$y <- rbinom(J * m, 1, plogis(0.3 * d2$x + rnorm(J, 0, 0.7)[d2$g]))
d2$y[d2$s == 1] <- 1L
f2 <- quiet(frm(y ~ x + s + (1 | g), family = bernoulli(), data = d2))
show("bernoulli, s separates, (1 | g)", f2, data.frame(x = c(0, 0.5, 1), s = c(0, 0, 1), g = factor(c(1, 2, 3))))
set.seed(9403)
d3 <- data.frame(g = factor(rep(1:8, each = 10)), x = rnorm(80))
d3$y <- 1 + 0.5 * d3$x + rnorm(80)
f3 <- quiet(frm(y ~ x + (1 | g), data = d3))
show("gaussian, sd(g) at 0", f3, d3[c(1, 21, 41), ])
set.seed(9404)
d4 <- data.frame(g = factor(rep(1:8, each = 15)), x = rnorm(120))
u <- rnorm(8, 0, 0.8)
d4$y <- rpois(120, exp(0.3 + (0.4 + 0.6 * u[d4$g]) * d4$x + u[d4$g]))
f4 <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d4))
show("poisson (1 + x | g), slope = intercept effect (cor 1)", f4, d4[c(1, 31, 61), ])
