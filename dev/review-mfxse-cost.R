# Review of lane mfxse, attack 6: where the time goes at 500 levels, and how
# much a restriction to the parameters a prediction reads would give back.
#   Rscript dev/review-mfxse-cost.R
# Data seed 7411 (dev/mfxse-time.R's design: poisson (1 + x | g), 40 rows per
# level). Minimum of 3 rounds per timing.
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lane, base, user))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
quiet <- function(e) suppressWarnings(suppressMessages(e))
el <- function(e) { t0 <- proc.time()[["elapsed"]]; force(e); proc.time()[["elapsed"]] - t0 }
best <- function(f, k = 3) min(replicate(k, el(f())))
J <- 500L; m <- 40L; n <- J * m
set.seed(7411)
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
invisible(quiet(get_vcov(fit)))
cf <- get_coef(fit)
cat("get_coef length", length(cf), "\n")
nd2 <- d[c(1, n), ]          # two rows, levels 1 and 500
cat("levels in nd2:", as.character(nd2$g), "\n")
p2 <- predictions(fit, newdata = nd2)
Jm <- attr(p2, "jacobian") %||% components(p2, "jacobian")
cat("jacobian columns", ncol(Jm), " nonzero columns", sum(colSums(abs(Jm)) > 0), "\n")
cat("nonzero:", paste(colnames(Jm)[colSums(abs(Jm)) > 0], collapse = ", "), "\n")
t_pred2 <- best(function() quiet(predictions(fit, newdata = nd2)))
t_pred2_nov <- best(function() quiet(predictions(fit, newdata = nd2, vcov = FALSE)))
t_gp2 <- best(function() for (i in 1:50) get_predict(fit, newdata = nd2)) / 50
t_set <- best(function() for (i in 1:50) set_coef(fit, cf)) / 50
t_get <- best(function() for (i in 1:50) get_coef(fit)) / 50
t_gpall <- best(function() for (i in 1:5) get_predict(fit, newdata = d)) / 5
cat(sprintf("predictions(nd2): %.2f s; vcov = FALSE %.3f s\n", t_pred2, t_pred2_nov))
cat(sprintf("per call: get_predict(2 rows) %.2f ms, set_coef %.2f ms, get_coef %.2f ms, get_predict(%d rows) %.2f ms\n",
            1e3 * t_gp2, 1e3 * t_set, 1e3 * t_get, n, 1e3 * t_gpall))
nz <- sum(colSums(abs(Jm)) > 0)
cat(sprintf("estimate with the unread b skipped (no get_predict for them): %.3f s\n",
            length(cf) * (t_set + t_get / length(cf)) + (nz + 1) * t_gp2))
# the analytic route: one frm_lp_basis() over all rows (what an exact
# Jacobian of avg_predictions() would need) against avg_predictions()
t_basis <- best(function() frm_lp_basis(fit, newdata = d), 2)
t_avg <- best(function() quiet(avg_predictions(fit)), 2)
cat(sprintf("frm_lp_basis(all %d rows) %.2f s; avg_predictions() %.2f s\n", n, t_basis, t_avg))
# interop_joint_from_outer() inverts Qrr dense: its time and size here
t_cl <- best(function() quiet(predictions(fit, newdata = nd2, vcov = ~ g)), 1)
cat(sprintf("predictions(nd2, vcov = ~ g): %.2f s\n", t_cl))
