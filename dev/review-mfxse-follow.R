# Review of lane mfxse, follow-ups of attack 5.
#   Rscript dev/review-mfxse-follow.R
# (a) avg_slopes()/avg_comparisons() at new levels gave no warning
#     (vcovwarn.log): what the slope there omits, and why the check is
#     not reached. Data seed 7411 (mre_data()).
# (b) an exact gp() at x = 9 read 1.0749 times fitted() with no warning:
#     the forward difference, or something else? Data seed 8150.
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(lane, base, user))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
u0 <- rnorm(J, 0, 0.8); u1 <- rnorm(J, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + (0.4 + u1[d$g]) * d$x + 0.3 * (d$f == "b") + u0[d$g]))
fit <- quiet(frm(cnt ~ x + f + (1 | g), family = poisson(), data = d))
ndn <- data.frame(x = c(0, 1, 2), f = "a", g = factor(c("3", "new", "new2")))
# trace the order of the methods marginaleffects calls, and the armed flag
log <- character(0)
gp0 <- ns$get_predict.frmtmb_fit; gv0 <- ns$get_vcov.frmtmb_fit
tr_gp <- function(model, newdata, type = "response", ...) {
  log <<- c(log, sprintf("get_predict rows=%d newlevel=%d armed=%s", nrow(newdata),
                         sum(!as.character(newdata$g) %in% levels(d$g)),
                         isTRUE(model$cache[["mfx_rows_armed"]])))
  gp0(model, newdata, type, ...)
}
tr_gv <- function(model, vcov = NULL, ...) { log <<- c(log, "get_vcov"); gv0(model, vcov, ...) }
registerS3method("get_predict", "frmtmb_fit", tr_gp, envir = asNamespace("marginaleffects"))
registerS3method("get_vcov", "frmtmb_fit", tr_gv, envir = asNamespace("marginaleffects"))
for (call in c("predictions", "avg_slopes", "avg_comparisons")) {
  log <- character(0)
  r <- quiet(switch(call,
    predictions = predictions(fit, newdata = ndn, allow_new_levels = TRUE),
    avg_slopes = avg_slopes(fit, newdata = ndn, variables = "x", allow_new_levels = TRUE),
    avg_comparisons = avg_comparisons(fit, newdata = ndn, variables = "f", allow_new_levels = TRUE)))
  cat("\n", call, ": first calls\n  ", paste(utils::head(log, 6), collapse = "\n   "), "\n")
}
registerS3method("get_predict", "frmtmb_fit", gp0, envir = asNamespace("marginaleffects"))
registerS3method("get_vcov", "frmtmb_fit", gv0, envir = asNamespace("marginaleffects"))
# what the slope at the new levels omits: its variance over the new level's
# effect u ~ N(0, sd_g^2), at the estimates, against the delta SE
sl <- quiet(slopes(fit, newdata = ndn, variables = "x", allow_new_levels = TRUE))
sdg <- exp(fit$estimates$theta[1])
bx <- fit$estimates$beta[2]
eta <- fit$estimates$beta[1] + bx * ndn$x
set.seed(7412)
u <- rnorm(1e5, 0, sdg)
sl_u <- sapply(seq_along(eta), function(i) stats::sd(bx * exp(eta[i] + u)))
cat("\nslopes() at rows 3 / new / new2:", sprintf("%.3f", sl$estimate), "\n")
cat("  delta SE:", sprintf("%.3f", sl$std.error), "\n")
cat("  SD of the slope over the new level's effect (1e5 draws, seed 7412):",
    sprintf("%.3f", sl_u), "(row 1 is a seen level)\n")
# (b) the far gp row
set.seed(8150)
mg <- 90
dg <- data.frame(x = runif(mg, 0, 6), g = factor(sample(1:6, mg, TRUE)))
dg$y <- sin(dg$x) + rnorm(6, 0, 0.5)[dg$g] + rnorm(mg, sd = 0.3)
fg <- quiet(frm(y ~ gp(x) + (1 | g), data = dg))
nd <- data.frame(x = c(7, 8, 9, 10, 12), g = factor("2", levels = levels(dg$g)))
ft <- quiet(fitted(fg, newdata = nd))[, "Est.Error"]
lb <- frm_lp_basis(fg, newdata = nd)
bk <- Filter(function(b) b$covstruct == "gp", fg$frame$re_blocks)[[1]]
v0 <- exp(fg$estimates$theta[bk$theta_idx[1]])^2
pf <- quiet(predictions(fg, newdata = nd))
pc <- quiet(predictions(fg, newdata = nd, numderiv = "fdcenter"))
cat("\ngp + (1|g) at x = 7, 8, 9, 10, 12; theta", sprintf("%.3f", fg$estimates$theta), "\n")
cat("  kriging share   ", sprintf("%.2e", lb$extra_var / v0), "\n")
cat("  sqrt(A V A')/fitted", sprintf("%.4f", sqrt(rowSums((as.matrix(lb$A) %*% lb$V) * as.matrix(lb$A))) / ft), "\n")
cat("  fdforward/fitted", sprintf("%.4f", pf$std.error / ft), "\n")
cat("  fdcenter/fitted ", sprintf("%.4f", pc$std.error / ft), "\n")
