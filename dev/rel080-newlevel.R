# Brief step 8 (mfxse x readers): avg_predictions() at a grouping level the
# fit did not see. The mean there is exp(eta + sd^2 / 2) on a poisson log
# link (readers), so its standard error reads the block's sd as well as
# the coefficients. Reference: the delta method by hand over the outer
# covariance, and a Monte Carlo of the outer parameters (seed 7).
#
#   Rscript dev/rel080-newlevel.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1L], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
set.seed(7411)
J <- 20; n <- 400
d <- data.frame(x = rnorm(n), g = factor(sample(J, n, TRUE)))
d$y <- rpois(n, exp(0.3 + 0.5 * d$x + rnorm(J, 0, 1)[d$g]))
fit <- frm(y ~ x + (1 | g), data = d, family = poisson())
nd <- data.frame(x = c(-1, 0, 1), g = "new")
ws <- character()
ap <- withCallingHandlers(
  avg_predictions(fit, newdata = nd, allow_new_levels = TRUE),
  warning = function(w) {
    ws <<- c(ws, class(w)[1L])
    invokeRestart("muffleWarning")
  })
fv <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
V <- vcov(fit, full = TRUE)
cat("outer parameters:", rownames(V), "\n")
est <- c(fixef(fit)[, "Estimate"])
sd_g <- VarCorr(fit)$g$sd[1, 1]
cat(sprintf("fitted() mean at the new level, averaged %.6f; avg_predictions %.6f\n",
            mean(fv[, "Estimate"]), ap$estimate))
# the outer vector: two coefficients, then the block's log sd
th <- grep("theta|sd", rownames(V))
stopifnot(length(th) == 1L, nrow(V) == 3L)
mu_of <- function(b0, b1, lsd) {
  rowMeans(exp(outer(b0, rep(1, 3)) + outer(b1, nd$x) + exp(2 * lsd) / 2))
}
p0 <- c(est[1:2], log(sd_g))
cat(sprintf("closed form at the estimates %.6f\n", mu_of(p0[1], p0[2], p0[3])))
gr <- vapply(1:3, function(k) {
  h <- 1e-5 * max(1, abs(p0[k])); a <- p0; b <- p0
  a[k] <- a[k] + h; b[k] <- b[k] - h
  (mu_of(a[1], a[2], a[3]) - mu_of(b[1], b[2], b[3])) / (2 * h)
}, 0)
se_full <- sqrt(drop(t(gr) %*% V %*% gr))
se_coef <- sqrt(drop(t(gr[1:2]) %*% V[1:2, 1:2] %*% gr[1:2]))
set.seed(7)
P <- MASS::mvrnorm(2e5, p0, V)
mc <- mu_of(P[, 1], P[, 2], P[, 3])
cat(sprintf("SE: avg_predictions %.6f | delta, coefficients and sd %.6f | delta, coefficients alone %.6f | MC sd %.6f\n",
            ap$std.error, se_full, se_coef, sd(mc)))
cat(sprintf("ratio avg_predictions / delta with sd %.5f\n", ap$std.error / se_full))
cat("warnings:", if (length(ws)) paste(ws, collapse = ", ") else "none", "\n")
t0 <- proc.time()[["elapsed"]]
invisible(avg_predictions(fit, newdata = nd, allow_new_levels = TRUE))
cat(sprintf("time of one avg_predictions() at the new level %.2f s\n",
            proc.time()[["elapsed"]] - t0))
cat("NEWLEVEL DONE\n")
