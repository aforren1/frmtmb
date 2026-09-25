# Probe: does simulate(re_formula = NA) redraw a POPULATION smooth's
# penalized coefficients, where predict(re_formula = NA) keeps the smooth?
#   Rscript dev/simnewdata-probe2.R
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
set.seed(21)
n <- 200
d <- data.frame(x = stats::runif(n, 0, 1))
d$y <- sin(2 * pi * d$x) * 2 + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ s(x)), data = d)
cat("blocks:", vapply(fit$frame$re_blocks, `[[`, "", "covstruct"), "\n")
cat("frm_linpred NA vs NULL, max abs diff:",
    max(abs(frm_linpred(fit, re_formula = NA) - frm_linpred(fit))), "\n")
s0 <- as.matrix(simulate(fit, nsim = 400, seed = 1))
s1 <- as.matrix(simulate(fit, nsim = 400, seed = 1, re_formula = NA))
cat(sprintf("sigma %.4f; sd of row draws, mean over rows: NULL %.4f, NA %.4f\n",
            sigma(fit), mean(apply(s0, 1, stats::sd)),
            mean(apply(s1, 1, stats::sd))))
cat(sprintf("RMS of (row mean of draws - fitted): NULL %.4f, NA %.4f\n",
            sqrt(mean((rowMeans(s0) - fitted(fit)[, 1])^2)),
            sqrt(mean((rowMeans(s1) - fitted(fit)[, 1])^2))))
