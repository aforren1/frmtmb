# Short check: why predread's 0.46 moved. The same design, seed 7403: ll, sd, lscale and SEs on each build; my closed form with z in X.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs(arm)
set.seed(7403)
m <- 80
g <- data.frame(x = sort(runif(m, 0, 6)), z = rnorm(80))
g$y <- sin(g$x) + 0.3 * g$z + rnorm(m, 0, 0.25)
f <- suppressWarnings(frm(bf(y ~ z + gp(x)), family = gaussian(), data = g))
s <- cf_setup(g$y, g$x, X = cbind(1, g$z), nug = if (arm == "lane") 1e-12 else 1e-6)
th <- f$estimates$theta
cf <- cf_max(s, extra = list(c(th[2], 2 * (f$estimates$betad[[1]] - th[1]))), ngrid = 30L)
cat(sprintf("RMY %s ll %.5f sd %.4f rho %.4f sigma %.4f | cf ll %.5f sd %.4f rho %.4f sigma %.4f | SE beta %s\n", arm, logLik(f), exp(th[1]), exp(th[2]), exp(f$estimates$betad[[1]]), cf$ll, cf$sd, cf$rho, cf$sigma,
            paste(signif(sqrt(diag(as.matrix(vcov(f)))), 4), collapse = " ")))
