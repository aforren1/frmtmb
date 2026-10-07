# Lane setier: the gp(x, k = 10) fit of test-brms-likelihood.R row 10:
# its two hyperparameters, their Hessian rows and the edge probes.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(as.integer(commandArgs(TRUE)[1] %||% 1))
n <- 60
d <- data.frame(x = sort(runif(n, 0, 10)))
d$y <- sin(d$x) + rnorm(n, 0, 0.3)
fit <- suppressWarnings(frm(bf(y ~ gp(x, k = 10)) + gaussian(), data = d))
print(fit$opt$par)
h <- fit$cache$hessian_fixed
print(signif(h$H, 4)); print(signif(h$E, 3))
print(ns$sdr_of(fit)$se_lost)
S <- h$H / sqrt(outer(diag(h$H), diag(h$H)))
print(eigen(S)$values)
p <- fit$opt$par
f0 <- fit$obj$fn(p)
for (j in 1:2) for (s in c(-20, -2, 2, 20)) {
  q <- p; q[j] <- q[j] + s
  cat("theta", j, "step", s, "dll", fit$obj$fn(q) - f0, "\n")
}
