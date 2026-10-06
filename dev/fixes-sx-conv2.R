# Lane fixes, punch round: the one smooth fit that reports singular
# convergence after the reparameterization (seed 1, s(x1, by = g) + g).
#   Rscript dev/fixes-sx-conv2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(1)
d <- mgcv::gamSim(eg = 6, n = 200, scale = 2, verbose = FALSE)
d$z <- runif(200)
d$g <- factor(sample(c("a", "b", "c"), 200, TRUE))
fit <- suppressWarnings(frm(bf(y ~ s(x1, by = g) + g), data = d))
print(fit$opt[c("convergence", "message", "iterations")])
p <- fit$opt$par
print(round(cbind(par = p, grad = as.numeric(fit$obj$gr(p))), 6))
print(VarCorr(fit))
fe <- fixef(fit)
print(round(fe, 5))
for (rs in c(1, 3)) {
  f2 <- suppressWarnings(frm(bf(y ~ s(x1, by = g) + g), data = d,
                             control = frmtmb_control(restarts = rs)))
  cat("restarts", rs, ": conv", f2$opt$convergence, f2$opt$message,
      " logLik", format(as.numeric(logLik(f2)), digits = 12), "\n")
}
