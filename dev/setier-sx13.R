# Lane setier: dev/setier-sx.R seed 13, form 7 (y ~ s(x0) + s(x1) +
# s(x2) + s(x3)): theta_2 is lost with the reference BLAS and kept with
# OpenBLAS. Prints the estimates and the Hessian rows of the thetas.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(13)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, scale = 2,
                                   verbose = FALSE))
fit <- frm(bf(y ~ s(x0) + s(x1) + s(x2) + s(x3)), data = d)
cat("theta:", format(fit$estimates$theta, digits = 8), "\n")
cat("logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")
h <- fit$cache$hessian_fixed
j <- grep("theta", names(fit$opt$par))
cat("theta row max:", format(apply(abs(h$H[j, , drop = FALSE]), 1, max),
                             digits = 3), "\n")
print(ns$sdr_of(fit)$se_lost)
