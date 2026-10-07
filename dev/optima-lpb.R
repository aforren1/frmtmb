# Lane optima: after mo_search(), does the objective's remembered best
# point (last.par.best, which sdreport() and the joint density read)
# still sit at the fit's optimum? test-brms-likelihood.R's check C row
# 3 variant, y ~ mo(inc):z + (1 | g).
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(3)
dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300),
                 g = factor(rep(1:20, 15)))
dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
fit <- frm(bf(y ~ mo(inc):z + (1 | g)) + gaussian(), data = dm,
           verbose = TRUE)
r <- fit$obj$env$random
lpb <- fit$obj$env$last.par.best[-r]
print(rbind(opt = fit$opt$par, last.par.best = lpb))
cat("objective at opt", format(fit$opt$objective, digits = 12),
    "value.best", format(fit$obj$env$value.best, digits = 12), "\n")
print(fit$opt$mo_search)
print(frmtmb:::mo_simplex(fit$estimates$zeta1))
