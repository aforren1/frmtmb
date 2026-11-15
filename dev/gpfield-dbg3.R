# Lane gpfield: gsigma step 40 from the constant-sigma model's maximum.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
r <- gf_row("ng:gsigma:step:40:1")
f0 <- gf_fit(y ~ gp(x), data = r$data, family = gaussian())$fit
cat("constant sigma ll", as.numeric(logLik(f0)), "theta", f0$estimates$theta,
    "betad", f0$estimates$betad, "\n")
f1 <- gf_fit(r$fml, data = r$data, family = r$family)$fit
cat("default ll", as.numeric(logLik(f1)), "pars", f1$opt$par, "\n")
q <- gf_fit(r$fml, data = r$data, family = r$family,
            start = list(theta = f0$estimates$theta))
cat("from its theta ll", as.numeric(logLik(q$fit)), "pars", q$fit$opt$par,
    "nwarn", length(q$warn), "\n")
bd <- f1$estimates$betad
bd[1] <- f0$estimates$betad[1]; bd[2] <- 0
q2 <- gf_fit(r$fml, data = r$data, family = r$family,
             start = list(theta = f0$estimates$theta, betad = bd,
                          beta = f0$estimates$beta))
cat("from its theta, beta, sigma ll", as.numeric(logLik(q2$fit)), "pars",
    q2$fit$opt$par, "nwarn", length(q2$warn), "\n")
q3 <- gf_fit(r$fml, data = r$data, family = r$family,
             start = list(theta = c(log(0.7552), log(0.7996))))
cat("from the better theta ll", as.numeric(logLik(q3$fit)), "pars",
    q3$fit$opt$par, "\n")
