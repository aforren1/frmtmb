# Probe 4b: frmtmb.eam's wiener simulator given the boundary aterm, and
# what brms's posterior_predict_wiener() returns on brms's prep.
# Log: dev/bin2port-out/probe4b.txt. Seed 409.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.eam))
w <- brmsfamily("wiener")
set.seed(409)
s <- sim_draw(list(family = w, dpars = list(mu = rnorm(5), bs = rep(1.5, 5),
                                            ndt = rep(0.3, 5),
                                            bias = rep(0.5, 5)),
                   aterms = list(dec = c(1, 0, 1, 1, 0)), n = 5,
                   autocor = NULL, extra = NULL, max_iter = NULL))
cat("wiener sim_draw with dec: class", class(s), " dim", dim(s), "\n")
print(s)
ns <- 5
nobs <- 3
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                   bs = rchisq(ns, 3), ndt = rep(0.5, ns),
                   bias = rbeta(ns, 1, 1))
prep$data <- list(Y = abs(rnorm(ns)) + 0.5, dec = c(1, 0, 1))
b <- tryCatch(brms:::posterior_predict_wiener(2, prep),
              error = conditionMessage)
cat("brms posterior_predict_wiener: class", class(b), " dim", dim(b), "\n")
print(b)
