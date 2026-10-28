# What brms 2.23.0's own helpers compute on the preps of the ported
# blocks "log_lik for ARMA models runs without errors" (log_lik:130,
# 133) and "posterior_predict for ARMA covariance models runs without
# errors" (posterior_predict:93, 97): the preps carry no acframe, so is
# the correlation brms builds the ARMA one or the identity, and does
# prep$data$se enter?
# Usage: Rscript dev/sigmainn-bin2-brmsprep.R
# Log: dev/sigmainn-out/bin2-brmsprep.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(130)
ns <- 20
nobs <- 15
prep <- structure(list(ndraws = ns), class = "brmsprep")
prep$dpars <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                   sigma = rchisq(ns, 3), nu = rgamma(ns, 5) + 15)
prep$ac <- list(ar = matrix(rbeta(ns, 0.5, 0.5), ncol = 1),
                ma = matrix(rbeta(ns, 0.2, 1), ncol = 1),
                begin_tg = 2, end_tg = 5)
prep$data <- list(Y = rnorm(nobs), se = rgamma(ns, 10))
cat("has_ac_class(NULL, 'arma'):",
    brms:::has_ac_class(prep$ac$acframe, "arma"), "\n")
S <- brms:::get_cov_matrix_ac(prep, 2:5)
cat("draw 1 covariance brms builds:\n")
print(S[1, , ])
cat("sigma^2 of draw 1:", prep$dpars$sigma[1]^2, "; se[2:5]^2:",
    prep$data$se[2:5]^2, "\n")
cat("off-diagonal max:", max(abs(S[1, , ][upper.tri(S[1, , ])])), "\n")
