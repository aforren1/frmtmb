# Minimal repro: Stan's neg_binomial_2_log_glm_lpmf against
# neg_binomial_2_log_lpmf and dnbinom at a large precision phi, one
# observation, for dev/upstream-bugs.md.
#   Rscript dev/nb2dm-stan-glm-bug.R > dev/nb2dm-out/stan-glm-bug.txt 2>&1
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
cat("rstan", format(packageVersion("rstan")), "StanHeaders",
    format(packageVersion("StanHeaders")), "\n")
code <- "functions {
  real glm_ll(array[] int y, matrix X, real a, vector beta, real phi) {
    return neg_binomial_2_log_glm_lpmf(y | X, a, beta, phi);
  }
  real plain_ll(array[] int y, vector eta, real phi) {
    return neg_binomial_2_log_lpmf(y | eta, phi);
  }
}
"
rstan::expose_stan_functions(rstan::stanc(model_code = code))
X <- matrix(0, 1, 1)
for (y in c(0L, 1L, 5L)) {
  for (eta in c(0, 3)) {
    for (phi in 10^c(4, 8, 12, 16, 20, 23)) {
      g <- glm_ll(y, X, eta, 0, phi)
      p <- plain_ll(y, eta, phi)
      r <- dnbinom(y, size = phi, mu = exp(eta), log = TRUE)
      cat(sprintf("y %d eta %g phi %.0e: glm %.12g plain %.12g dnbinom %.12g\n",
                  y, eta, phi, g, p, r))
    }
  }
}
