# Why brms 2.23.0's negbinomial2 fit of dev/nb2dm-brmsfit.R (case nb2)
# did not mix: per-chain summaries, and Stan's GLM primitive
# neg_binomial_2_log_glm_lpmf against neg_binomial_2_log_lpmf at a draw
# of the stuck chain.
#   Rscript dev/nb2dm-brms-stuck.R > dev/nb2dm-out/brms-stuck.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
b <- readRDS("dev/nb2dm-out/brmsfit-nb2-glmbug.rds")
dr <- posterior::as_draws_df(b)
cat("per-chain means:\n")
print(aggregate(cbind(b_Intercept, b_x, sigma, lp__) ~ .chain, data = dr,
                FUN = mean), digits = 5)
cat("per-chain divergences:\n")
np <- nuts_params(b)
print(tapply(np$Value[np$Parameter == "divergent__"],
             np$Chain[np$Parameter == "divergent__"], sum))
d <- b$data
code <- "functions {
  real glm_ll(array[] int y, matrix X, real a, vector beta, real sigma) {
    return neg_binomial_2_log_glm_lpmf(y | X, a, beta, inv(sigma));
  }
  real plain_ll(array[] int y, vector eta, real sigma) {
    return neg_binomial_2_log_lpmf(y | eta, inv(sigma));
  }
}
"
rstan::expose_stan_functions(rstan::stanc(model_code = code))
X <- matrix(d$x - mean(d$x), ncol = 1)
pick <- function(ch) {
  s <- dr[dr$.chain == ch, ]
  s[nrow(s), ]
}
for (ch in 1:4) {
  s <- pick(ch)
  a <- s$b_Intercept + s$b_x * mean(d$x)   # brms's centered Intercept
  eta <- a + X[, 1] * s$b_x
  g <- glm_ll(as.integer(d$y), X, a, s$b_x, s$sigma)
  p <- plain_ll(as.integer(d$y), eta, s$sigma)
  r <- sum(dnbinom(d$y, size = 1 / s$sigma, mu = exp(eta), log = TRUE))
  cat(sprintf("chain %d last draw: b_Intercept %.4g b_x %.4g sigma %.4g | glm %.10g plain %.10g dnbinom %.10g\n",
              ch, s$b_Intercept, s$b_x, s$sigma, g, p, r))
}
