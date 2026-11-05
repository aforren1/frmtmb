# Item 1: exp(0.8x), seed 5701, 100 rows. The z tape started at the closed form's maximum,
# with the intercept profiled. Does it stay there?
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
set.seed(5701); n <- 100; x <- sort(runif(n, 0, 10)); e0 <- rnorm(n); y <- exp(0.8 * x) + 0.3 * e0
s <- cf_setup(y, x)
m <- cf_max(s, ngrid = 30L)
cat(sprintf("cf ll %.5f 256 bits %.5f mu %.1f sd %.1f rho %.4f sigma %.5f\n", m$ll, cf_ll_mpfr(s, m$beta, m$sigma, m$sd, m$rho), m$beta, m$sd, m$rho, m$sigma))
r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x), start = list(theta = c(log(m$sd), log(m$rho)), betad = log(m$sigma), beta = m$beta[1]))
g <- rev_gp_read(r$fit)
cat(sprintf("from cf: tape ll %.5f 256 bits %.5f mu %.1f sd %.1f rho %.4f sigma %.5f conv %d w %d\n", g$ll, cf_ll_mpfr(s, g$beta, g$sigma, g$sd, g$rho), g$beta, g$sd, g$rho, g$sigma, g$conv, length(r$warn)))
if (length(r$warn)) cat(substr(r$warn[1], 1, 200), "\n")
cat(sprintf("tape objective at the cf point, before any step: %.5f\n", -r$fit$obj$fn(c(m$beta, log(m$sigma), log(m$sd), log(m$rho)))))
