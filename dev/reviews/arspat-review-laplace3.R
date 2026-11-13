# Review of lane arspat: fcor() + mi() on student() with nu pinned by a
# constant prior: logLik against the exact margin and the closed form.
source("dev/reviews/arspat-review-common.R")
n <- 60; m <- 20
M <- 0.5^abs(outer(1:n, 1:n, "-"))
set.seed(4)
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
miss <- sort(sample(n, m)); d$y[miss] <- NA
for (pr in c("constant(4)", "gamma(2, 0.1)")) chk(pr, {
  f <- frm(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
           family = student(), prior = set_prior(pr, class = "nu"))
  dp <- ns$eval_dpars(f)[[1L]]; nu <- dp$nu[1]; ok <- !is.na(d$y)
  ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok], dp$sigma[1]^2 * M[ok, ok],
                      df = nu, log = TRUE)
  cf <- lgamma((nu + n) / 2) - lgamma((nu + n - m) / 2) + m / 2 * log(2 / (nu + n))
  cat(sprintf("nu %.4g logLik-exact %+.4f closed form %+.4f tk %s\n", nu,
              as.numeric(logLik(f)) - ex, cf, format(f$mi_laplace$tk %||% NA)))
  NULL
})
