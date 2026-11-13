# Review of lane arspat: does a fit of fcor() + mi() on student() whose
# nu comes out moderate say that its logLik carries the Laplace error?
# The error at the estimate is compared to the closed form
# lgamma((nu+n)/2) - lgamma((nu+n-m)/2) + m/2 log(2/(nu+n)).
source("dev/reviews/arspat-review-common.R")
n <- 60; m <- 20
M <- 0.5^abs(outer(1:n, 1:n, "-"))
for (seed in 1:12) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  # one multivariate-t draw, the model's own, df 2
  d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n)) *
    sqrt(2 / rchisq(1, 2))
  miss <- sort(sample(n, m))
  d$y[miss] <- NA
  warns <- character(0)
  f <- withCallingHandlers(
    frm(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
        family = student()),
    warning = function(w) {
      warns <<- c(warns, paste(class(w)[1], substr(conditionMessage(w), 1, 60)))
      invokeRestart("muffleWarning")
    })
  dp <- ns$eval_dpars(f)[[1L]]
  nu <- dp$nu[1]
  ok <- !is.na(d$y)
  ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok],
                      dp$sigma[1]^2 * M[ok, ok], df = nu, log = TRUE)
  cf <- lgamma((nu + n) / 2) - lgamma((nu + n - m) / 2) + m / 2 * log(2 / (nu + n))
  cat(sprintf("seed %2d nu %.4g laplace-exact %+.4f closed form %+.4f tk %s warnings: %s\n",
              seed, nu, as.numeric(logLik(f)) - ex, cf,
              format(f$mi_laplace$tk %||% NA, digits = 3),
              if (length(warns)) paste(warns, collapse = " | ") else "none"))
}
