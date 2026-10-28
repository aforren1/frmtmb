# frmtmb.sample test-loglik-autocor.R: sigma is the innovation sd under
# ar(), ma() and arma() with cov = TRUE, as in brms.
f <- "extensions/frmtmb.sample/tests/testthat/test-loglik-autocor.R"
x <- paste(readLines(f), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop("not once (", n, "): ", old)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("  acf_of <- function(a, m) function(L) {
    v <- stats::ARMAacf(ar = a[a != 0], ma = m[m != 0], lag.max = max(L))
    matrix(v[L + 1], nrow(L))
  }",
"  # the autocovariance at unit innovation variance: the correlation of
  # stats::ARMAacf() times the process variance from the psi weights,
  # since sigma is the innovation sd, as in brms
  acf_of <- function(a, m) function(L) {
    v <- stats::ARMAacf(ar = a[a != 0], ma = m[m != 0], lag.max = max(L))
    psi <- c(1, stats::ARMAtoMA(ar = a[a != 0], ma = m[m != 0],
                                lag.max = 5000))
    matrix(v[L + 1] * sum(psi^2), nrow(L))
  }")
rep1("  # brms's sigma is the innovation sd and frmtmb's the marginal one;
  # g0 is the ARMA(1, 1) variance at unit innovations
  g0 <- function(a, m) (1 + 2 * a * m + m^2) / (1 - a^2)
", "  # brms's sigma is the innovation sd, and so is frmtmb's: the draws'
  # sigma goes to brms unchanged
")
rep1("                  sigma = sg[s] / sqrt(g0(cs$a[s], cs$m[s])))",
     "                  sigma = sg[s])")
rep1("  # sigma ~ x: brms's log sigma is frmtmb's plus log sqrt(1 - ar^2)",
     "  # sigma ~ x: brms's log sigma, the log innovation sd, is frmtmb's")
rep1("                Intercept_sigma = bs0[s] + bs1[s] * mean(d$x) +
                  0.5 * log(1 - a[s]^2))",
     "                Intercept_sigma = bs0[s] + bs1[s] * mean(d$x))")
rep1("test_that(\"ar(cov = TRUE) under the default priors: log_lik() of the
           draws reads sigma as the marginal sd\", {
  skip_sampler()
  # the sigma prior is read at the innovation sd, sigma sqrt(1 - ar^2),
  # but the draws' sigma column is frmtmb's marginal sd, and the
  # conditional is of the marginal covariance sigma^2 ar^|lag|",
"test_that(\"ar(cov = TRUE) under the default priors: log_lik() of the
           draws reads sigma as the innovation sd\", {
  skip_sampler()
  # the draws' sigma column is the innovation sd, as brms's is, and the
  # conditional is of the covariance sigma^2 ar^|lag| / (1 - ar^2).
  # frmtmb 0.76.0's sigma was the marginal sd, with covariance
  # sigma^2 ar^|lag|; this block pinned that reading")
rep1("                       rep(m[s, \"sigma\"], nrow(d)), function(L) a1^L)",
     "                       rep(m[s, \"sigma\"], nrow(d)),
                       function(L) a1^L / (1 - a1^2))")
writeLines(x, f)
