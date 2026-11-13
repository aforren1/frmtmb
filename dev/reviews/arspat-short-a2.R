# Short check part A2: with (1 | h) the correction is off. How large is
# the error left, and does any check report it? The RE sd at its boundary
# (0) gives the no-RE model, whose margin is closed form; a real RE is
# checked against brute-force quadrature over the 6 RE levels is out of
# reach, so the second case reads tk and the warnings only.
source("dev/reviews/arspat-review-common.R")
n <- 60; m <- 20
M <- 0.5^abs(outer(1:n, 1:n, "-"))
cf <- function(nu, n, m) lgamma((nu + n) / 2) - lgamma((nu + n - m) / 2) +
  m / 2 * log(2 / (nu + n))
mk <- function(re_sd, seed = 4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)))
  d$y <- 1 + 0.5 * d$x + rnorm(6, 0, re_sd)[d$h] +
    as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
  miss <- sort(sample(n, m)); d$y[miss] <- NA
  d$z <- 0.3 + 0.7 * ifelse(is.na(d$y), 1, d$y) + rnorm(n, 0, 0.5)
  d
}
pr <- set_prior("gamma(2, 0.1)", class = "nu")
for (re_sd in c(0, 1.5)) {
  chk(paste("(1 | h), data RE sd", re_sd), {
    d <- mk(re_sd)
    ok <- !is.na(d$y)
    warns <- character(0)
    f <- withCallingHandlers(
      frm(y | mi() ~ x + (1 | h) + fcor(M), data = d, data2 = list(M = M),
          family = student(), prior = pr),
      warning = function(w) {
        warns <<- c(warns, paste(class(w)[1], substr(conditionMessage(w), 1, 80)))
        invokeRestart("muffleWarning")
      }, message = function(mm) invokeRestart("muffleMessage"))
    dp <- ns$eval_dpars(f)[[1L]]; nu <- dp$nu[1]
    vc <- tryCatch(VarCorr(f), error = function(e) NULL)
    sdh <- tryCatch(as.numeric(vc$h$sd[1, 1]), error = function(e) NA)
    cat(sprintf("nu %.4g sd_h %s tk %s mi_laplace_exact %s\n", nu,
                format(sdh, digits = 4), format(f$mi_laplace$tk %||% NA, digits = 4),
                ns$mi_laplace_exact(f$frame)))
    cat("warnings:", if (length(warns)) paste(warns, collapse = " | ") else "none", "\n")
    # with sd_h at 0 the fit is the no-RE model: exact margin closed form
    mu0 <- fixef(f)[1, 1] + fixef(f)[2, 1] * d$x
    ex <- mvtnorm::dmvt(d$y[ok], mu0[ok], dp$sigma[1]^2 * M[ok, ok], df = nu,
                        log = TRUE)
    cat(sprintf("logLik - exact(no-RE margin at these params) %+.4f, c %+.4f\n",
                as.numeric(logLik(f)) - ex, cf(nu, n, m)))
    NULL
  })
}
chk("student fcor y | mi() read by z ~ mi(y), resp prior", {
  d <- mk(0)
  f <- frm(bf(y | mi() ~ x + fcor(M), family = student()) +
             bf(z ~ mi(y), family = gaussian()) + set_rescor(FALSE),
           data = d, data2 = list(M = M),
           prior = set_prior("gamma(2, 0.1)", class = "nu", resp = "y"))
  cat("seam_t_missing:", length(ns$autocor_seam_t_missing(f$frame)),
      " mi_laplace_exact:", ns$mi_laplace_exact(f$frame),
      " tk:", format(f$mi_laplace$tk %||% NA), "\n")
  NULL
})
