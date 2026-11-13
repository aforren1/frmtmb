# Short check of lane arspat's fix round, part A: the student correction
# of a registered term with mi(). Exactness on the reviewer's prior case,
# the cases where the correction must be off, and what the fit says there.
source("dev/reviews/arspat-review-common.R")
n <- 60; m <- 20
M <- 0.5^abs(outer(1:n, 1:n, "-"))
cf <- function(nu, n, m) lgamma((nu + n) / 2) - lgamma((nu + n - m) / 2) +
  m / 2 * log(2 / (nu + n))
set.seed(4)
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)))
d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
miss <- sort(sample(n, m)); d$y[miss] <- NA
d$z <- 0.3 + 0.7 * ifelse(is.na(d$y), 1, d$y) + rnorm(n, 0, 0.5)
ok <- !is.na(d$y)
pr <- set_prior("gamma(2, 0.1)", class = "nu")

fit_w <- function(...) {
  warns <- character(0)
  f <- withCallingHandlers(frm(...), warning = function(w) {
    warns <<- c(warns, paste0(class(w)[1], ": ", substr(conditionMessage(w), 1, 90)))
    invokeRestart("muffleWarning")
  })
  attr(f, "warns") <- warns
  f
}
show_w <- function(f) cat("warnings:", if (length(attr(f, "warns")))
  paste(attr(f, "warns"), collapse = "\n  ") else "none", "\n")

chk("A1 prior gamma(2, 0.1), no random effect: now exact?", {
  f <- fit_w(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
             family = student(), prior = pr)
  dp <- ns$eval_dpars(f)[[1L]]; nu <- dp$nu[1]
  ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok], dp$sigma[1]^2 * M[ok, ok],
                      df = nu, log = TRUE)
  cat(sprintf("nu %.4g logLik - exact %+.3e (c would be %+.4f)\n", nu,
              as.numeric(logLik(f)) - ex, cf(nu, n, m)))
  cat("seam_t_missing:", unlist(ns$autocor_seam_t_missing(f$frame)),
      " mi_laplace_exact:", ns$mi_laplace_exact(f$frame), "\n")
  show_w(f)
  NULL
})
chk("A2 ML, no prior: still exact", {
  f <- fit_w(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
             family = student())
  dp <- ns$eval_dpars(f)[[1L]]; nu <- dp$nu[1]
  ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok], dp$sigma[1]^2 * M[ok, ok],
                      df = nu, log = TRUE)
  cat(sprintf("nu %.4g logLik - exact %+.3e\n", nu, as.numeric(logLik(f)) - ex))
  NULL
})
chk("A3 with (1 | h): correction off; is the error reported?", {
  f <- fit_w(y | mi() ~ x + (1 | h) + fcor(M), data = d,
             data2 = list(M = M), family = student(), prior = pr)
  cat("seam_t_missing:", length(ns$autocor_seam_t_missing(f$frame)),
      " mi_laplace_exact:", ns$mi_laplace_exact(f$frame),
      " tk:", format(f$mi_laplace$tk %||% NA), "\n")
  show_w(f)
  # the objective with the random-effect sd pushed to ~0 is the no-RE
  # model, whose exact margin is closed form: the error left there is c
  lab <- ns$outer_par_labels(f)
  cat("labels:", lab, "\n")
  p <- f$opt$par
  isd <- grep("^sd_|sd\\(|h__", lab)
  cat("sd index:", isd, "\n")
  if (length(isd) == 1L) {
    p[isd] <- -12
    lap <- -f$obj$fn(p)
    beta <- unname(p[which(ns$outer_par_map(f)$comp == "beta")])
    sig <- exp(unname(p[grep("^sigma", lab)]))
    nu <- 1 + exp(unname(p[grep("^nu", lab)]))
    mu <- beta[1] + beta[2] * d$x
    ex <- mvtnorm::dmvt(d$y[ok], mu[ok], sig^2 * M[ok, ok], df = nu, log = TRUE)
    cat(sprintf("at sd_h = exp(-12), nu %.4g: objective - exact %+.4f, c %+.4f\n",
                nu, lap - ex, cf(nu, n, m)))
  }
  cat("summary lines with Laplace/approx:",
      sum(grepl("Laplace|approx", capture.output(print(summary(f))))), "\n")
  NULL
})
chk("A4 student fcor y | mi() read by z ~ mi(y): correction must be off", {
  f <- fit_w(bf(y | mi() ~ x + fcor(M), family = student()) +
               bf(z ~ mi(y), family = gaussian()) + set_rescor(FALSE),
             data = d, data2 = list(M = M), prior = pr)
  cat("seam_t_missing:", length(ns$autocor_seam_t_missing(f$frame)),
      " mi_laplace_exact:", ns$mi_laplace_exact(f$frame),
      " tk:", format(f$mi_laplace$tk %||% NA), "\n")
  show_w(f)
  NULL
})
chk("A5 nu fixed: constant() prior", {
  f <- fit_w(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
             family = student(), prior = set_prior("constant(4)", class = "nu"))
  NULL
})
chk("A6 nu held by a tight prior normal(4, 0.001)", {
  f <- fit_w(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M),
             family = student(), prior = set_prior("normal(4, 0.001)", class = "nu"))
  dp <- ns$eval_dpars(f)[[1L]]; nu <- dp$nu[1]
  ex <- mvtnorm::dmvt(d$y[ok], rep_len(dp$mu, n)[ok], dp$sigma[1]^2 * M[ok, ok],
                      df = nu, log = TRUE)
  cat(sprintf("nu %.6g logLik - exact %+.3e\n", nu, as.numeric(logLik(f)) - ex))
  NULL
})
chk("A7 gaussian fcor + mi: unchanged (no correction)", {
  f <- frm(y | mi() ~ x + fcor(M), data = d, data2 = list(M = M))
  dp <- ns$eval_dpars(f)[[1L]]
  ex <- mvtnorm::dmvnorm(d$y[ok], rep_len(dp$mu, n)[ok], dp$sigma[1]^2 * M[ok, ok],
                         log = TRUE)
  cat(sprintf("logLik - exact %+.3e  seam_t_missing %d\n",
              as.numeric(logLik(f)) - ex,
              length(ns$autocor_seam_t_missing(f$frame))))
  NULL
})
