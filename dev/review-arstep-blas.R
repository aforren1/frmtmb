# Review of lane arstep: the same fits under whichever BLAS this R carries.
# Prints, per row, the log-likelihood, whether the AR marginal was used,
# the second descent's record and the estimates, to 10 digits, for a diff
# across reference BLAS, OpenBLAS 0.3.26 and 0.3.32 (both with LAPACK).
# Usage: <Rscript of the build> review-arstep-blas.R <tag>
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
tag <- commandArgs(TRUE)[1]
rv_libs("lane")
cat("BLAS", tag, sum(chol2inv(chol(crossprod(matrix(sin(1:400), 20)) + diag(20)))), "\n",
    file = stdout())
cat("blas sig", sprintf("%.17g", sum(chol2inv(chol(crossprod(matrix(sin(1:400), 20)) + diag(20))))), "\n")
step_ar <- function(n, seed, phi = 0.4) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  data.frame(x = x, t = seq_len(n),
             y = 2 * (x > 5) - 1 + as.numeric(arima.sim(list(ar = phi), n, sd = 0.03)))
}
ng <- function(case, cn, n) {
  curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x,
                 two = function(x) sin(x) + 0.4 * sin(10 * x),
                 step = function(x) 2 * (x > 5) - 1)
  set.seed(74000 + 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- curves[[cn]](x)
  switch(case,
    bern = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))),
    pois_mid = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))),
    binom5 = data.frame(x = x, tr = 5, y = rbinom(n, 5, plogis(eta - mean(eta)))),
    gamma = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))))
}
rows <- list(
  list("gar step 160 s1", y ~ gp(x) + ar(time = t), step_ar(160, 1), gaussian()),
  list("gar step 80 s22", y ~ gp(x) + ar(time = t), step_ar(80, 22), gaussian()),
  list("gar step 80 s23", y ~ gp(x) + ar(time = t), step_ar(80, 23), gaussian()),
  list("gar step 160 s3", y ~ gp(x) + ar(time = t), step_ar(160, 3), gaussian()),
  list("bern mix 40", y ~ gp(x), ng("bern", "mix", 40), bernoulli()),
  list("bern step 40", y ~ gp(x), ng("bern", "step", 40), bernoulli()),
  list("binom5 mix 40", y | trials(tr) ~ gp(x), ng("binom5", "mix", 40), binomial()),
  list("binom5 two 40", y | trials(tr) ~ gp(x), ng("binom5", "two", 40), binomial()),
  list("pois_mid two 80", y ~ gp(x), ng("pois_mid", "two", 80), poisson()),
  list("gamma two 40", y ~ gp(x), ng("gamma", "two", 40), Gamma(link = "log")))
for (r in rows) {
  q <- rv_quiet(frm(r[[2]], data = r[[3]], family = r[[4]]))
  fit <- q$value
  if (inherits(fit, "error")) { cat("BL", tag, r[[1]], "ERROR", conditionMessage(fit), "\n"); next }
  e <- fit$estimates
  cat(sprintf("BL %s %-16s ll %.10f conv %d w %d marg %s second %s theta %s ac %s\n",
              tag, r[[1]], as.numeric(logLik(fit)), fit$opt$convergence, length(q$warnings),
              !is.null(fit$obj$env$.frmtmb_gp_marginal),
              if (is.null(fit$opt$gp_second)) "-" else paste(sprintf("%.6g", fit$opt$gp_second), collapse = ","),
              paste(sprintf("%.8f", e$theta), collapse = ","),
              if (length(e$thetaac)) paste(sprintf("%.8f", e$thetaac), collapse = ",") else "-"))
  flush.console()
}
cat("DONE\n")
