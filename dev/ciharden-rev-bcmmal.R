# Reviewer: what a user sees from the Malingering_2 fit, and the test's
# new tolerance at the fitted point.
.libPaths(c("C:/Users/adf44/source/r/wt-ciharden-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
suppressPackageStartupMessages(library(testthat))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden/tests/testthat"
env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-stan.R", "helper-bcm.R", "test-bcm-latent-mixtures.R")) {
  if (file.exists(file.path(wt, f))) {
    ex <- parse(file.path(wt, f))
    for (e in ex) {
      if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
    }
  }
}
d <- env$bcm_malingering_data()
w <- character()
fit <- withCallingHandlers(
  frm(bf(k | trials(n) ~ 1),
      family = mixture(beta_binomial, beta_binomial), data = d,
      start = env$bcm_mix_start(c(0.55, 0.95), phi = c(15, 15), theta = 0)),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
cat("BLAS threads", Sys.getenv("OPENBLAS_NUM_THREADS"), "R", R.home(), "\n")
p <- env$bcm_betamix_pars(fit)
cat(sprintf("phi1 %.6g phi2 %.6g mu2 %.10f\n", p$phi1, p$phi2, p$mu2))
cat("warnings:", length(w), "\n")
for (x in w) cat("  -", substr(gsub("\n", " ", x), 1, 300), "\n")
size <- function(mu, phi) {
  a <- mu * phi
  b <- (1 - mu) * phi
  sum(abs(lgamma(d$k + a)) + abs(lgamma(45 - d$k + b)) +
        abs(lgamma(a + b + 45)) + abs(lgamma(a)) + abs(lgamma(b)) +
        abs(lgamma(a + b)))
}
S <- size(p$mu1, p$phi1) + size(p$mu2, p$phi2)
lp <- env$frm_joint_lp(fit)
cat(sprintf("S %.4g lp %.6g old bound %.4g new bound %.4g\n", S, lp,
            1e-6 * max(1, abs(lp)), max(1e-6 * max(1, abs(lp)),
                                       4 * .Machine$double.eps * S)))
sm <- tryCatch(summary(fit), error = function(e) NULL)
cf <- tryCatch(fixef(fit), error = function(e) NULL)
print(tryCatch(sqrt(diag(vcov(fit))), error = function(e) conditionMessage(e)))
print(tryCatch(fit$cache$se_explained, error = function(e) NULL))
