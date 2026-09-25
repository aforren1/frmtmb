# Is check's examples NOTE on ?profile.frmtmb_fit (9.08 s user) the
# lane's doing? Time the example's body, interleaved per round, with an
# arithmetic control that must not move between arms.
#   PREDFIX_ARM=base|lane Rscript dev/predfix-extime.R
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
body <- function() {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  rownames(confint(fit))
  pr <- profile(fit, "theta_1")
  identical(profile(fit, "(Intercept)"), profile(fit, "Intercept"))
  confint(pr)
  prs <- profile(fit, c("x", "theta_1"))
  stopifnot(is.null(fit$par_units))
  invisible(NULL)
}
control <- function() sum(sort(runif(3e6)))
for (r in 1:3) {
  a <- system.time(body())[["user.self"]]
  b <- system.time(control())[["user.self"]]
  cat(sprintf("%s round %d: example %.2f s user, control %.2f s user\n",
              PREDFIX_ARM, r, a, b))
}
