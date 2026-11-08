# A lognormal design with me() on mu and sigma whose fit converges cleanly,
# for the test of the not-finite mean. Seeds 9203 to 9212.
source("dev/readers-setup.R")
for (seed in 9203:9212) {
  set.seed(seed)
  n <- 600
  d <- data.frame(zo = rnorm(n), sz = 0.5)
  zl <- d$zo + rnorm(n, 0, 0.5)
  d$y <- exp(0.5 + 0.4 * zl + rnorm(n, 0, exp(-0.5 + 0.5 * zl)))
  w <- character()
  fit <- withCallingHandlers(
    frm(bf(y ~ me(zo, sz), sigma ~ me(zo, sz)), data = d,
        family = lognormal()),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  nd <- data.frame(zo = c(0, 1.5), sz = 0.5)
  f <- suppressWarnings(fitted(fit, newdata = nd))
  cat(seed, "warnings:", length(w), " Estimate:", f[, "Estimate"],
      " Q:", signif(f[, "Q2.5"], 3), signif(f[, "Q97.5"], 3), "\n")
}
