# Lane parity2, fix round: which fit of test-ce-default-effects.R warns
# under OpenBLAS (run with the OpenBLAS R copy). Data as the test's.
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- 150
d <- data.frame(x = stats::runif(n, 0.1, 2), z = stats::rnorm(n),
                k = sample(1:4, n, TRUE), w = stats::rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)))
d$y <- 1 + d$x * d$z + 0.3 * d$k + stats::rnorm(n)
d$f <- factor(rep(c("a", "b", "c"), length.out = n))
fs <- list(bf(y ~ exp(x):z), bf(y ~ t2(x, z)), bf(y ~ mo(k):z),
           bf(y ~ x, sigma ~ exp(z)), bf(y ~ x + (1 + z | g), sigma ~ w),
           bf(y ~ x * z * w), bf(y ~ poly(x, 2) + z + x:w),
           bf(y ~ a * exp(b * x), a ~ 1 + w, b ~ z, nl = TRUE),
           bf(y ~ f * x))
for (f in fs) {
  w <- character()
  withCallingHandlers(
    suppressMessages(conditional_effects(frm(f, family = gaussian(),
                                             data = d), resolution = 5)),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  cat(deparse1(f$formula), " pforms:", length(f$pforms), " warnings:",
      length(w), if (length(w)) substr(w[1], 1, 90), "\n")
}
