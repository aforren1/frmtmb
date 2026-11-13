# Short check of 0.82.0, probe 3: the asym_laplace warning (which
# parameter, what words), and zi ~ x seed 7 at new x on the side where
# the undetermined slope pushes zi up; seed 5's verdicts.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
q <- function(expr) withCallingHandlers(expr,
  warning = function(w) {
    cat("  W:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("  M:", substr(conditionMessage(m), 1, 300), "\n")
    invokeRestart("muffleMessage")
  })
set.seed(5)
x <- stats::rnorm(100)
dd <- data.frame(x = x, y = 1 + 0.5 * x + stats::rnorm(100))
cat("## asym_laplace y ~ x: the fit\n")
fa <- q(frm(y ~ x, data = dd, family = asym_laplace()))
print(fixef(fa)[, 1:2])
print(fa$cache$sdr$se_lost %||% "no sdr yet")
cat("## fitted()\n")
ft <- q(fitted(fa))
cat("## summary of fixed effects after fitted()\n")
print(q(fixef(fa))[, 1:2])
print(fa$cache$sdr$se_lost)
zib <- function(seed) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
  d
}
f7 <- q(frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = zib(7)))
cat("## zib seed 7, new x where zi_x pushes zi up (slope -0.32)\n")
for (x in c(-5, -10, -30, -60, -100)) {
  cat("x =", x, "\n")
  print(q(fitted(f7, newdata = data.frame(x = x, N = 10))))
}
cat("## zib seed 5: verdicts\n")
f5 <- q(frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(),
            data = zib(5)))
print(q(fixef(f5))[, 1:2])
print(f5$cache$sdr$se_lost)
cat("range of x:", range(zib(5)$x), "\n")
