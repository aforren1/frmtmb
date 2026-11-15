# Short check: the warning texts of the ALPHAX seed 11 regression fit.
.libPaths(c("C:/Users/adf44/source/r/wt-quadseam2-lib", "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
sn <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) + sqrt(1 - dl^2) * stats::rnorm(n)
}
set.seed(1011); stats::rnorm(135); stats::rnorm(60); sn(50, 1.5)
d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
d$y <- 0.5 * d$x + stats::rnorm(16, 0, 1)[d$g] + sn(80, 2)
withCallingHandlers(frm(bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x),
  family = skew_normal(), data = d, quadrature = TRUE,
  start = list(betad = c(0, 0, 0))),
  warning = function(w) { cat("W", class(w)[1], conditionMessage(w), "\n\n"); invokeRestart("muffleWarning") })
