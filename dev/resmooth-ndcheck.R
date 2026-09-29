# Lane wt-resmooth. Is predict(newdata = the fitted rows) equal to the
# in-sample prediction? Probed because resmooth-before.txt showed a
# max|nd-ins| of 0.101 on y ~ s(x), which would be a defect of its own.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)

chk <- function(lab, f) {
  fit <- suppressWarnings(frm(f, data = d))
  a <- fitted(fit, re_formula = NA)[, "Estimate"]
  b <- predict(fit, newdata = d, re_formula = NA)[, "Estimate"]
  c2 <- predict(fit, re_formula = NA)[, "Estimate"]
  a1 <- fitted(fit)[, "Estimate"]
  b1 <- predict(fit, newdata = d)[, "Estimate"]
  cat(sprintf("%-14s NA: fitted vs predict(nd) %9.3g | fitted vs predict() %9.3g | NULL: %9.3g\n",
              lab, max(abs(a - b)), max(abs(a - c2)), max(abs(a1 - b1))))
}
chk("y ~ x", bf(y ~ x))
chk("y ~ s(x)", bf(y ~ s(x)))
chk("y ~ s(x)+(1|g)", bf(y ~ s(x) + (1 | g)))
