# The vigport row "brms's bf(y ~ 1, mu1 ~ x1, mu2 ~ x1) for a mixture"
# on the lane build (argument: "before" reads rellib-r8 alone).
# Output: dev/famsimplex-mixmu1-<arm>.log
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (!identical(arm, "before")) "C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3)
n <- 400
dd <- data.frame(x1 = rnorm(n))
k <- rbinom(n, 1, 0.4)
dd$y <- ifelse(k == 1, rnorm(n, -2 + 0.5 * dd$x1, 1), rnorm(n, 2 - 0.3 * dd$x1, 1))
r <- tryCatch({
  f1 <- frm(bf(y ~ 1, mu1 ~ x1, mu2 ~ x1), family = mixture(gaussian(), gaussian()),
            data = dd)
  f0 <- frm(bf(y ~ x1), family = mixture(gaussian(), gaussian()), data = dd)
  c(logLik_mu1_mu2 = as.numeric(logLik(f1)), logLik_main = as.numeric(logLik(f0)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
