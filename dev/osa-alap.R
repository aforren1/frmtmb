# Lane osa: is the 2.5e-6 the audit reads on asym_laplace() the
# reference integral or the residual? The audit integrates the density
# from -Inf to y in one call, across the kink at mu. Here the same fit
# (dev/osa-audit.R's real-line design, data seed 9, n = 150) is compared
# with qnorm(F(y)) from the closed-form CDF of brms's asymmetric Laplace,
# F(y) = p exp((1 - p)(y - mu) / sigma) for y < mu and
# 1 - (1 - p) exp(-p (y - mu) / sigma) above.
#
#   Rscript dev/osa-alap.R
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(9)
n <- 150
x <- rnorm(n)
zi <- runif(n)  # the audit draws its zero-inflation flags here
y <- 1 + 0.7 * x + rnorm(n)
fit <- frm(y ~ x, family = asym_laplace(), data = data.frame(x = x, y = y))
dp <- lapply(frmtmb:::eval_dpars(fit)[[1]], as.numeric)
p <- dp$quantile
z <- (y - dp$mu) / dp$sigma
Fy <- ifelse(z < 0, p * exp((1 - p) * z), 1 - (1 - p) * exp(-p * z))
r <- residuals(fit, type = "osa")[, "Estimate"]
cat("BEGIN osa-alap\n")
cat(sprintf("max|r - qnorm(F(y))| closed form: %.3g\n",
            max(abs(r - qnorm(Fy)))))
cat("END osa-alap\n")
