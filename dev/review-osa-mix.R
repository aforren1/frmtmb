# Review lane osa: mixture() families now take the "cdf" route by default;
# the lane's audit did not fit one. Reference: the two-component CDF from
# ppois()/pnorm() at the fitted dpars, TMB's uniforms (set.seed(123)).
#   Rscript dev/review-osa-mix.R
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-mix\n")
set.seed(501)
n <- 160
z <- runif(n) < 0.4
d <- data.frame(cnt = ifelse(z, rpois(n, 1), rpois(n, 8)),
                yc = ifelse(z, rnorm(n, -2), rnorm(n, 2)))
U <- local({ set.seed(123); runif(n) })
fp <- frm(bf(cnt ~ 1), family = mixture(poisson(), poisson()), data = d)
dp <- lapply(frmtmb:::eval_dpars(fp)[[1]], as.numeric)
cat("dpars:", paste(names(dp), collapse = " "), "\n")
th1 <- plogis(dp$theta1)  # softmax against theta2 = 0 (summary: 0.39)
Fm <- function(q) th1 * ppois(q, dp$mu1) + (1 - th1) * ppois(q, dp$mu2)
pm <- function(q) th1 * dpois(q, dp$mu1) + (1 - th1) * dpois(q, dp$mu2)
r <- residuals(fp, type = "osa")[, "Estimate"]
ref <- qnorm(Fm(d$cnt) - U * pm(d$cnt))
cat(sprintf("mixture(poisson, poisson): max|r - ref| %.3g\n", max(abs(r - ref))))
fg <- frm(bf(yc ~ 1), family = mixture(gaussian(), gaussian()), data = d)
dg <- lapply(frmtmb:::eval_dpars(fg)[[1]], as.numeric)
t1 <- plogis(dg$theta1)
rg <- residuals(fg, type = "osa")[, "Estimate"]
refg <- qnorm(t1 * pnorm(d$yc, dg$mu1, dg$sigma1) +
                (1 - t1) * pnorm(d$yc, dg$mu2, dg$sigma2))
cat(sprintf("mixture(gaussian, gaussian): max|r - ref| %.3g\n", max(abs(rg - refg))))
cat("END review-osa-mix\n")
