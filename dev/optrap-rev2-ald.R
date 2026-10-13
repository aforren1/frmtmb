# Reviewer of lane optrap, re-check: test-v14.R's asym_laplace fit
# (quantile 0.25, seed 13), the eleventh worse fit of the lane's fit
# diff. Its exact maximum: the betas of quantreg::rq() and sigma the
# mean check loss, for the brms density tau (1 - tau) / s exp(-rho / s).
#   Rscript dev/optrap-rev2-ald.R <lib or none>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(13)
n <- 500
x <- rnorm(n)
dd <- data.frame(x = x, y = 1 + 0.5 * x + rnorm(n) * (1 + 0.4 * x^2))
fq <- suppressWarnings(frm(bf(y ~ x, quantile = 0.25) + asym_laplace(),
                           data = dd))
rq <- quantreg::rq(y ~ x, tau = 0.25, data = dd)
r <- residuals(rq)
tau <- 0.25
rho <- r * (tau - (r < 0))
s <- mean(rho)
ll_ex <- sum(log(tau * (1 - tau)) - log(s) - rho / s)
cat("lib", dirname(find.package("frmtmb")), "| code", fq$opt$convergence,
    "| logLik", format(as.numeric(logLik(fq)), digits = 10), "| exact",
    format(ll_ex, digits = 10), "| gap", signif(as.numeric(logLik(fq)) - ll_ex, 3),
    "| units", if (is.null(fq$par_units)) "none" else signif(fq$par_units, 3),
    "\n")
