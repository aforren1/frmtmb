# Reviewer of lane optstart, round 2: test-v14.R's asym_laplace fit
# (quantile 0.25, seed 13, n 500) against its exact maximum: mu from
# quantreg::rq(tau = 0.25), sigma = the mean check loss of its
# residuals, logLik = n (log p + log(1 - p) - log sigma - 1).
#   Rscript dev/optstart-rev2-v14.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
set.seed(13)
n <- 500
x <- rnorm(n)
dd <- data.frame(x = x, y = 1 + 0.5 * x + rnorm(n) * (1 + 0.4 * x^2))
w <- character()
fq <- withCallingHandlers(
  frm(bf(y ~ x, quantile = 0.25) + asym_laplace(), data = dd),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
p <- 0.25
rq <- quantreg::rq(y ~ x, tau = p, data = dd)
r <- residuals(rq)
s <- mean(0.5 * (abs(r) + (2 * p - 1) * r))
ll <- n * (log(p) + log(1 - p) - log(s) - 1)
cat(sprintf("frmtmb logLik %.10f code %d (%s), %d warnings: %s\n",
            as.numeric(logLik(fq)), fq$opt$convergence, fq$opt$message,
            length(w), paste(substr(w, 1, 60), collapse = " | ")))
cat(sprintf("exact maximum %.10f; frmtmb short by %.3g\n", ll,
            ll - as.numeric(logLik(fq))))
cat("mu frmtmb", format(unname(fixef_by_dpar(fq)$mu), digits = 10),
    " rq", format(unname(coef(rq)), digits = 10), "\n")
