# Review: an integer mo() variable with two values (1 and 2), 0.75.0
# against 0.76.0: fitted, or refused?
#   Rscript dev/review-rel076-mo2.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(2)
d <- data.frame(k = sample(1:2, 100, TRUE)); d$y <- 0.5 * d$k + rnorm(100)
r <- tryCatch({f <- frm(y ~ mo(k), data = d); paste("FITTED, logLik", format(as.numeric(logLik(f))), "| se of bsp", format(sqrt(diag(vcov(f)))[2]))},
  error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("frmtmb", format(packageVersion("frmtmb")), ":", r, "\n")
