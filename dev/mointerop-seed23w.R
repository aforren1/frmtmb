## The warnings of the review's degenerate (mo(x) | g) gap fit, seed 23
.libPaths(c("C:/Users/adf44/source/r/wt-mointerop-lib", "C:/Users/adf44/source/r/rellib-r14", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(23); n <- 300
d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n), g = factor(sample(1:10, n, TRUE)))
d$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d$g]) * d$x + 0.5 * d$z + rnorm(n)
w <- character(0)
fit <- withCallingHandlers(frm(y ~ mo(x) + z + (mo(x) | g), data = d), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat(substr(w, 1, 150), sep = "\n")
