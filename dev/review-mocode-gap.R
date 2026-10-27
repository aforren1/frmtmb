## Reviewer probe: ML predictions at an empty category inside the range
## (values 1, 2, 5; increasing response), seed 7.
##   Rscript dev/review-mocode-gap.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), arm, "\n")
set.seed(7)
d <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
d$y <- d$x + rnorm(200)
w <- character(0)
fit <- withCallingHandlers(frm(y ~ mo(x), data = d), warning = function(c) {
  w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
})
cat("fit warnings:", length(w), "\n")
for (m in w) cat("  ", substr(m, 1, 300), "\n")
mi <- fit$frame$linpreds[[1]]$mo[[1]]
cat("D =", mi$D, " min =", format(mi$min), "\n")
w <- character(0)
p <- withCallingHandlers(fitted(fit, newdata = data.frame(x = 1:5)),
                         warning = function(c) {
  w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
})
print(signif(p, 6))
cat("predict warnings:", length(w), "\n")
for (m in w) cat("  ", substr(m, 1, 300), "\n")
