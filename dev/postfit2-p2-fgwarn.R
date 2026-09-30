.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(47)
d <- data.frame(x = stats::runif(200), g = factor(rep(1:5, 40)))
d$y <- sin(2 * pi * d$x) + stats::rnorm(5, 0, 0.5)[d$g] + stats::rnorm(200, 0, 0.3)
w <- NULL
fg <- withCallingHandlers(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | g)), family = gaussian(), data = d),
  warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
print(w)
