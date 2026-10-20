# drop1() on an ordinary gaussian fit with log(x2), on rellib-r9: does
# the refit find x2? Output: dev/catmulti-out/drop1-base.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
f <- function() {
  set.seed(1)
  d <- data.frame(x = rnorm(100), x2 = exp(rnorm(100)))
  d$y <- d$x + log(d$x2) + rnorm(100)
  fit <- frm(bf(y ~ x + log(x2)), data = d)
  print(tryCatch(drop1(fit), error = function(e) conditionMessage(e)))
  print(names(fit$frame$data_frame))
}
f()
