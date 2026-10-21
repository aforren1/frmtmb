# Are RTMB::ppois()'s third derivatives in the mean finite? Seed 3.
.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
set.seed(3)
bad <- 0; tot <- 0
for (k in c(0, 1, 3, 10, 40)) {
  f <- function(p) log(RTMB::ppois(k, exp(p[1])))
  T3 <- MakeTape(f, 0)$jacfun()$jacfun()
  for (lp in runif(100, -3, 5)) {
    tot <- tot + 1
    if (!all(is.finite(T3$jacobian(lp)))) bad <- bad + 1
  }
}
cat("log(RTMB::ppois) third derivative non-finite:", bad, "of", tot, "\n")
