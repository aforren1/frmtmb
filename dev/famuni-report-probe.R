.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
inner <- function(p) { gevxi <- p * 2; RTMB::REPORT(gevxi); p^2 }
f <- function(pr) inner(pr$p) + sum(pr$b^2)
obj <- MakeADFun(f, list(p = 1.5, b = c(0.1, 0.2)), silent = TRUE)
print(obj$report(c(3, 0, 0)))
obj2 <- MakeADFun(f, list(p = 1.5, b = c(0.1, 0.2)), random = "b", silent = TRUE)
obj2$fn(4); print(obj2$report())
print(obj2$report(c(5, 0, 0)))
