.libPaths(c("C:/Users/adf44/source/r/rellib-r14", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
f <- function(p) { z <- p$z; e <- exp(c(0, z)); w <- e / sum(e); (p$b - 2)^2 + sum(w) }
obj <- MakeADFun(f, list(b = 0, z = numeric(0)), silent = TRUE)
print(obj$par); print(obj$fn(obj$par)); print(obj$gr(obj$par)); print(obj$he(obj$par))
o <- nlminb(obj$par, obj$fn, obj$gr); print(o$par)
print(sdreport(obj))
print(obj$env$parList())
