args <- commandArgs(trailingOnly=TRUE); n <- as.integer(args[1]); m <- args[2]
library(RTMB); library(RTMBode)
times <- c(0,1,2,3); set.seed(1)
obs <- exp(-0.3*rep(times[-1], each=n)) + rnorm(n*3, 0, 0.1)
f <- function(p) { getAll(p)
  sol <- RTMBode::ode(y=rep(1,n), times=times, func=function(t,y,q) list(-q*y),
                      parms=exp(mu+u), method=m)
  -sum(dnorm(u,0,1,log=TRUE)) + sum((obs - as.vector(t(sol[-1,-1,drop=FALSE])))^2) }
o <- MakeADFun(f, list(mu=log(0.3), u=numeric(n)), random="u", silent=TRUE)
g <- suppressWarnings(tryCatch(o$gr(log(0.3)), error=function(e) NaN))
cat(sprintf("n=%2d method=%-6s gr = %s\n", n, m, format(g, digits=6)))
