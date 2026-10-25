.libPaths(c("C:/Users/adf44/source/r/rellib-r12", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
W <- Matrix::rsparsematrix(5, 5, 0.4)
Wd <- as.matrix(W)
f <- function(p) {
  r <- p - 0.3 * as.vector(W %*% p)
  sum(r^2)
}
g <- function(p) {
  r <- p - 0.3 * as.vector(Wd %*% p)
  sum(r^2)
}
x <- rnorm(5)
tp <- MakeTape(f, x); print(c(tp(x), f(x)))
tp2 <- MakeTape(g, x); print(c(tp2(x), g(x)))
print(all.equal(tp$jacobian(x), tp2$jacobian(x)))
