# Which dense matrix operations RTMB tapes: solve(), determinant().
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
print(packageVersion("RTMB"))
print(grep("det|chol|solve|inv", ls(getNamespace("RTMB")), value = TRUE))
f <- function(p) {
  A <- RTMB::matrix(c(p[1], 0.3, 0.3, p[2]), 2, 2)
  v <- solve(A, c(p[1], 1))
  sum(v) + determinant(A, logarithm = TRUE)$modulus
}
print(f(c(2, 3)))
tp <- MakeTape(f, c(2, 3))
print(tp(c(2, 3)))
print(tp(c(2.5, 3.5)))
print(f(c(2.5, 3.5)))
print(tp$jacobian(c(2, 3)))
