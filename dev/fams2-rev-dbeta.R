# Reviewer: when is RTMB::dbeta()'s gradient NaN at large shapes?
.libPaths(c("C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
set.seed(1); y <- rbeta(50, 1500, 3500)
f1 <- MakeTape(function(p) sum(RTMB::dbeta(y, p[1], p[2], log = TRUE)), c(1, 1))
f2 <- MakeTape(function(p) sum(RTMB::dbeta(0.3, p[1], p[2], log = TRUE)), c(1, 1))
f3 <- MakeTape(function(p) sum(RTMB::dbeta(y + 0 * p[1], p[1], p[2], log = TRUE)), c(1, 1))
f4 <- MakeTape(function(p) sum(RTMB::dbeta(rep(0.3, 2), p[1], p[2], log = TRUE)), c(1, 1))
for (s in c(1000, 1500, 5000, 2e4)) {
  p <- c(0.3 * s, 0.7 * s)
  cat(s, "vector y const:", f1$jacobian(p), " scalar 0.3 const:", f2$jacobian(p),
      " vector y taped:", f3$jacobian(p), " rep(0.3, 2):", f4$jacobian(p), "\n")
}
