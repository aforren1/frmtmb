# Which elementary functions tape in RTMB 2.0 for the three new densities.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
f <- function(p) {
  a <- p[1]; b <- p[2]
  c(log1p(a * 0.3), expm1(b), lgamma(1 - 1 / (2 + b^2)),
    exp(-log1p(a) / b), (2 + a^2)^b, -expm1(-a^2), log(-expm1(-a^2 - 0.1)),
    gamma(1.5 + a^2))
}
tp <- tryCatch(MakeTape(f, c(0.2, 0.5)), error = function(e) e)
print(tp)
if (!inherits(tp, "error")) {
  print(tp(c(0.2, 0.5)))
  print(f(c(0.2, 0.5)))
  print(tp$jacobian(c(0.2, 0.5)))
}
