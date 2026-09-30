# Punch 2, n1: where RTMB::pbeta()'s derivatives are NaN near the mean.
# Tried at x = a / (a + b), at the mode (a - 1) / (a + b - 2), and at
# (a + 1) / (a + b + 1), the mean of the (a + 1, b) beta, each exactly
# and one relative 1e-15 and 1e-12 away. Taped in (x, a, b).
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
F <- MakeTape(function(p) log(RTMB::pbeta(p[1], p[2], p[3])), c(0.3, 200, 400))
J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
fin <- function(p) c(g = all(is.finite(J(p))), h = all(is.finite(H(p))),
                     t = all(is.finite(T3(p))))
for (ab in list(c(300, 700), c(150, 150), c(151, 1e7), c(2e4, 8e4),
                c(1e7, 1e7), c(268.941421, 731.058579))) {
  a <- ab[1]; b <- ab[2]; s <- a + b
  pts <- c(mean = a / s, mode = (a - 1) / (s - 2), mean_a1 = (a + 1) / (s + 1))
  for (nm in names(pts)) for (dl in c(0, 1e-15, 1e-12)) {
    f <- fin(c(pts[[nm]] * (1 + dl), a, b))
    cat(sprintf("a %-10g b %-10g %-8s rel %-6g grad %s hess %s third %s\n",
                a, b, nm, dl, f[["g"]], f[["h"]], f[["t"]]))
  }
}
