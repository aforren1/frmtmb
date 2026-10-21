.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
fd <- function(f, x0, h = 1e-6) vapply(seq_along(x0), function(j) { e <- replace(numeric(length(x0)), j, h); (f(x0 + e) - f(x0 - e)) / (2 * h) }, 0)
f1 <- function(p) RTMB::pbeta(p[1], p[2], p[3])
for (x0 in list(c(0.3, 2, 5), c(0.3, 2.0001, 5), c(0.6, 3, 1.5), c(0.3, 0.5, 0.5), c(0.1, 20, 30))) {
  F <- MakeTape(f1, x0)
  cat(sprintf("x=%g a=%g b=%g  tape %s  fd %s\n", x0[1], x0[2], x0[3],
    paste(signif(F$jacobian(x0), 8), collapse=","), paste(signif(fd(function(p) stats::pbeta(p[1],p[2],p[3]), x0), 8), collapse=",")))
}
# tape at one point, evaluate at another
F <- MakeTape(f1, c(0.5, 1, 1))
x0 <- c(0.3, 2, 5); cat("retaped-at-(.5,1,1) jac at (0.3,2,5):", F$jacobian(x0), "\n")
f2 <- function(p) RTMB::pbeta(0.3, p[1], 5)
F <- MakeTape(f2, 2); cat("pbeta shape1 only:", F$jacobian(2), " fd", fd(function(p) stats::pbeta(0.3, p, 5), 2), "\n")
f3 <- function(p) RTMB::pbeta(0.3, p[1], p[2])
F <- MakeTape(f3, c(2, 5)); cat("pbeta both shapes:", F$jacobian(c(2, 5)), "\n")
f4 <- function(p) RTMB::pbeta(0.3, exp(p[1]), exp(p[2]))
F <- MakeTape(f4, log(c(2, 5))); cat("pbeta exp shapes:", F$jacobian(log(c(2, 5))), " fd", fd(function(p) stats::pbeta(0.3, exp(p[1]), exp(p[2])), log(c(2,5))), "\n")
print(packageVersion("TMB"))
