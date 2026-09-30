# Punch 1, B2: at large shapes, how far from m (in standard deviations of
# the beta) must x be for the 50-step continued fraction to be exact?
# Random shapes with min(a, b) from 30 to 1e7, seed 1; error against
# stats::pbeta(log.p = TRUE), floored at one.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cf <- function(x, a, b, N) frmtmb:::log_ibeta_half(x, a, b, N = N)
set.seed(1)
M <- 6000
lmin <- runif(M, log(30), log(1e7))
lrat <- runif(M, 0, log(1e4))
a <- exp(lmin); b <- exp(lmin + lrat)
sw <- runif(M) < 0.5
tmp <- a; a[sw] <- b[sw]; b[sw] <- tmp[sw]
m <- (a + 1) / (a + b + 2)
sd <- sqrt(a * b / ((a + b)^2 * (a + b + 1)))
k <- runif(M, -12, 12)
x <- m + k * sd
ok <- x > 0 & x < 0.5
x <- x[ok]; a <- a[ok]; b <- b[ok]; k <- k[ok]
ref <- pbeta(x, a, b, log.p = TRUE)
keep <- is.finite(ref) & ref > -600
e <- sapply(which(keep), function(i) {
  abs(cf(x[i], a[i], b[i], 50L) - ref[i]) / max(1, abs(ref[i]))
})
kb <- cut(abs(k[keep]), c(0, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12))
sb <- cut(pmin(a, b)[keep], c(0, 100, 1e3, 1e4, 1e5, 1e6, 1e8))
cat("points", sum(keep), "\n")
cat("log10 worst error by |x - m| / sd (rows) and min(a, b) (cols):\n")
print(round(log10(pmax(tapply(e, list(kb, sb), max), 1e-17)), 1))
