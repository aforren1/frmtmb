# Punch 1, B2: where is the 50-step continued fraction accurate, and where
# is RTMB::pbeta() safe to differentiate three times? A grid over the two
# shapes separately, x at and near m = (a + 1) / (a + b + 2) and far from
# it, x < 1/2 as xbeta needs. Error against stats::pbeta(log.p = TRUE),
# floored at one.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
cf <- function(x, a, b, N) frmtmb:::log_ibeta_half(x, a, b, N = N)
Tp <- MakeTape(function(p) log(RTMB::pbeta(p[1], exp(p[2]), exp(p[3]))),
               c(0.2, 0, 0))
T3 <- Tp$jacfun()$jacfun()$jacfun()
lg <- c(-2, -1, 0, 1, 2, 3, 4, 5, 6)
rows <- list()
for (la in lg) for (lb in lg) {
  a <- 10^la; b <- 10^lb
  m <- (a + 1) / (a + b + 2)
  sd <- sqrt(a * b / ((a + b)^2 * (a + b + 1)))
  xs <- c(m, m * (1 - 1e-3), m * (1 + 1e-3), m * 0.99, m * 1.01,
          m - sd, m + sd, m - 4 * sd, m + 4 * sd, m / 3, 0.45)
  xs <- xs[xs > 0 & xs < 0.5]
  for (x in xs) {
    ref <- pbeta(x, a, b, log.p = TRUE)
    if (!is.finite(ref) || ref < -600) next
    e50 <- abs(cf(x, a, b, 50L) - ref) / max(1, abs(ref))
    e200 <- abs(cf(x, a, b, 200L) - ref) / max(1, abs(ref))
    t3 <- all(is.finite(T3(c(x, log(a), log(b)))))
    rows[[length(rows) + 1L]] <- data.frame(la, lb, x, near = abs(x / m - 1) <= 0.011,
                                            e50, e200, t3)
  }
}
res <- do.call(rbind, rows)
cat("points", nrow(res), "\n")
cat("50 steps, worst error by (log10 a, log10 b):\n")
print(round(log10(pmax(tapply(res$e50, list(res$la, res$lb), max), 1e-17)), 1))
cat("200 steps:\n")
print(round(log10(pmax(tapply(res$e200, list(res$la, res$lb), max), 1e-17)), 1))
cat("RTMB::pbeta third derivative non-finite count by (log10 a, log10 b):\n")
print(tapply(!res$t3, list(res$la, res$lb), sum))
