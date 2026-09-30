# Punch 1, B2: RTMB::pbeta()'s third derivatives where the candidate
# reads it: both shapes at least 150 (the clamp), x within 6.5 standard
# deviations of the mean. Seed 11, 20000 points, shapes log-uniform on
# [150, 1e7] each.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
F <- MakeTape(function(p) log(RTMB::pbeta(p[1], exp(p[2]), exp(p[3]))),
              c(0.2, 6, 6))
T3 <- F$jacfun()$jacfun()$jacfun()
set.seed(11)
M <- 20000
a <- exp(runif(M, log(150), log(1e7)))
b <- exp(runif(M, log(150), log(1e7)))
s <- a + b
m <- a / s
sd <- sqrt(a * b / (s * s * (s + 1)))
x <- m + runif(M, -6.5, 6.5) * sd
bad <- 0
for (i in seq_len(M)) {
  if (!all(is.finite(T3(c(x[i], log(a[i]), log(b[i])))))) bad <- bad + 1
}
cat("points", M, " non-finite third derivatives", bad, "\n")
