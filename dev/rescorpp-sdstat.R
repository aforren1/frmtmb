## Lane rescorpp: is the sd-ratio z of dev/rescorpp-compare.R centered?
## A null simulation with no frmtmb code: e = sigma_s * z, 2000 draws x
## 200 rows, sigma_s ~ N(1, 0.05). Seed 7, 300 replicates.
set.seed(7)
S <- 2000; n <- 200
sg <- rnorm(S, 1, 0.05)
tg <- sqrt(mean(sg^2))
z <- replicate(300, {
  e <- sg * matrix(rnorm(S * n), S, n)
  d <- apply(e, 2, sd) - tg
  mean(d) / (sd(d) / sqrt(n))
})
cat("null z: mean", mean(z), " sd", sd(z), " share < 0", mean(z < 0), "\n")
