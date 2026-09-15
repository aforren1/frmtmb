source("frailty-lib.R")
fp <- function(u, eps2 = 1e-14) {
  sg <- sign(u); s <- sqrt(u * u + eps2); w <- abs(u) + s
  0.25 * ((1 + sg) * w + (1 - sg) * (eps2 / w))
}
for (u in c(1, 0.1, 0.01, 1e-3, 1e-4, 0, -0.01, -0.2, -1)) {
  cat(sprintf("slope %8.4g  log(floor) %10.4f\n", u, log(fp(u))))
}
