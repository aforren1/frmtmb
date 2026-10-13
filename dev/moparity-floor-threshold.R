# Where sd_floored() and sd2_floored() equal exp() bit for bit.
f <- function(c, lo, hi) {
  # bisection on log x for the smallest x with x + c == x, checked on a
  # dense grid above
  for (i in 1:200) {
    mid <- (lo + hi) / 2
    if (exp(mid) + c == exp(mid)) hi <- mid else lo <- mid
  }
  grid <- seq(hi, 50, length.out = 2e6)
  cat(sprintf("c = %g: smallest identical log x %.6f (x %.4g); grid above: %s\n",
              c, hi, exp(hi), all(exp(grid) + c == exp(grid))))
}
f(1e-150, -400, -200)
f(1e-300, -800, -500)
th <- seq(-306.6, 50, length.out = 2e6)
cat("sd_floored identical on [-306.6, 50]:",
    identical(exp(th) + 1e-150, exp(th)), "\n")
