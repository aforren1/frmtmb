# Why mu_residuals() hands back the RESPONSE and not the centred
# response when mu is an intercept alone: sd(y - mean(y)) is not
# bitwise sd(y), and a start that moves in its last bit moves every
# iterate after it.
set.seed(7)
n <- 0L
for (i in 1:5000) {
  y <- rnorm(200, 50, 3)
  if (!identical(stats::sd(y), stats::sd(y - mean(y)))) n <- n + 1L
}
cat("sd(y) differs from sd(y - mean(y)) on", n, "of 5000 samples",
    "(n = 200, mean 50, sd 3)\n")
