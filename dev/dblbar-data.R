# Shared data for lane dblbar's probes. Seed 20261013.
dblbar_data <- function(n = 360, seed = 20261013) {
  set.seed(seed)
  g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
  f <- factor(sample(c("a", "b", "c"), n, TRUE))
  h <- factor(sample(c("u", "v"), n, TRUE))
  g1 <- factor(sample(sprintf("p%d", 1:6), n, TRUE))
  g2 <- factor(sample(sprintf("q%d", 1:5), n, TRUE))
  m1 <- factor(sample(sprintf("m%02d", 1:12), n, TRUE))
  m2 <- factor(sample(sprintf("m%02d", 1:12), n, TRUE))
  x <- round(rnorm(n), 3)
  z <- round(runif(n, -1, 1), 3)
  xo <- sample(0:3, n, TRUE)
  ug <- rnorm(24, 0, 0.8)[as.integer(g)]
  uf <- matrix(rnorm(24 * 3, 0, 0.6), 24)[cbind(as.integer(g),
                                                  as.integer(f))]
  ux <- rnorm(24, 0, 0.5)[as.integer(g)] * x
  y <- round(1 + 0.5 * x - 0.3 * z + c(0, 0.4, -0.2)[as.integer(f)] +
               ug + uf + ux + rnorm(n, 0, 0.7), 3)
  yo <- cut(y, c(-Inf, quantile(y, c(0.3, 0.6, 0.85)), Inf),
            labels = FALSE)
  hg <- factor(ifelse(as.integer(g) <= 12L, "u", "v"))
  data.frame(y, yo, x, z, xo, f, h, hg, g, g1, g2, m1, m2)
}
