# Extracted from test-v19.R:84

# prequel ----------------------------------------------------------------------
sim_lcm <- function(seed = 23, ng = 50, m = 8) {
  set.seed(seed)
  cls <- rbinom(ng, 1, 0.4)
  bg <- rnorm(ng, 0, 0.5)
  g <- rep(seq_len(ng), each = m)
  y <- rnorm(ng * m, (c(-1, 2)[cls + 1] + bg)[g], 0.8)
  list(dd = data.frame(y = y, g = factor(g)), cls = cls, m = m,
       ng = ng)
}
lcm_llk <- function(Ym, mu, tau, sig) {
  m <- nrow(Ym)
  s2 <- sig^2
  t2 <- tau^2
  ld <- m * log(s2) + log1p(m * t2 / s2)
  r <- Ym - mu
  q <- colSums(r^2) / s2 -
    (t2 / (s2 * (s2 + m * t2))) * colSums(r)^2
  -0.5 * (m * log(2 * pi) + ld + q)
}

# test -------------------------------------------------------------------------
dd <- data.frame(y = rnorm(60), x = rnorm(60),
                   g = factor(rep(1:6, 10)))
gp <- get_prior(bf(y ~ x + (x | g), sigma ~ x) + gaussian(),
                  data = dd)
expect_setequal(unique(gp$class), c("Intercept", "b", "sd", "theta"))
expect_true("x" %in% gp$coef[gp$class == "b" & gp$dpar == ""])
expect_true("x" %in% gp$coef[gp$class == "b" & gp$dpar == "sigma"])
expect_true("g" %in% gp$group[gp$class == "sd"])
expect_equal(sum(gp$class == "theta" & nzchar(gp$coef)), 4L)
