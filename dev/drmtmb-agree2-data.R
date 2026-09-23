# Data for the second batch of shared models: the remaining families
# both packages implement, the random-effect SD formula, and mi().
sim_batch2 <- function(seed, ng = 40, nper = 10) {
  set.seed(seed)
  n <- ng * nper
  g <- factor(rep(seq_len(ng), each = nper))
  u <- rnorm(ng, 0, 0.5)[g]
  wg <- rnorm(ng)
  x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(g = g, x = x, z = z, w = wg[g])
  d$nt <- sample(5:15, n, TRUE)
  p <- rbeta(n, plogis(0.3 * x + u) * 8, (1 - plogis(0.3 * x + u)) * 8)
  d$ys <- rbinom(n, d$nt, p)
  # Tweedie by its compound Poisson-gamma construction, power 1.5.
  mu <- exp(0.5 + 0.3 * x + u); phi <- 0.8; pw <- 1.5
  lam <- mu^(2 - pw) / (phi * (2 - pw))
  shp <- (2 - pw) / (pw - 1); scl <- phi * (pw - 1) * mu^(pw - 1)
  N <- rpois(n, lam)
  d$ytw <- vapply(seq_len(n),
                  function(i) sum(rgamma(N[i], shp, scale = scl[i])), 0)
  d$yln <- exp(0.2 + 0.3 * x + u + rnorm(n, 0, 0.5))
  d$ygam <- rgamma(n, shape = 4, scale = exp(0.3 + 0.3 * x + u) / 4)
  nb <- rnbinom(n, mu = exp(1 + 0.3 * x), size = 2)
  d$yzi <- ifelse(runif(n) < 0.25, 0L, nb)
  d$ypos <- NA_integer_
  for (i in seq_len(n)) {
    repeat { k <- rnbinom(1, mu = exp(0.8 + 0.3 * x[i] + u[i]), size = 2)
             if (k > 0) break }
    d$ypos[i] <- k
  }
  m <- plogis(0.2 + 0.4 * x)
  yb <- rbeta(n, m * 10, (1 - m) * 10)
  r <- runif(n)
  d$yzob <- ifelse(r < 0.1, 0, ifelse(r < 0.2, 1, yb))
  # Group SD depends on the group-level covariate w: log sd = -0.5 + 0.5 w.
  ug <- rnorm(ng, 0, exp(-0.5 + 0.5 * wg))[g]
  d$y <- 1 + 0.5 * x + ug + rnorm(n, 0, 0.7)
  d$w2 <- rnorm(n)
  d$xm <- 0.6 * d$w2 + rnorm(n, 0, 0.8)
  d$y <- d$y + 0.4 * d$xm
  d
}
