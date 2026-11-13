# Lane names: the data and formulas both arms read (dev/names-brms.R for
# brms 2.23.0, dev/names-frm.R for frmtmb). Each case is a list of the
# formula, family and data, built with bf() and families that frmtmb and
# brms both export, so the one definition serves both.
names_cases <- function(bf, fam) {
  set.seed(20261111)
  n <- 120
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(sprintf("g%02d", 1:12), each = 10)),
                  h = factor(rep(1:8, length.out = n)),
                  ga = factor(rep(1:6, length.out = n)))
  ug <- rnorm(12, 0, 0.5)[d$g]
  d$y <- 1 + 0.5 * d$x + ug + rnorm(n)
  # categorical with three levels
  eta <- cbind(0, 0.3 + 0.6 * d$x + ug, -0.2 - 0.4 * d$x)
  p <- exp(eta) / rowSums(exp(eta))
  d$yc <- factor(apply(p, 1, function(pp) sample(c("a", "b", "c"), 1,
                                                 prob = pp)))
  d$ys <- d$y + rnorm(n, 0, 0.3)
  d$ydist <- 1 + 0.5 * d$x + ug + rnorm(n, 0, exp(0.2 * d$z))
  d$k <- sample(1:4, n, replace = TRUE)
  d$ymo <- 1 + 0.4 * d$k + rnorm(n)
  d$f <- factor(rep(c("fa", "fb"), length.out = n))
  d$ygp <- sin(d$x) + rnorm(n, 0, 0.3)
  d$sx <- 0.3
  d$sz <- runif(n, 0.2, 0.4)
  d$xg <- rnorm(12)[d$g]
  d$wg <- rnorm(12)[d$g]
  d$sg <- 0.25
  d$yme <- 1 + 0.5 * d$xg + 0.3 * d$z - 0.2 * d$wg + rnorm(n)
  d$t <- rep(1:10, 12)
  d$ycount <- rpois(n, exp(0.5 + 0.3 * d$x))
  d$y2 <- 0.5 * d$y + rnorm(n)
  d$yf <- 1 + 0.5 * d$x + rnorm(12)[d$g] * (1 + (d$f == "fb")) + rnorm(n)
  list(
    cat = list(f = bf(yc ~ x + (1 | g)), fam = fam$categorical, d = d),
    skew = list(f = bf(ys ~ (1 | h), alpha ~ (1 | g), sigma ~ (1 | ga)),
                fam = fam$skew_normal, d = d),
    dist = list(f = bf(ydist ~ x + (1 | g), sigma ~ z + (1 | g)),
                fam = fam$gaussian, d = d),
    distpop = list(f = bf(ydist ~ x, sigma ~ 1), fam = fam$gaussian,
                   d = d),
    twoblk = list(f = bf(y ~ x + (1 | g) + (0 + x | g)),
                  fam = fam$gaussian, d = d),
    smooth = list(f = bf(ygp ~ s(x)), fam = fam$gaussian, d = d),
    mo = list(f = bf(ymo ~ mo(k)), fam = fam$gaussian, d = d),
    gpby = list(f = bf(ygp ~ gp(x, by = f)), fam = fam$gaussian, d = d),
    me = list(f = bf(yme ~ me(xg, sg, gr = g) + me(z, sz) +
                       me(wg, sg, gr = g)),
              fam = fam$gaussian, d = d),
    latres = list(f = bf(ycount ~ x + ar(time = t, gr = g, p = 1)),
                  fam = fam$poisson, d = d),
    mv = list(f = bf(y ~ x + (1 | g)) + bf(y2 ~ z + (1 | g)),
              fam = fam$gaussian, d = d),
    dblbar = list(f = bf(yf ~ x + (1 + f || g)), fam = fam$gaussian, d = d)
  )
}
