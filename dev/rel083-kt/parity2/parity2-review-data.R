# Review of lane parity2: data and formulas both the brms and frmtmb arms
# of dev/reviews/parity2-review-ce*.R read. Seed 21.
set.seed(21)
n <- 200
rd <- data.frame(x = runif(n, 0.1, 2), z = rnorm(n), w = rnorm(n),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 h = factor(sample(c("p", "q"), n, TRUE)),
                 k = sample(1:4, n, TRUE),
                 g = factor(sample(letters[1:8], n, TRUE)),
                 o = runif(n, 0, 0.5))
rd$xm <- rd$x + rnorm(n, sd = 0.2)
rd$sdx <- 0.2
rd$y <- 1 + rd$x * rd$z + 0.3 * rd$k + 0.5 * (rd$f == "b") + rnorm(n)
rd$y2 <- 0.5 * rd$x + rnorm(n)
rd$yc <- rpois(n, exp(0.3 + 0.4 * rd$x))
rd$yc[sample(n, 40)] <- 0L
shapes <- list(
  sby = list(f = bf(y ~ s(x, by = f) + f), fam = gaussian()),
  sbynum = list(f = bf(y ~ s(x, by = z)), fam = gaussian()),
  t2w = list(f = bf(y ~ t2(x, z) + w), fam = gaussian()),
  gp = list(f = bf(y ~ gp(x) + z), fam = gaussian()),
  mox = list(f = bf(y ~ mo(k) * x), fam = gaussian()),
  me = list(f = bf(y ~ me(xm, sdx) + z), fam = gaussian()),
  mez = list(f = bf(y ~ me(xm, sdx) * z), fam = gaussian()),
  mvbind = list(f = bf(mvbind(y, y2) ~ x + z) + set_rescor(FALSE),
                fam = gaussian()),
  bfbf = list(f = bf(y ~ x + w) + bf(y2 ~ z) + set_rescor(FALSE),
              fam = gaussian()),
  mix = list(f = bf(y ~ x, mu2 ~ z), fam = mixture(gaussian, gaussian)),
  dpint = list(f = bf(y ~ x, sigma ~ z * w), fam = gaussian(),
               dpar = "sigma"),
  nlfsig = list(f = bf(y ~ x) + nlf(sigma ~ a * z) + lf(a ~ 1 + w),
                fam = gaussian(), dpar = "sigma"),
  fac = list(f = bf(y ~ f * h), fam = gaussian()),
  off = list(f = bf(y ~ x + offset(o)), fam = gaussian()),
  fx = list(f = bf(y ~ f * x), fam = gaussian()),
  zip = list(f = bf(yc ~ x, zi ~ w), fam = zero_inflated_poisson(),
             dpar = "zi"),
  icpt = list(f = bf(y ~ 0 + Intercept + x, sigma ~ z), fam = gaussian()),
  mozw = list(f = bf(y ~ mo(k):z + w), fam = gaussian())
)
