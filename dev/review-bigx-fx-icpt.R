# What does a bound on class Intercept hold: the natural intercept or
# brms's centered one? x far from 0 against sd, no basis involved (m = 10).
# Then a non-binding intercept bound at 1e6: does the ill-design warning's
# fit stop short? REV_ARM selects the build.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
set.seed(11); n <- 200
z <- rnorm(n); d <- data.frame(z = z, x = z + 10)
d$y <- 1 + 0.4 * d$x + rnorm(n)        # natural intercept 1, centered 5
for (lb in c(3, -100)) {
  pr <- set_prior("normal(0, 100)", class = "Intercept", lb = lb)
  r <- rfit(bf(y ~ x), data = d, prior = pr)
  b <- quiet(fixef(r$f))
  cat(sprintf("x = z + 10, Intercept lb = %g: natural intercept %.4f slope %.4f ll %.4f warn: %s\n",
              lb, b[1, 1], b[2, 1], ll_of(r$f), short_w(r)))
}
cat("lm:", coef(lm(y ~ x, d)), "\n")
for (s in 1:3) {
  set.seed(s); n <- 150; G <- 15
  g <- factor(rep(seq_len(G), length.out = n)); z <- rnorm(n)
  eta <- 0.5 + 0.4 * z + rnorm(G, 0, 0.7)[g]
  d2 <- data.frame(y = rpois(n, exp(eta)), z = z, w = z + 1e6, g = g)
  pr <- set_prior("normal(0, 1e7)", class = "Intercept", lb = -1e7)
  fz <- quiet(frm(bf(y ~ z + (1 | g)), family = poisson(), data = d2, prior = pr))
  r <- rfit(bf(y ~ w + (1 | g)), family = poisson(), data = d2, prior = pr)
  bw <- quiet(fixef(r$f)); bz <- quiet(fixef(fz))
  cat(sprintf("1e6, Intercept lb = -1e7 (binds neither), seed %d: dll %+.3e shift %+.2f SE, natural icpt %.1f\n   warn: %s\n",
              s, ll_of(r$f) - ll_of(fz), (bw["w", 1] - bz["z", 1]) / bz["z", 2], bw[1, 1],
              short_w(r)))
}
