# v0.24 method-surface residue: deviance residuals across the family set
# and joint delta-method standard errors for the expected response.

# glm's IRLS stops on a relative deviance change, which leaves the
# coefficients good to ~1e-4 only; both sides need pinning before an
# exactness claim means anything.
tight_glm <- function() stats::glm.control(epsilon = 1e-14, maxit = 200)

res_env <- new.env()
res_env$dd <- local({
  set.seed(11)
  n <- 120
  x <- rnorm(n)
  data.frame(
    x = x,
    gau = 1 + 0.5 * x + rnorm(n),
    cnt = rpois(n, exp(0.2 + 0.4 * x)),
    ber = rbinom(n, 1, plogis(-0.2 + 0.9 * x)),
    k = rbinom(n, 5, plogis(-0.2 + 0.9 * x)),
    m = 5,
    gam = rgamma(n, 3, rate = 3 / exp(0.5 + 0.3 * x)),
    ig = RTMBdist::rinvgauss(n, mean = exp(0.5 + 0.2 * x), shape = 5),
    w = runif(n, 0.5, 2)
  )
})

test_that("deviance residuals reproduce stats::glm exactly", {
  dd <- res_env$dd
  # start at the glm solution so both optima are the same point; the
  # unit deviance is what is under test, not the optimizer
  cmp <- function(fit, ref) {
    expect_vector_equal(residuals(fit, type = "deviance"),
                        residuals(ref, type = "deviance"), tol = 1e-8)
  }

  g <- glm(gau ~ x, data = dd, control = tight_glm())
  f <- frm(bf(gau ~ x) + gaussian(), data = dd,
           start = list(beta = unname(coef(g)),
                        betad = log(sqrt(mean(residuals(g)^2)))))
  cmp(f, g)

  g <- glm(cnt ~ x, data = dd, family = poisson(), control = tight_glm())
  cmp(frm(bf(cnt ~ x) + poisson(), data = dd,
          start = list(beta = unname(coef(g)))), g)

  g <- glm(ber ~ x, data = dd, family = binomial(), control = tight_glm())
  cmp(frm(bf(ber ~ x) + bernoulli(), data = dd,
          start = list(beta = unname(coef(g)))), g)

  # counts out of trials(): the same unit deviance glm() computes from
  # the proportion and the prior weight
  g <- glm(cbind(k, m - k) ~ x, data = dd, family = binomial(),
           control = tight_glm())
  cmp(frm(bf(k | trials(m) ~ x) + binomial(), data = dd,
          start = list(beta = unname(coef(g)))), g)

  # a fixed dispersion isolates the mean: the Gamma and inverse gaussian
  # unit deviances do not involve the shape at all
  g <- glm(gam ~ x, data = dd, family = Gamma(link = "log"),
           control = tight_glm())
  cmp(frm(bf(gam ~ x, shape = 2) + Gamma(link = "log"), data = dd,
          start = list(beta = unname(coef(g)))), g)

  g <- glm(ig ~ x, data = dd, family = inverse.gaussian(link = "log"),
           control = tight_glm())
  cmp(frm(bf(ig ~ x, shape = 5) + stats::inverse.gaussian(link = "log"),
          data = dd, start = list(beta = unname(coef(g)))), g)
})

test_that("weights() multiply the unit deviance the glm way", {
  dd <- res_env$dd
  g <- glm(cnt ~ x, data = dd, family = poisson(), weights = dd$w,
           control = tight_glm())
  f <- frm(bf(cnt | weights(w) ~ x) + poisson(), data = dd,
           start = list(beta = unname(coef(g))))
  expect_vector_equal(residuals(f, type = "deviance"),
                      residuals(g, type = "deviance"), tol = 1e-8)
  # the weight enters the unit deviance, so the residual scales by
  # sqrt(w), not by w
  mu <- fitted(f)
  d <- 2 * (ifelse(dd$cnt > 0, dd$cnt * log(dd$cnt / mu), 0) -
              (dd$cnt - mu))
  expect_vector_equal(residuals(f, type = "deviance"),
                      sign(dd$cnt - mu) * sqrt(dd$w * d), tol = 1e-12)
})

test_that("negative-binomial deviance residuals match glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(11)
  n <- 120
  dn <- data.frame(x = rnorm(n))
  dn$y <- rnbinom(n, size = 2, mu = exp(1 + 0.4 * dn$x))

  f2 <- frm(bf(y ~ x) + negbinomial(), data = dn)
  g2 <- glmmTMB::glmmTMB(y ~ x, data = dn, family = glmmTMB::nbinom2)
  expect_vector_equal(residuals(f2, type = "deviance"),
                      as.numeric(residuals(g2, type = "deviance")),
                      tol = 1e-5)
  # at a shared mu and shape the two definitions agree to machine noise
  mu <- fitted(f2)
  sh <- unique(round(predict(f2, dpar = "shape", type = "response"), 12))
  dv <- stats::family(g2)$dev.resids(dn$y, mu, 1, theta = sh)
  expect_vector_equal(residuals(f2, type = "deviance"),
                      sign(dn$y - mu) * sqrt(pmax(dv, 0)), tol = 1e-10)

  f1 <- frm(bf(y ~ x) + nbinom1(), data = dn)
  g1 <- glmmTMB::glmmTMB(y ~ x, data = dn, family = glmmTMB::nbinom1)
  expect_vector_equal(residuals(f1, type = "deviance"),
                      as.numeric(residuals(g1, type = "deviance")),
                      tol = 1e-5)
  mu1 <- fitted(f1)
  ph1 <- unique(round(predict(f1, dpar = "phi", type = "response"), 12))
  dv1 <- stats::family(g1)$dev.resids(dn$y, mu1, 1, phi = ph1)
  expect_vector_equal(residuals(f1, type = "deviance"),
                      sign(dn$y - mu1) * sqrt(pmax(dv1, 0)), tol = 1e-10)

  # geometric is negbinomial with the shape pinned at 1
  fg <- frm(bf(y ~ x) + geometric(), data = dn)
  mug <- fitted(fg)
  dvg <- 2 * (ifelse(dn$y > 0, dn$y * log(dn$y / mug), 0) -
                (dn$y + 1) * log((dn$y + 1) / (mug + 1)))
  expect_vector_equal(residuals(fg, type = "deviance"),
                      sign(dn$y - mug) * sqrt(pmax(dvg, 0)), tol = 1e-12)
})

test_that("beta and tweedie deviances follow the saturated likelihood", {
  set.seed(11)
  n <- 120
  db <- data.frame(x = rnorm(n))
  db$y <- rbeta(n, 2, 3)
  fb <- frm(bf(y ~ x) + Beta(), data = db)
  mu <- fitted(fb)
  ph <- unique(round(predict(fb, dpar = "phi", type = "response"), 12))
  ll <- function(y, m) stats::dbeta(y, m * ph, (1 - m) * ph, log = TRUE)
  dev <- 2 * (ll(db$y, db$y) - ll(db$y, mu))
  expect_vector_equal(residuals(fb, type = "deviance"),
                      sign(db$y - mu) * sqrt(pmax(dev, 0)), tol = 1e-10)

  set.seed(12)
  dt <- data.frame(x = rnorm(n))
  dt$y <- ifelse(runif(n) < 0.3, 0,
                 rgamma(n, 2, rate = 2 / exp(0.3 + 0.2 * dt$x)))
  ft <- frm(bf(y ~ x) + tweedie(), data = dt)
  mu <- fitted(ft)
  pht <- unique(round(predict(ft, dpar = "phi", type = "response"), 12))
  pwr <- unique(round(predict(ft, dpar = "power", type = "response"), 12))
  # d = 2 phi (ll_sat - ll_fit); dtweedie's series is good to ~1e-7
  llt <- function(y, m) RTMB::dtweedie(y, m, pht, pwr, log = TRUE)
  dev <- 2 * pht * (llt(dt$y, ifelse(dt$y > 0, dt$y, 1e-8)) -
                      llt(dt$y, mu))
  rt <- residuals(ft, type = "deviance")
  expect_vector_equal(rt, sign(dt$y - mu) * sqrt(pmax(dev, 0)), tol = 1e-4)
  # the zero rows have a closed form: d = 2 mu^(2-p) / (2-p)
  i0 <- dt$y == 0
  expect_true(any(i0))
  expect_vector_equal(rt[i0],
                      -sqrt(2 * mu[i0]^(2 - pwr) / (2 - pwr)), tol = 1e-10)
})

test_that("deviance residuals condition on the random-effect modes", {
  dd <- sim_pois_glmm(seed = 77, n_g = 12, n_per = 12)
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  mu <- fitted(fit)   # conditional on the modes, the glmmTMB convention
  dev <- 2 * (ifelse(dd$y > 0, dd$y * log(dd$y / mu), 0) - (dd$y - mu))
  expect_vector_equal(residuals(fit, type = "deviance"),
                      sign(dd$y - mu) * sqrt(pmax(dev, 0)), tol = 1e-12)
  # deviance() is a different quantity and stays -2 logLik (lme4)
  expect_equal(deviance(fit), -2 * as.numeric(logLik(fit)))
  expect_gt(abs(deviance(fit) -
                  sum(residuals(fit, type = "deviance")^2)), 1)
})

test_that("families without a unit deviance are refused by name", {
  dd <- res_env$dd
  fs <- frm(bf(gau ~ x) + student(), data = dd)
  expect_error(residuals(fs, type = "deviance"),
               "no standard unit deviance")
  # the refusal names the families that do have one
  expect_error(residuals(fs, type = "deviance"), "inverse.gaussian")

  set.seed(4)
  dz <- data.frame(x = rnorm(120))
  dz$y <- rbinom(120, 1, 0.7) * rpois(120, 2)
  fz <- suppressWarnings(frm(bf(y ~ x) + zero_inflated_poisson(),
                             data = dz))
  expect_error(residuals(fz, type = "deviance"), "zero_inflated_poisson")

  do <- data.frame(x = rnorm(120), y = sample(1:3, 120, TRUE))
  fo <- frm(bf(y ~ x) + cumulative(), data = do)
  expect_error(residuals(fo, type = "deviance"), "cumulative")

  fm <- frm(bf(gau ~ x) + mixture(gaussian(), gaussian()), data = dd)
  expect_error(residuals(fm, type = "deviance"), "no standard unit")
})

test_that("trunc() and cens() responses refuse deviance residuals", {
  set.seed(5)
  n <- 200
  dt <- data.frame(x = rnorm(n))
  dt$y <- rnorm(n, 1 + 0.5 * dt$x, 1)
  dtr <- dt[dt$y > 0, ]
  ftr <- frm(bf(y | trunc(lb = 0) ~ x) + gaussian(), data = dtr)
  expect_error(residuals(ftr, type = "deviance"), "trunc")

  dc <- dt
  dc$cen <- as.numeric(dc$y > 2)
  dc$y <- pmin(dc$y, 2)
  fc <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dc)
  expect_error(residuals(fc, type = "deviance"), "cens")
})

# --- se.fit for the expected response ---------------------------------

se_env <- new.env()
se_env$dd <- local({
  set.seed(3)
  n <- 300
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(1:20, each = 15)))
  re <- rnorm(20, 0, 0.5)
  d$y <- rbinom(n, 1, 1 - plogis(-0.6 + 0.7 * d$z)) *
    rpois(n, exp(0.4 + 0.5 * d$x + re[d$g]))
  d
})

test_that("expected-response se.fit matches glmmTMB (zi and plain)", {
  skip_if_not_installed("glmmTMB")
  dd <- se_env$dd
  # plain poisson: the mean IS mu, so this exercises the unchanged
  # single-predictor path against the same reference
  fp <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  gp <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = dd, family = poisson)
  pp <- predict(fp, type = "response", se.fit = TRUE)
  gpp <- predict(gp, type = "response", se.fit = TRUE)
  expect_lt(max(abs(pp$fit - gpp$fit) / gpp$fit), 1e-3)
  expect_lt(max(abs(pp$se.fit - gpp$se.fit) / gpp$se.fit), 1e-3)

  fz <- suppressWarnings(
    frm(bf(y ~ x + (1 | g), zi ~ z) + zero_inflated_poisson(), data = dd))
  gz <- glmmTMB::glmmTMB(y ~ x + (1 | g), ziformula = ~z, data = dd,
                         family = poisson)
  pz <- predict(fz, type = "response", se.fit = TRUE)
  gzz <- predict(gz, type = "response", se.fit = TRUE)
  expect_lt(max(abs(pz$fit - gzz$fit) / gzz$fit), 1e-3)
  expect_lt(max(abs(pz$se.fit - gzz$se.fit) / gzz$se.fit), 1e-3)

  # newdata reproduces the in-sample answer exactly
  pn <- predict(fz, newdata = dd, type = "response", se.fit = TRUE)
  expect_equal(unname(pn$fit), unname(pz$fit))
  expect_equal(unname(pn$se.fit), unname(pz$se.fit))
})

test_that("the FD gradients reproduce the analytic ones (zi, lognormal)", {
  dd <- se_env$dd
  fz <- suppressWarnings(
    frm(bf(y ~ x, zi ~ z) + zero_inflated_poisson(), data = dd))
  p <- predict(fz, type = "response", se.fit = TRUE)
  mu <- predict(fz, type = "conditional")
  zi <- predict(fz, type = "zprob")
  # m = (1 - p) mu, log link on mu and logit on zi:
  # dm/deta_mu = (1 - p) mu, dm/deta_zi = -mu p (1 - p)
  G <- cbind(((1 - zi) * mu) * stats::model.matrix(~x, dd),
             (-mu * zi * (1 - zi)) * stats::model.matrix(~z, dd))
  V <- as.matrix(vcov(fz, full = TRUE))
  expect_vector_equal(p$se.fit, sqrt(rowSums((G %*% V) * G)), tol = 1e-8)

  set.seed(9)
  dl <- data.frame(x = rnorm(200), z = rnorm(200))
  dl$y <- rlnorm(200, 0.3 + 0.4 * dl$x, exp(-0.5 + 0.3 * dl$z))
  fl <- frm(bf(y ~ x, sigma ~ z) + frmtmb::lognormal(), data = dl)
  pl <- predict(fl, type = "response", se.fit = TRUE)
  mul <- predict(fl, dpar = "mu", type = "response")
  sgl <- predict(fl, dpar = "sigma", type = "response")
  m <- exp(mul + sgl^2 / 2)
  # identity link on mu, log link on sigma: dm/deta_mu = m,
  # dm/deta_sigma = m sigma^2
  Gl <- cbind(m * stats::model.matrix(~x, dl),
              (m * sgl^2) * stats::model.matrix(~z, dl))
  Vl <- as.matrix(vcov(fl, full = TRUE))
  expect_equal(unname(pl$fit), unname(m))
  expect_vector_equal(pl$se.fit, sqrt(rowSums((Gl %*% Vl) * Gl)),
                      tol = 1e-8)
})

test_that("se.fit covers trials-binomial and truncated responses", {
  set.seed(6)
  db <- data.frame(x = rnorm(150), m = 8)
  db$k <- rbinom(150, 8, plogis(0.2 + 0.6 * db$x))
  fb <- frm(bf(k | trials(m) ~ x) + binomial(), data = db)
  pb <- predict(fb, type = "response", se.fit = TRUE)
  pc <- predict(fb, type = "conditional", se.fit = TRUE)
  # the mean is trials * p, so both the value and its SE scale by trials
  expect_vector_equal(pb$fit, 8 * pc$fit, tol = 1e-12)
  expect_vector_equal(pb$se.fit, 8 * pc$se.fit, tol = 1e-10)

  set.seed(7)
  dt <- data.frame(x = rnorm(300))
  dt$y <- rnorm(300, 1 + 0.5 * dt$x, 1)
  dt <- dt[dt$y > 0, ]
  ft <- frm(bf(y | trunc(lb = 0) ~ x) + gaussian(), data = dt)
  pt <- predict(ft, type = "response", se.fit = TRUE)
  expect_true(all(is.finite(pt$se.fit)) && all(pt$se.fit > 0))
  # the point predictions are the ones the no-se path gives
  expect_equal(unname(pt$fit),
               unname(predict(ft, type = "response")))
  # and the truncated mean is above the untruncated one here
  expect_true(all(pt$fit > predict(ft, type = "conditional")))
})

test_that("the identity-mean se.fit path is untouched", {
  dd <- se_env$dd
  fp <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  pl <- predict(fp, type = "link", se.fit = TRUE)
  pr <- predict(fp, type = "response", se.fit = TRUE)
  # plain one-predictor delta method, not the joint one
  expect_equal(pr$se.fit, exp(pl$fit) * pl$se.fit)
})

test_that("se.fit still refuses a nonlinear predictor's response mean", {
  set.seed(8)
  dn <- data.frame(x = runif(80, 0, 5))
  dn$y <- 3 * exp(-0.8 * dn$x) + rnorm(80, 0, 0.1)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) +
               gaussian(), data = dn)
  expect_error(predict(fit, se.fit = TRUE), "se.fit is not supported")
})
