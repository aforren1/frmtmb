# Agreement with drmTMB, a frequentist distributional-regression package
# built on a compiled TMB template. The two packages share a large model
# space, so each shared model is fitted in both and compared.
#
# Each comparison evaluates BOTH packages' objectives at ONE point, one
# standard error away from drmTMB's optimum, in addition to comparing the
# two optima. Agreement at the optima alone cannot reject a slightly wrong
# parameter map, because the gradient vanishes there; the off-optimum
# point can, and it is how drmTMB's 0.999999 cap on correlations was
# found. Every map between the two parameter vectors is written out.
#
# Tolerances are ratios: likelihood gaps to |logLik|, estimate gaps to
# the standard error, standard errors to each other. The gap between the
# two optima carries both optimizers' stopping error, so it gets 1e-9;
# the gap at one shared point off the optimum carries neither, so it
# gets 1e-11, which is what lets it see the correlation cap.
# dev/drmtmb-test-margins.R measures all five margins over the 24
# comparisons below: 19.1x for the gap between the optima, 26.4x at
# drmTMB's optimum, 72.1x off the optimum, 67.2x for estimates and
# 54.8x for standard errors.
#
# drmTMB fits are read through $obj, $opt and $sdr, which are not drmTMB
# API. A change there fails these tests loudly rather than skipping them.

# Three gates, on the pattern helper-brms.R's skip_unless_brms_fit()
# sets. Every test here FITS about 50 models in two packages, so it is a
# fit tier, not a structural one: without the environment variable an
# ordinary check job would install drmTMB and run all of them on every
# platform. drmTMB stays in Suggests because the file calls `drmTMB::`
# and the declaration is what keeps R CMD check quiet.
#   Sys.setenv(FRMTMB_DRMTMB_FIT_TESTS = "true", NOT_CRAN = "true")
skip_unless_drmtmb <- function() {
  skip_on_cran()
  skip_if_not_installed("drmTMB", "0.7.0")
  if (!identical(Sys.getenv("FRMTMB_DRMTMB_FIT_TESTS"), "true")) {
    skip("set FRMTMB_DRMTMB_FIT_TESTS=true to run drmTMB agreement tests")
  }
}

# drmTMB finds a structured effect's matrix or tree (A, tree) through the
# frame that calls drm_formula(), so a plain wrapper hides them. Both
# wrappers evaluate the real call in the test's own frame; drmTMB is
# only in Suggests, so it cannot be bound at load time.
drm_call_in_caller <- function(fn) {
  function(...) {
    cl <- sys.call()
    cl[[1]] <- fn
    eval(cl, parent.frame())
  }
}
drm_fit <- drm_call_in_caller(quote(drmTMB::drmTMB))
dbf <- drm_call_in_caller(quote(drmTMB::drm_formula))

# One row of a parameter map: frmtmb index f, drmTMB index d, the map
# from drmTMB's value to frmtmb's, and its derivative, which carries a
# standard error across the map.
drm_row <- function(f, d, fun, dfun) {
  list(f = f, d = d, fun = fun, dfun = dfun)
}
drm_lin <- function(f, d, a = 1, b = 0) {
  drm_row(f, d, function(x) a * x + b, function(x) a + 0 * x)
}
drm_ids <- function(k) lapply(seq_len(k), function(i) drm_lin(i, i))
# drmTMB: nu = 2 + exp(d). frmtmb, as brms: nu = 1 + exp(f).
drm_nu <- function(f, d) {
  drm_row(f, d, function(x) log1p(exp(x)),
          function(x) exp(x) / (1 + exp(x)))
}
# drmTMB: rho = 0.999999 tanh(d), for rho12 and for random-effect
# correlations. frmtmb: rho = f / sqrt(1 + f^2).
drm_rho <- function(f, d, cc = 0.999999) {
  drm_row(f, d,
          function(x) { r <- cc * tanh(x); r / sqrt(1 - r^2) },
          function(x) {
            r <- cc * tanh(x)
            cc * (1 - tanh(x)^2) / (1 - r^2)^1.5
          })
}

drm_measure <- function(fd, ff, map) {
  pd <- fd$opt$par
  pf <- ff$opt$par
  fi <- as.integer(vapply(map, `[[`, 0, "f"))
  di <- as.integer(vapply(map, `[[`, 0, "d"))
  stopifnot(setequal(fi, seq_along(pf)), setequal(di, seq_along(pd)))
  to_f <- function(p) {
    q <- pf
    for (m in map) q[m$f] <- m$fun(p[m$d])
    unname(q)
  }
  sed <- sqrt(diag(fd$sdr$cov.fixed))
  sef <- sqrt(diag(vcov(ff, full = TRUE)))
  # A fixed alternating sign pattern, not a draw: the off-optimum point
  # has to be the same on every run, and seeding it here would move the
  # global RNG stream that the surrounding tests simulate from.
  off <- pd + sed * rep(c(1, -1), length.out = length(pd))
  est_d <- vapply(map, function(m) m$fun(pd[m$d]), 0)
  se_d <- vapply(map, function(m) abs(m$dfun(pd[m$d])) * sed[m$d], 0)
  ll_f <- as.numeric(logLik(ff))
  list(
    ll_frm = ll_f,
    ll_drm = as.numeric(logLik(fd)),
    gap_at_drm = fd$obj$fn(unname(pd)) - ff$obj$fn(to_f(pd)),
    gap_off = fd$obj$fn(unname(off)) - ff$obj$fn(to_f(off)),
    diff_over_se = max(abs(pf[fi] - est_d) / sef[fi]),
    log_se_ratio = max(abs(log(sef[fi] / se_d)))
  )
}

expect_same_model <- function(fd, ff, map) {
  m <- drm_measure(fd, ff, map)
  scale <- abs(m$ll_frm)
  expect_lt(abs(m$ll_frm - m$ll_drm) / scale, 1e-9)
  expect_lt(abs(m$gap_at_drm) / scale, 1e-9)
  expect_lt(abs(m$gap_off) / scale, 1e-11)
  expect_lt(m$diff_over_se, 1e-2)
  expect_lt(m$log_se_ratio, 1e-3)
  invisible(m)
}

drm_sim <- function(seed = 101, ng = 40, nper = 10) {
  set.seed(seed)
  n <- ng * nper
  g <- factor(rep(seq_len(ng), each = nper))
  u <- rnorm(ng, 0, 0.6)[g]
  v <- rnorm(ng, 0, 0.3)[g]
  x <- rnorm(n)
  z <- rnorm(n)
  d <- data.frame(g = g, x = x, z = z, w = rnorm(ng)[g])
  d$y <- 1 + 0.5 * x + u + rnorm(n, 0, exp(-0.2 + 0.3 * z + v))
  m <- plogis(0.2 + 0.4 * x + u)
  d$yb <- rbeta(n, m * 20, (1 - m) * 20)
  d$yc <- rnbinom(n, mu = exp(1 + 0.3 * x + u), size = 3)
  d$yt <- 1 + 0.5 * x + u + 0.8 * rt(n, df = 5)
  lat <- 0.8 * x + u + rlogis(n)
  d$yo <- factor(cut(lat, c(-Inf, -1, 0.5, 2, Inf), labels = FALSE),
                 ordered = TRUE)
  d$ysn <- 1 + 0.5 * x + u + 0.8 * (abs(rnorm(n)) - sqrt(2 / pi))
  e1 <- d$y - 1 - 0.5 * x - u
  u2 <- 0.5 * u + rnorm(ng, 0, 0.4)[g]
  d$y2 <- 0.3 + 0.2 * x + u2 + 0.5 * e1 + rnorm(n, 0, 0.8)
  d$yln <- exp(0.2 + 0.3 * x + u + rnorm(n, 0, 0.5))
  d$ygam <- rgamma(n, shape = 4, scale = exp(0.3 + 0.3 * x + u) / 4)
  d$nt <- sample(5:15, n, TRUE)
  d$ys <- rbinom(n, d$nt, rbeta(n, m * 8, (1 - m) * 8))
  d$yzi <- ifelse(runif(n) < 0.25, 0L, rnbinom(n, mu = exp(1 + 0.3 * x),
                                                size = 2))
  d
}

test_that("gaussian location-scale with a mu random intercept, ML and REML", {
  skip_unless_drmtmb()
  d <- drm_sim()
  expect_same_model(
    drm_fit(dbf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d),
    frm(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d),
    drm_ids(5))
  expect_same_model(
    drm_fit(dbf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d,
            REML = TRUE),
    frm(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d,
        REML = TRUE),
    drm_ids(3))
})

test_that("random intercepts in mu and sigma, independent and correlated", {
  skip_unless_drmtmb()
  d <- drm_sim()
  expect_same_model(
    drm_fit(dbf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
            data = d),
    frm(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
        data = d),
    drm_ids(6))
  # drmTMB orders sd_mu, cor, sd_sigma; frmtmb sd_mu, sd_sigma, cor.
  expect_same_model(
    drm_fit(dbf(y ~ x + (1 | p | g), sigma ~ z + (1 | p | g)),
            family = gaussian(), data = d),
    frm(bf(y ~ x + (1 | p | g), sigma ~ z + (1 | p | g)), family = gaussian(),
        data = d),
    c(drm_ids(4), list(drm_lin(5, 5), drm_lin(6, 7), drm_rho(7, 6))))
})

test_that("beta and NB2: drmTMB's sigma is 1 / sqrt(phi or shape)", {
  skip_unless_drmtmb()
  d <- drm_sim()
  # log phi = -2 log sigma, and log shape = -2 log sigma.
  sig <- list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3, -2), drm_lin(4, 4))
  expect_same_model(
    drm_fit(dbf(yb ~ x + (1 | g), sigma ~ 1), family = drmTMB::beta(),
            data = d),
    frm(bf(yb ~ x + (1 | g), phi ~ 1), family = Beta(), data = d), sig)
  expect_same_model(
    drm_fit(dbf(yc ~ x + (1 | g), sigma ~ 1), family = drmTMB::nbinom2(),
            data = d),
    frm(bf(yc ~ x + (1 | g), shape ~ 1), family = negbinomial(), data = d),
    sig)
  # A random intercept in log sigma is one of SD s in log shape of SD 2s.
  expect_same_model(
    drm_fit(dbf(yc ~ x, sigma ~ 1 + (1 | g)), family = drmTMB::nbinom2(),
            data = d),
    frm(bf(yc ~ x, shape ~ 1 + (1 | g)), family = negbinomial(), data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3, -2),
         drm_lin(4, 4, 1, log(2))))
})

test_that("bivariate gaussian: rho12, correlated group effects, REML", {
  skip_unless_drmtmb()
  d <- drm_sim()
  expect_same_model(
    drm_fit(dbf(mu1 = y ~ x, mu2 = y2 ~ x, rho12 = ~ 1),
            family = drmTMB::biv_gaussian(), data = d),
    frm(mvbf(bf(y ~ x), bf(y2 ~ x)) + set_rescor(TRUE), family = gaussian(),
        data = d),
    c(drm_ids(6), list(drm_rho(7, 7))))
  re_map <- c(drm_ids(6), list(drm_rho(10, 7), drm_lin(7, 8), drm_lin(8, 9),
                               drm_rho(9, 10)))
  expect_same_model(
    drm_fit(dbf(mu1 = y ~ x + (1 | p | g), mu2 = y2 ~ x + (1 | p | g),
                rho12 = ~ 1), family = drmTMB::biv_gaussian(), data = d),
    frm(mvbf(bf(y ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
          set_rescor(TRUE), family = gaussian(), data = d),
    re_map)
  expect_same_model(
    drm_fit(dbf(mu1 = y ~ x + (1 | p | g), mu2 = y2 ~ x + (1 | p | g),
                rho12 = ~ 1), family = drmTMB::biv_gaussian(), data = d,
            REML = TRUE),
    frm(mvbf(bf(y ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
          set_rescor(TRUE), family = gaussian(), data = d, REML = TRUE),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_rho(6, 3), drm_lin(3, 4),
         drm_lin(4, 5), drm_rho(5, 6)))
})

test_that("animal model: drmTMB animal(A = A) is frmtmb gr(cov = A)", {
  skip_unless_drmtmb()
  set.seed(202)
  # Three generations of random mating from 20 founders; A by the tabular
  # method, so that neither package is compared against its own builder.
  ped <- data.frame(id = paste0("f", 1:20), sire = NA, dam = NA)
  parents <- ped$id
  for (k in 1:3) {
    ids <- paste0("g", k, "_", 1:40)
    half <- length(parents) %/% 2
    ped <- rbind(ped, data.frame(
      id = ids, sire = sample(parents[seq_len(half)], 40, TRUE),
      dam = sample(parents[-seq_len(half)], 40, TRUE)))
    parents <- ids
  }
  n <- nrow(ped)
  A <- matrix(0, n, n, dimnames = list(ped$id, ped$id))
  for (i in seq_len(n)) {
    s <- match(ped$sire[i], ped$id)
    m <- match(ped$dam[i], ped$id)
    for (j in seq_len(i - 1)) {
      A[i, j] <- A[j, i] <- ((if (is.na(s)) 0 else A[j, s]) +
                               (if (is.na(m)) 0 else A[j, m])) / 2
    }
    A[i, i] <- 1 + if (!is.na(s) && !is.na(m)) A[s, m] / 2 else 0
  }
  a <- as.vector(t(chol(A)) %*% rnorm(n)) * 0.7
  dd <- data.frame(id = factor(rep(ped$id, each = 3), levels = ped$id),
                   x = rnorm(3 * n))
  dd$y <- 1 + 0.4 * dd$x + a[as.integer(dd$id)] + rnorm(3 * n, 0, 0.8)
  expect_same_model(
    drm_fit(dbf(y ~ x + animal(1 | id, A = A), sigma ~ 1),
            family = gaussian(), data = dd),
    frm(bf(y ~ x + (1 | gr(id, cov = A))), family = gaussian(), data = dd,
        data2 = list(A = A)),
    drm_ids(4))
  expect_same_model(
    drm_fit(dbf(y ~ x + animal(1 | id, A = A), sigma ~ 1),
            family = gaussian(), data = dd, REML = TRUE),
    frm(bf(y ~ x + (1 | gr(id, cov = A))), family = gaussian(), data = dd,
        data2 = list(A = A), REML = TRUE),
    drm_ids(2))
})

test_that("phylo(tree) rescales the tree to unit height", {
  skip_unless_drmtmb()
  skip_if_not_installed("ape")
  set.seed(404)
  tree <- ape::rcoal(50)
  tree$tip.label <- paste0("sp", 1:50)
  Vt <- ape::vcv(tree)
  h <- Vt[1, 1]
  # rcoal's height is random; a unit-height tree would hide the rescale.
  expect_gt(abs(log(h)), 0.1)
  dp <- data.frame(sp = factor(rep(tree$tip.label, each = 4),
                               levels = tree$tip.label), x = rnorm(200))
  a <- as.vector(t(chol(Vt)) %*% rnorm(50)) * 0.8
  dp$y <- 0.5 + 0.3 * dp$x + a[as.integer(dp$sp)] + rnorm(200, 0, 0.6)
  fd <- drm_fit(dbf(y ~ x + phylo(1 | sp, tree = tree), sigma ~ 1),
                family = gaussian(), data = dp)
  # Against the raw covariance, the SD moves by sqrt(h).
  expect_same_model(fd,
    frm(bf(y ~ x + (1 | gr(sp, cov = Vt))), family = gaussian(), data = dp,
        data2 = list(Vt = Vt)),
    c(drm_ids(3), list(drm_lin(4, 4, 1, -log(h) / 2))))
  C <- ape::vcv(tree, corr = TRUE)
  expect_same_model(fd,
    frm(bf(y ~ x + (1 | gr(sp, cov = C))), family = gaussian(), data = dp,
        data2 = list(C = C)),
    drm_ids(4))
})

test_that("student and cumulative logit with a mu random intercept", {
  skip_unless_drmtmb()
  d <- drm_sim()
  expect_same_model(
    drm_fit(dbf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1),
            family = drmTMB::student(), data = d),
    frm(bf(yt ~ x + (1 | g), sigma ~ 1, nu ~ 1), family = student(),
        data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3), drm_nu(4, 4),
         drm_lin(5, 5)))
  # Both: first cutpoint raw, then log increments. frmtmb puts the SD
  # first.
  expect_same_model(
    drm_fit(dbf(yo ~ x + (1 | g)), family = drmTMB::cumulative_logit(),
            data = d),
    frm(bf(yo ~ x + (1 | g)), family = cumulative("logit"), data = d),
    list(drm_lin(1, 1), drm_lin(3, 2), drm_lin(4, 3), drm_lin(5, 4),
         drm_lin(2, 5)))
})

test_that("more families: lognormal, Gamma, beta-binomial, ZI NB2", {
  skip_unless_drmtmb()
  d <- drm_sim()
  expect_same_model(
    drm_fit(dbf(yln ~ x, sigma ~ 1 + (1 | g)), family = drmTMB::lognormal(),
            data = d),
    frm(bf(yln ~ x, sigma ~ 1 + (1 | g)), family = lognormal(), data = d),
    drm_ids(4))
  # Gamma sigma is the CV: shape = 1 / sigma^2.
  expect_same_model(
    drm_fit(dbf(ygam ~ x, sigma ~ 1 + (1 | g)), family = Gamma(link = "log"),
            data = d),
    frm(bf(ygam ~ x, shape ~ 1 + (1 | g)), family = Gamma(link = "log"),
        data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3, -2),
         drm_lin(4, 4, 1, log(2))))
  expect_same_model(
    drm_fit(dbf(cbind(ys, nt - ys) ~ x + (1 | g), sigma ~ 1),
            family = drmTMB::beta_binomial(), data = d),
    frm(bf(ys | trials(nt) ~ x + (1 | g), phi ~ 1), family = beta_binomial(),
        data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3, -2), drm_lin(4, 4)))
  expect_same_model(
    drm_fit(dbf(yzi ~ x, sigma ~ 1, zi ~ 1), family = drmTMB::nbinom2(),
            data = d),
    frm(bf(yzi ~ x, shape ~ 1, zi ~ 1), family = zero_inflated_negbinomial(),
        data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3, -2), drm_lin(4, 4)))
})

test_that("meta-analysis: meta_V(V = vi) is se(sei, sigma = TRUE)", {
  skip_unless_drmtmb()
  set.seed(303)
  vi <- runif(60, 0.02, 0.3)
  m <- data.frame(x = rnorm(60), vi = vi, sei = sqrt(vi))
  m$yi <- 0.3 + 0.2 * m$x + rnorm(60, 0, 0.25) + rnorm(60, 0, sqrt(vi))
  expect_same_model(
    drm_fit(dbf(yi ~ x + meta_V(V = vi), sigma ~ 1),
            family = gaussian(), data = m),
    frm(bf(yi | se(sei, sigma = TRUE) ~ x), family = gaussian(), data = m),
    drm_ids(3))
})

test_that("sd(g) ~ w is the nl model b0 + exp(c1 * w) * zz", {
  skip_unless_drmtmb()
  set.seed(505)
  ng <- 40
  g <- factor(rep(seq_len(ng), each = 10))
  wg <- rnorm(ng)
  d <- data.frame(g = g, x = rnorm(400), w = wg[g])
  d$y <- 1 + 0.5 * d$x + rnorm(ng, 0, exp(-0.5 + 0.5 * wg))[g] +
    rnorm(400, 0, 0.7)
  # log sd_g = c0 + c1 w_g in drmTMB; frmtmb has no such grammar, and the
  # SD of zz supplies c0 = log s.
  expect_same_model(
    drm_fit(dbf(y ~ x + (1 | g), sigma ~ 1, sd(g) ~ w), family = gaussian(),
            data = d),
    frm(bf(y ~ b0 + exp(lsd) * zz, b0 ~ x, lsd ~ 0 + w, zz ~ 0 + (1 | g),
           nl = TRUE), family = gaussian(), data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(4, 3), drm_lin(5, 4),
         drm_lin(3, 5)))
})

test_that("mi(): a gaussian missing predictor with its own model", {
  skip_unless_drmtmb()
  set.seed(506)
  n <- 400
  d <- data.frame(z = rnorm(n), w2 = rnorm(n))
  d$xm <- 0.6 * d$w2 + rnorm(n, 0, 0.8)
  d$y <- 1 + 0.3 * d$z + 0.4 * d$xm + rnorm(n, 0, 0.7)
  d$xm[sample(n, 60)] <- NA
  expect_same_model(
    drm_fit(dbf(y ~ z + mi(xm), sigma ~ 1), family = gaussian(), data = d,
            impute = list(xm = xm ~ w2),
            missing = drmTMB::miss_control(predictor = "model")),
    frm(bf(y ~ z + mi(xm)) + bf(xm | mi() ~ w2) + set_rescor(FALSE),
        family = gaussian(), data = d),
    list(drm_lin(1, 1), drm_lin(2, 2), drm_lin(3, 3), drm_lin(6, 4),
         drm_lin(4, 5), drm_lin(5, 6), drm_lin(7, 7)))
})

test_that("skew_normal: same likelihood; frmtmb's default start stalls", {
  skip_unless_drmtmb()
  set.seed(1)
  n <- 200
  xs <- -abs(rnorm(n)) * 3
  d <- data.frame(xs = xs, y = xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 +
                    rnorm(n, 0, 0.3))
  fd <- drm_fit(dbf(y ~ xs, sigma ~ 1, nu ~ 1),
                family = drmTMB::skew_normal(), data = d)
  # drmTMB's nu is brms's alpha, so the map is the identity. Started on
  # the right side of alpha = 0, frmtmb reaches drmTMB's optimum.
  expect_same_model(fd,
    frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = d,
        start = list(betad = c(0, 2))),
    drm_ids(4))
  # KNOWN DEFECT, pinned on purpose. frmtmb's alpha start takes its sign
  # from the skewness of the RAW response, which here is opposite to the
  # residual skewness, and alpha = 0 is a stationary point of the
  # mean-parameterized skew normal. The default fit stops there with
  # convergence code 0, about 25 log-likelihood units low. When the start
  # is fixed this expectation fails: flip it to expect_same_model(), do
  # not delete it.
  ff <- frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(),
            data = d)
  expect_identical(ff$opt$convergence, 0L)
  # Measured: 24.7 units on a logLik of -284, and alpha within 2e-4 of
  # its standard error from zero.
  expect_gt((as.numeric(logLik(fd)) - as.numeric(logLik(ff))) /
              abs(as.numeric(logLik(fd))), 0.01)
  se_alpha <- sqrt(diag(vcov(ff, full = TRUE)))[[4]]
  expect_lt(abs(ff$opt$par[[4]]) / se_alpha, 0.01)
})

test_that("REML with a sigma random effect: not the same criterion", {
  skip_unless_drmtmb()
  # drmTMB integrates beta_sigma out as well as beta_mu once sigma has a
  # random effect, and not otherwise; frmtmb integrates beta_mu only in
  # both. So there is no agreement to pin. What is pinned is frmtmb's
  # nesting: sigma ~ z is the sd -> 0 limit of sigma ~ z + (1 | g), and
  # its REML criterion must not fall when the zero-variance effect is
  # added.
  set.seed(1)
  g <- factor(rep(1:30, each = 8))
  d <- data.frame(g = g, x = rnorm(240), z = rnorm(240))
  d$y <- 1 + 0.5 * d$x + rnorm(30, 0, 0.6)[g] +
    rnorm(240, 0, exp(-0.2 + 0.3 * d$z))
  fA <- frm(bf(y ~ x + (1 | g), sigma ~ z), family = gaussian(), data = d,
            REML = TRUE)
  fB <- frm(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)), family = gaussian(),
            data = d, REML = TRUE)
  dA <- drm_fit(dbf(y ~ x + (1 | g), sigma ~ z), family = gaussian(),
                data = d, REML = TRUE)
  dB <- drm_fit(dbf(y ~ x + (1 | g), sigma ~ z + (1 | g)),
                family = gaussian(), data = d, REML = TRUE)
  scale <- abs(as.numeric(logLik(fA)))
  expect_lt(abs(as.numeric(logLik(fA)) - as.numeric(logLik(dA))) / scale,
            1e-9)
  expect_gt((as.numeric(logLik(fB)) - as.numeric(logLik(fA))) / scale,
            -1e-9)
  rand <- function(fit) unique(names(fit$obj$env$par)[fit$obj$env$random])
  expect_false("beta_sigma" %in% rand(dA))
  expect_true("beta_sigma" %in% rand(dB))
})
