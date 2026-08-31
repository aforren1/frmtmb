# v0.15: tier-2 sweep - covstruct variants, sequential ordinal
# families, mo(), influence(), CE prediction intervals and condition
# sets, vint()/vreal(), frm_multiple(), and the draws method surface.

test_that("hetar1/homcs/homtoep match glmmTMB (or a self-consistency)", {
  skip_if_not_installed("glmmTMB")
  dd <- local({
    set.seed(51)
    n_g <- 40
    n_t <- 5
    g <- factor(rep(seq_len(n_g), each = n_t))
    tim <- factor(rep(seq_len(n_t), n_g))
    u <- replicate(n_g, {
      e <- rnorm(n_t)
      for (t in 2:n_t) e[t] <- 0.5 * e[t - 1] + rnorm(1, 0, 0.9)
      e
    })
    data.frame(y = 1 + as.vector(u) + rnorm(n_g * n_t, 0, 0.4),
               g = g, tim = tim)
  })
  f1 <- frm(bf(y ~ 1 + hetar1(tim + 0 | g)) + gaussian(), data = dd)
  g1 <- suppressWarnings(
    glmmTMB::glmmTMB(y ~ 1 + hetar1(tim + 0 | g), data = dd)
  )
  expect_loglik_equal(f1, g1, tol = 1e-5)

  f2 <- frm(bf(y ~ 1 + homcs(tim + 0 | g)) + gaussian(), data = dd)
  g2 <- suppressWarnings(
    glmmTMB::glmmTMB(y ~ 1 + homcs(tim + 0 | g), data = dd)
  )
  expect_loglik_equal(f2, g2, tol = 1e-5)

  # glmmTMB homtoep false-converges here (like toep before); ours must
  # nest above ar1, whose correlation pattern it generalizes
  f3 <- suppressWarnings(   # PD-boundary probes, same as toep
    frm(bf(y ~ 1 + homtoep(tim + 0 | g)) + gaussian(), data = dd)
  )
  fa <- suppressWarnings(
    frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd)
  )
  expect_gte(as.numeric(logLik(f3)) - as.numeric(logLik(fa)), -1e-5)
})

test_that("spatial exp/gau/mat covstructs work over num_factor(x, y)", {
  set.seed(9)
  np <- 25
  reps <- 8
  xc <- round(runif(np), 2) * 10
  yc <- round(runif(np), 2) * 10
  pos <- num_factor(xc, yc)
  D <- as.matrix(dist(cbind(xc, yc)))
  u <- drop(crossprod(chol(exp(-D / 2) + diag(1e-8, np)), rnorm(np)))
  sp <- data.frame(y = rep(u, reps) + rnorm(np * reps, 2, 0.7),
                   pos = rep(pos, reps),
                   grp = factor(rep(1, np * reps)))

  fe <- frm(bf(y ~ 1 + exp(pos + 0 | grp)) + gaussian(), data = sp)
  fg <- frm(bf(y ~ 1 + gau(pos + 0 | grp)) + gaussian(), data = sp)
  fm <- frm(bf(y ~ 1 + mat(pos + 0 | grp)) + gaussian(), data = sp)
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    ge <- glmmTMB::glmmTMB(y ~ 1 + exp(pos + 0 | grp), data = sp)
    gg <- glmmTMB::glmmTMB(y ~ 1 + gau(pos + 0 | grp), data = sp)
    expect_loglik_equal(fe, ge, tol = 1e-5)
    expect_loglik_equal(fg, gg, tol = 1e-5)
  }
  # matern converges (glmmTMB's does not on this data) and is sane:
  # at least as good as the exponential special case region
  expect_true(is.finite(as.numeric(logLik(fm))))
  expect_gte(as.numeric(logLik(fm)) - as.numeric(logLik(fe)), -0.5)
  # exp() the covstruct did not shadow exp() the function
  fx <- frm(bf(y ~ exp(as.numeric(pos) / 25) + (1 | grp)) + gaussian(),
            data = sp)
  expect_length(fixef(fx)$mu, 2L)
})

test_that("sratio/cratio/acat match direct ML and collapse at K = 2", {
  set.seed(4)
  n <- 300
  x <- rnorm(n)
  y2 <- 1L + (runif(n) < plogis(0.8 * x - 0.3))
  d2 <- data.frame(y = y2, x = x)
  g <- stats::glm(I(y == 2) ~ x, binomial, d2)
  for (fam in list(sratio(), cratio(), acat())) {
    f <- frm(bf(y ~ x) + fam, data = d2)
    expect_lt(abs(as.numeric(logLik(f)) - as.numeric(logLik(g))), 1e-6)
  }

  # K = 4 sratio vs a hand-rolled ML reference
  tau_t <- c(-1, 0.2, 1.4)
  P <- local({
    K1 <- length(tau_t)
    Fm <- stats::plogis(outer(-0.8 * x, tau_t, `+`))
    P <- matrix(0, n, K1 + 1)
    surv <- rep(1, n)
    for (k in seq_len(K1)) {
      P[, k] <- Fm[, k] * surv
      surv <- surv * (1 - Fm[, k])
    }
    P[, K1 + 1] <- surv
    P
  })
  y4 <- apply(P, 1, function(p) sample.int(4, 1, prob = p))
  d4 <- data.frame(y = y4, x = x)
  f4 <- frm(bf(y ~ x) + sratio(), data = d4)
  nll <- function(p) {
    K1 <- 3
    Fm <- stats::plogis(outer(-p[1] * x, p[2:4], `+`))
    Pm <- matrix(0, n, 4)
    surv <- rep(1, n)
    for (k in seq_len(K1)) {
      Pm[, k] <- Fm[, k] * surv
      surv <- surv * (1 - Fm[, k])
    }
    Pm[, 4] <- surv
    -sum(log(Pm[cbind(seq_len(n), y4)]))
  }
  op <- stats::optim(c(0.8, tau_t), nll, method = "BFGS",
                     control = list(reltol = 1e-12))
  expect_lt(abs(as.numeric(logLik(f4)) + op$value), 1e-6)
})

test_that("mo() matches direct ML and predicts monotonically", {
  set.seed(6)
  n <- 400
  inc <- sample(0:3, n, replace = TRUE)
  x <- rnorm(n)
  y <- 1 + c(0, 0.5, 0.8, 1)[inc + 1] * 2 + 0.3 * x + rnorm(n, 0, 0.8)
  dd <- data.frame(y = y, inc = inc, x = x)
  fit <- frm(bf(y ~ x + mo(inc)) + gaussian(), data = dd)

  nll <- function(p) {
    zr <- exp(c(0, p[4:6]))
    cz0 <- c(0, cumsum(zr / sum(zr)))
    -sum(stats::dnorm(y, p[1] + p[2] * x + p[3] * 3 * cz0[inc + 1],
                      exp(p[7]), log = TRUE))
  }
  op <- stats::optim(c(1, 0.3, 0.6, 0, 0, 0, log(0.8)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)

  p_new <- predict(fit, newdata = data.frame(x = 0, inc = 0:3))
  expect_true(all(diff(p_new) >= -1e-8))
  expect_equal(unname(fitted(fit)),
               unname(predict(fit, newdata = dd, type = "response")),
               tolerance = 1e-10)
  ps <- predict(fit, newdata = data.frame(x = 0, inc = 0:3),
                se.fit = TRUE)
  expect_true(all(is.finite(ps$se.fit)))
  # ordered factors work and interactions are refused
  dd$incf <- factor(inc, levels = 0:3, ordered = TRUE)
  ff <- frm(bf(y ~ x + mo(incf)) + gaussian(), data = dd)
  expect_loglik_equal(ff, fit, tol = 1e-6)
  expect_error(frm(bf(y ~ mo(inc) * x) + gaussian(), data = dd),
               "standalone")
})

test_that("influence() flags a distorted group", {
  set.seed(23)
  dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
  dd$y <- 1 + 0.5 * dd$x + rnorm(12, 0, 0.4)[dd$g] + rnorm(120, 0, 0.5)
  dd$y[dd$g == "7"] <- dd$y[dd$g == "7"] + 4 * dd$x[dd$g == "7"]
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  infl <- influence(fit, groups = "g")
  expect_equal(nrow(infl$fixed), 12L)
  cd <- cooks.distance(infl)
  expect_equal(names(which.max(cd)), "7")
  expect_output(print(infl), "cooks_d")
})

test_that("CE prediction intervals and condition sets", {
  set.seed(29)
  dd <- data.frame(x = rnorm(150), f = factor(rep(c("a", "b"), 75)))
  dd$y <- rnorm(150, 1 + 0.5 * dd$x + (dd$f == "b"), 1)
  fit <- frm(bf(y ~ x + f) + gaussian(), data = dd)

  ce_e <- conditional_effects(fit, effects = "x")
  ce_p <- conditional_effects(fit, effects = "x", method = "predict",
                              ndraws = 800)
  # prediction band strictly wider than the epred band
  expect_true(all(ce_p$x$upper__ - ce_p$x$lower__ >
                    ce_e$x$upper__ - ce_e$x$lower__))
  # roughly sigma-wide: half-width near 1.96 * sd(y|x)
  hw <- mean(ce_p$x$upper__ - ce_p$x$estimate__)
  expect_lt(abs(hw - 1.96 * sigma(fit)), 0.35)

  cond <- data.frame(f = c("a", "b"), row.names = c("A", "B"))
  ce_c <- conditional_effects(fit, effects = "x", conditions = cond)
  expect_equal(sort(unique(ce_c$x$cond__)), c("A", "B"))
  expect_equal(nrow(ce_c$x), 200L)
  tmp <- file.path(tempdir(), "frmtmb-ce-cond.pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  expect_no_error(plot(ce_c))
})

test_that("vint()/vreal() reach a custom family", {
  set.seed(31)
  n <- 200
  size <- sample(3:6, n, replace = TRUE)
  x <- rnorm(n)
  y <- rbinom(n, size, plogis(0.3 + 0.6 * x))
  dd <- data.frame(y = y, size = size, x = x)
  fam <- custom_family(
    "vbinom", dpars = "mu", links = list(mu = "logit"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
    },
    type = "discrete"
  )
  fv <- frm(bf(y | vint(size) ~ x) + fam, data = dd)
  fr <- frm(bf(y | trials(size) ~ x) + binomial(), data = dd)
  expect_loglik_equal(fv, fr, tol = 1e-6)
  expect_vector_equal(fixef(fv)$mu, fixef(fr)$mu, tol = 1e-5)
  expect_error(frm(bf(y | vint(x) ~ 1) + fam, data = dd), "integers")
})

test_that("frm_multiple pools by Rubin's rules", {
  set.seed(37)
  n <- 150
  x <- rnorm(n)
  y <- rnorm(n, 1 + 0.5 * x, 1)
  x_mis <- x
  x_mis[sample(n, 30)] <- NA
  # crude stochastic regression imputations
  imps <- lapply(1:4, function(i) {
    xi <- x_mis
    xi[is.na(xi)] <- rnorm(sum(is.na(xi)), mean(x_mis, na.rm = TRUE),
                           stats::sd(x_mis, na.rm = TRUE))
    data.frame(y = y, x = xi)
  })
  mfit <- frm_multiple(bf(y ~ x) + gaussian(), data = imps)
  expect_s3_class(mfit, "frmtmb_multiple")
  tab <- mfit$pooled
  expect_equal(nrow(tab), 3L)
  # pooled se exceeds the average within-imputation se
  within_se <- mean(vapply(mfit$fits,
                           function(f) sqrt(vcov(f)["x", "x"]),
                           numeric(1)))
  expect_gt(tab["x", "se"], within_se)
  expect_true(all(tab$fmi >= 0 & tab$fmi <= 1))
  # identical datasets: pooling reduces to a single fit
  m0 <- frm_multiple(bf(y ~ x) + gaussian(),
                     data = list(imps[[1]], imps[[1]]))
  f0 <- frm(bf(y ~ x) + gaussian(), data = imps[[1]])
  expect_vector_equal(m0$pooled$estimate,
                      c(fixef(f0)$mu, fixef(f0)$sigma), tol = 1e-6)
  expect_output(print(mfit), "Rubin")
})

test_that("wiener diffusion as a custom family (vint decision)", {
  # Navarro-Fuss small-time density; dec = 1 hits the upper boundary
  wiener_lpdf <- function(y, dpars, aterms) {
    dec <- aterms$vint1
    a <- dpars$bs
    w0 <- dpars$bias
    w <- w0 + dec * (1 - 2 * w0)
    v <- dpars$mu * (1 - 2 * dec)
    tt <- (y - dpars$ndt) / a^2
    s <- 0
    for (k in -4:4) {
      s <- s + (w + 2 * k) * exp(-(w + 2 * k)^2 / (2 * tt))
    }
    -v * a * w - v^2 * (y - dpars$ndt) / 2 - 2 * log(a) -
      0.5 * log(2 * pi) - 1.5 * log(tt) + log(s)
  }
  set.seed(41)
  # simulate first-passage times by Euler walk
  sim_wiener <- function(n, v, a, w, ndt, dt = 5e-4) {
    out <- matrix(NA_real_, n, 2)
    for (i in seq_len(n)) {
      z <- a * w
      t <- 0
      while (z > 0 && z < a) {
        z <- z + v * dt + sqrt(dt) * rnorm(1)
        t <- t + dt
      }
      out[i, ] <- c(ndt + t, as.numeric(z >= a))
    }
    out
  }
  sw <- sim_wiener(150, v = 1, a = 1.5, w = 0.5, ndt = 0.25)
  dd <- data.frame(rt = sw[, 1], dec = sw[, 2])
  fam <- custom_family(
    "wiener", dpars = c("mu", "bs", "ndt", "bias"),
    links = list(mu = "identity", bs = "log", ndt = "log",
                 bias = "logit"),
    lpdf = wiener_lpdf, type = "continuous",
    init_dpars = list(
      mu = function(y, aterms) 0.5,
      bs = function(y, aterms) 1,
      ndt = function(y, aterms) min(y) / 2,
      bias = function(y, aterms) 0.5
    )
  )
  fit <- frm(bf(rt | vint(dec) ~ 1) + fam, data = dd)

  # direct ML over the same density
  nll <- function(p) {
    dp <- list(mu = p[1], bs = exp(p[2]), ndt = exp(p[3]),
               bias = stats::plogis(p[4]))
    if (dp$ndt >= min(dd$rt)) return(1e10)
    -sum(wiener_lpdf(dd$rt, dp, list(vint1 = dd$dec)))
  }
  op <- stats::optim(c(0.5, 0, log(min(dd$rt) / 2), 0), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-12, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-4)
  # parameters near the simulation truth (Euler + n = 150: loose)
  expect_lt(abs(fixef(fit)$mu[[1]] - 1), 0.6)
  expect_lt(abs(exp(fixef(fit)$bs[[1]]) - 1.5), 0.4)
})

test_that("the draws surface runs the model machinery per draw", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(43)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.6)[dd$g], 0.8)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 600,
                                    refresh = 0, seed = 1))

  s <- summary(ds)
  expect_true(all(c("mean", "sd", "Rhat") %in% colnames(s)))
  fe <- fixef(ds)
  expect_equal(rownames(fe),
               c("(Intercept)", "x", "sigma_(Intercept)"))
  expect_lt(abs(fe["x", "Estimate"] - fixef(fit)$mu[["x"]]), 0.15)

  vc <- VarCorr(ds)
  expect_true(all(c("estimate", "lwr", "upr") %in% names(vc)))

  ep <- posterior_epred(ds, ndraws = 25)
  expect_equal(dim(ep), c(25L, 80L))
  expect_lt(max(abs(colMeans(ep) - fitted(fit))), 0.5)
  pp <- posterior_predict(ds, ndraws = 25)
  expect_gt(mean(apply(pp, 2, stats::sd)),
            mean(apply(ep, 2, stats::sd)))
  ep2 <- posterior_epred(ds, newdata = data.frame(x = 0:1, g = "3"),
                         ndraws = 10)
  expect_equal(dim(ep2), c(10L, 2L))

  h <- hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")
  expect_s3_class(h, "frmtmb_hypothesis")
  expect_true(h$lwr > 0 && h$upr < 1)
  expect_equal(dim(attr(h, "draws")), c(nrow(ds$draws), 1L))

  if (requireNamespace("posterior", quietly = TRUE)) {
    dm <- posterior::as_draws(ds)
    expect_true(inherits(dm, "draws"))
  }
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    expect_s3_class(bayesplot::pp_check(ds, ndraws = 10), "ggplot")
  }
})
