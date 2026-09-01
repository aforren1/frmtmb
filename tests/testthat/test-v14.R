# v0.14: anova labels, se() meta-analysis, proportion responses,
# family additions, ranef condVar, profile control, frm_allfit,
# frm_simulate.

test_that("anova labels distinguish distributional fits", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  f1 <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
            data = sleepstudy)
  f2 <- frm(bf(Reaction ~ Days + (Days | Subject), sigma ~ Days) +
              gaussian(), data = sleepstudy)
  a <- anova(f1, f2)
  expect_equal(nrow(a), 2L)
  expect_false(anyDuplicated(rownames(a)) > 0)
  expect_match(rownames(a)[2], "sigma ~ Days", fixed = TRUE)
})

test_that("se() reproduces fixed- and random-effects meta-analysis", {
  set.seed(41)
  k <- 40
  sei <- runif(k, 0.1, 0.5)
  yi <- rnorm(k, 0.3 + rnorm(k, 0, 0.2), sei)
  dd <- data.frame(yi = yi, sei = sei, obs = factor(seq_len(k)))

  # fixed-effect MA: inverse-variance weighted mean, sigma mapped out
  ff <- frm(bf(yi | se(sei) ~ 1) + gaussian(), data = dd)
  w <- 1 / sei^2
  expect_equal(unname(fixef(ff)$mu), sum(w * yi) / sum(w),
               tolerance = 1e-6)
  expect_equal(sigma(ff), 0)
  expect_equal(unname(sqrt(vcov(ff)[1, 1])), sqrt(1 / sum(w)),
               tolerance = 1e-5)

  # random-effects MA (ML): observation-level RE gives tau^2
  fr <- frm(bf(yi | se(sei) ~ 1 + (1 | obs)) + gaussian(), data = dd)
  nll <- function(p) {
    -sum(stats::dnorm(yi, p[1], sqrt(exp(2 * p[2]) + sei^2), log = TRUE))
  }
  op <- stats::optim(c(0, log(0.2)), nll, method = "BFGS")
  expect_lt(abs(as.numeric(logLik(fr)) + op$value), 1e-4)
  expect_equal(unname(fixef(fr)$mu), op$par[1], tolerance = 1e-3)

  # se(sigma = TRUE): estimated sigma added in quadrature
  fs <- frm(bf(yi | se(sei, sigma = TRUE) ~ 1) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fs)) + op$value), 1e-4)
  expect_gt(sigma(fs), 0)

  expect_error(frm(bf(yi | se(sei) ~ 1) + poisson(), data = dd),
               "gaussian and student")
})

test_that("proportion response with trials() matches the counts form", {
  skip_if_not_installed("lme4")
  data(cbpp, package = "lme4")
  cbpp$prop <- cbpp$incidence / cbpp$size
  f1 <- frm(bf(incidence | trials(size) ~ period + (1 | herd)) +
              binomial(), data = cbpp)
  f2 <- frm(bf(prop | trials(size) ~ period + (1 | herd)) + binomial(),
            data = cbpp)
  expect_loglik_equal(f1, f2, tol = 1e-8)
  expect_vector_equal(fixef(f1)$mu, fixef(f2)$mu, tol = 1e-8)
})

test_that("weibull and exponential match survreg", {
  skip_if_not_installed("survival")
  set.seed(3)
  n <- 300
  x <- rnorm(n)
  dd <- data.frame(x = x, y = rweibull(n, 1.5, exp(1 + 0.5 * x)))
  fw <- frm(bf(y ~ x) + weibull(), data = dd)
  rw <- survival::survreg(survival::Surv(y) ~ x, data = dd,
                          dist = "weibull")
  expect_lt(abs(as.numeric(logLik(fw)) - as.numeric(logLik(rw))), 1e-4)
  # both are log-linear; the mean parameterization shifts the intercept
  expect_equal(fixef(fw)$mu[["x"]], unname(coef(rw)["x"]),
               tolerance = 1e-3)

  fe <- frm(bf(y ~ x) + exponential(), data = dd)
  re <- survival::survreg(survival::Surv(y) ~ x, data = dd,
                          dist = "exponential")
  expect_lt(abs(as.numeric(logLik(fe)) - as.numeric(logLik(re))), 1e-4)
  expect_vector_equal(fixef(fe)$mu, coef(re), tol = 1e-3)
})

test_that("bernoulli and geometric reduce to their parents", {
  set.seed(7)
  n <- 300
  x <- rnorm(n)
  dd <- data.frame(x = x, yb = rbinom(n, 1, plogis(0.3 + 0.8 * x)),
                   yc = rnbinom(n, size = 1, mu = exp(0.5 + 0.4 * x)))
  fb <- frm(bf(yb ~ x) + bernoulli(), data = dd)
  g <- stats::glm(yb ~ x, binomial, dd)
  expect_lt(abs(as.numeric(logLik(fb)) - as.numeric(logLik(g))), 1e-6)
  expect_vector_equal(fixef(fb)$mu, coef(g), tol = 1e-5)
  expect_error(frm(bf(yc ~ x) + bernoulli(), data = dd), "0/1")

  fg <- frm(bf(yc ~ x) + geometric(), data = dd)
  fn <- frm(bf(yc ~ x, shape = 1) + negbinomial(), data = dd)
  expect_loglik_equal(fg, fn, tol = 1e-8)
  expect_vector_equal(fixef(fg)$mu, fixef(fn)$mu, tol = 1e-6)
})

test_that("shifted_lognormal matches a direct ML reference", {
  set.seed(5)
  y <- 0.3 + stats::rlnorm(500, 0.2, 0.5)
  dd <- data.frame(y = y)
  f <- frm(bf(y ~ 1) + shifted_lognormal(), data = dd)
  nll <- function(p) {
    ndt <- exp(p[3])
    if (ndt >= min(y)) return(1e10)
    -sum(stats::dlnorm(y - ndt, p[1], exp(p[2]), log = TRUE))
  }
  op <- stats::optim(c(0.2, log(0.5), log(0.3)), nll,
                     method = "Nelder-Mead",
                     control = list(reltol = 1e-12, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(f)) + op$value), 1e-3)
})

test_that("hurdle and zero-inflated additions match references", {
  set.seed(11)
  n <- 400
  x <- rnorm(n)
  nz <- rbinom(n, 1, 0.7)
  dd <- data.frame(
    x = x,
    yg = nz * rgamma(n, shape = 2, scale = exp(0.5 + 0.3 * x) / 2),
    yl = nz * stats::rlnorm(n, 0.2 + 0.3 * x, 0.6),
    zb = rbinom(n, 1, 0.8) * rbinom(n, 10, plogis(0.2 + 0.5 * x)),
    bt = ifelse(rbinom(n, 1, 0.8) == 1,
                stats::rbeta(n, 2 * 5, (1 - 0.4) * 5 / 0.4), 0)
  )
  dd$bt <- pmin(dd$bt, 0.999)

  fh <- frm(bf(yg ~ x) + hurdle_gamma(), data = dd)
  hu_hat <- unname(plogis(fixef(fh)$hu))
  expect_lt(abs(hu_hat - mean(dd$yg == 0)), 0.05)
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    gt <- glmmTMB::glmmTMB(yg ~ x, ziformula = ~1,
                           family = glmmTMB::ziGamma(link = "log"),
                           data = dd)
    expect_lt(abs(as.numeric(logLik(fh)) - as.numeric(logLik(gt))), 1e-4)
    expect_vector_equal(fixef(fh)$mu, glmmTMB::fixef(gt)$cond,
                        tol = 1e-3)
  }

  fl <- frm(bf(yl ~ x) + hurdle_lognormal(), data = dd)
  nll <- function(p) {
    hu <- plogis(p[4])
    i0 <- dd$yl == 0
    -sum(ifelse(i0, log(hu),
                log(1 - hu) + stats::dlnorm(dd$yl, p[1] + p[2] * x,
                                            exp(p[3]), log = TRUE)))
  }
  op <- stats::optim(c(0, 0, 0, 0), nll, method = "BFGS")
  expect_lt(abs(as.numeric(logLik(fl)) + op$value), 1e-4)

  fz <- frm(bf(zb | trials(10) ~ x) + zero_inflated_binomial(),
            data = dd)
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    gz <- glmmTMB::glmmTMB(cbind(zb, 10 - zb) ~ x, ziformula = ~1,
                           family = stats::binomial(), data = dd)
    expect_lt(abs(as.numeric(logLik(fz)) - as.numeric(logLik(gz))), 1e-4)
  }

  fbt <- frm(bf(bt ~ x) + zero_inflated_beta(), data = dd)
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    gb <- glmmTMB::glmmTMB(bt ~ x, ziformula = ~1,
                           family = glmmTMB::beta_family(), data = dd)
    expect_lt(abs(as.numeric(logLik(fbt)) - as.numeric(logLik(gb))),
              1e-4)
  }
})

test_that("asym_laplace reproduces quantile regression", {
  skip_if_not_installed("quantreg")
  set.seed(13)
  n <- 500
  x <- rnorm(n)
  dd <- data.frame(x = x, y = 1 + 0.5 * x + rnorm(n) * (1 + 0.4 * x^2))
  # the check-loss kink can trigger benign false-convergence warnings
  fq <- suppressWarnings(
    frm(bf(y ~ x, quantile = 0.25) + asym_laplace(), data = dd)
  )
  rq <- quantreg::rq(y ~ x, tau = 0.25, data = dd)
  expect_vector_equal(fixef(fq)$mu, coef(rq), tol = 0.02)
})

test_that("ranef condVar and the tidy data-frame forms", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  r <- ranef(fit, condVar = TRUE)
  S <- attr(r[[1]], "condSD")
  expect_equal(dim(S), dim(r[[1]]))
  expect_true(all(S > 0))

  # frmtmb follows the TMB/glmmTMB convention (uncertainty in the
  # fixed effects propagates into the conditional SDs); lme4's postVar
  # is the narrower purely-conditional quantity
  skip_if_not_installed("glmmTMB")
  gt <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject),
                         data = sleepstudy)
  gdf <- as.data.frame(glmmTMB::ranef(gt))
  expect_vector_equal(S[, "(Intercept)"],
                      gdf$condsd[gdf$term == "(Intercept)"], tol = 0.05)
  expect_vector_equal(S[, "Days"], gdf$condsd[gdf$term == "Days"],
                      tol = 0.05)

  df <- as.data.frame(r)
  expect_named(df, c("grp", "term", "level", "condval", "condsd"))
  expect_equal(nrow(df), 36L)

  vc <- as.data.frame(VarCorr(fit))
  expect_named(vc, c("grp", "var1", "var2", "vcov", "sdcor"))
  expect_equal(nrow(vc), 3L)   # two SDs and one correlation
  expect_equal(vc$sdcor[1]^2, vc$vcov[1], tolerance = 1e-10)
})

test_that("control profile = TRUE reproduces the plain fit", {
  dd <- sim_pois_glmm()
  dd$x2 <- rnorm(nrow(dd))
  form <- bf(y ~ x + x2 + (1 | g)) + poisson()
  f0 <- frm(form, data = dd)
  fp <- frm(form, data = dd, control = frmtmb_control(profile = TRUE))

  # like glmmTMB profile=TRUE / glmer nAGQ=0, this is an approximation
  # (the profiled betas ride inside the Laplace step)
  expect_loglik_equal(fp, f0, tol = 0.1)
  expect_equal(stats::AIC(fp), stats::AIC(f0), tolerance = 1e-4)
  expect_vector_equal(fixef(fp)$mu, fixef(f0)$mu, tol = 0.05)
  s <- summary(fp)$coefficients$mu
  s0 <- summary(f0)$coefficients$mu
  expect_vector_equal(s[, "Std. Error"], s0[, "Std. Error"], tol = 0.02)
  # wald confint under profile covers the remaining outer parameters
  cw <- confint(fp)
  expect_true(all(rownames(cw) %in% rownames(confint(f0))))
  # vcov and hypothesis run off the joint precision
  expect_equal(rownames(vcov(fp)), rownames(vcov(f0)))
  expect_vector_equal(sqrt(diag(vcov(fp))), sqrt(diag(vcov(f0))),
                      tol = 0.02)
  h <- hypothesis(fp, "x - x2")
  expect_equal(h$estimate,
               unname(fixef(fp)$mu["x"] - fixef(fp)$mu["x2"]),
               tolerance = 1e-10)

  expect_error(confint(fp, parm = "x", method = "uniroot"),
               "Unknown parameter|profile")
  expect_error(hypothesis(fp, "x", method = "profile"), "profile")
  expect_error(frm(form, data = dd, REML = TRUE,
                   control = frmtmb_control(profile = TRUE)), "REML")
})

test_that("frm_allfit agrees across optimizers", {
  dd <- sim_pois_glmm(n_g = 10, n_per = 10)
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  af <- frm_allfit(fit)
  expect_s3_class(af, "frmtmb_allfit")
  ok <- !vapply(af$fits, is.null, TRUE)
  expect_gte(sum(ok), 2L)
  # every offered optimizer must actually produce a fit: frm_allfit
  # swallows refit errors, so a broken wrapper is invisible otherwise
  expect_true(all(ok))
  ll <- vapply(af$fits[ok], function(f) as.numeric(logLik(f)), 0)
  expect_lt(diff(range(ll)), 1e-4)
  expect_output(print(af), "logLik spread")
})

test_that("frm_allfit drives the nloptr optimizer", {
  skip_if_not_installed("nloptr")
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  af <- frm_allfit(fit)
  expect_true("nloptr_lbfgs" %in% names(af$fits))
  expect_false(is.null(af$fits$nloptr_lbfgs))
  expect_equal(af$fits$nloptr_lbfgs$opt$convergence, 0L)
  expect_equal(as.numeric(logLik(af$fits$nloptr_lbfgs)),
               as.numeric(logLik(af$fits$nlminb)), tolerance = 1e-4)
})

test_that("frm_simulate simulates de novo and recovers parameters", {
  set.seed(19)
  dd <- data.frame(x = rnorm(300), g = factor(rep(1:15, 20)), y = 0)
  form <- bf(y ~ x + (1 | g)) + gaussian()
  tpl <- frm(form, data = dd, dry_run = "frame")$par_template
  expect_named(tpl, c("beta", "betad", "b", "theta"))

  s <- frm_simulate(form, dd, nsim = 5, seed = 1,
                    newparams = list(beta = c(1, 0.5),
                                     betad = log(0.5),
                                     theta = log(0.8)))
  expect_equal(dim(s), c(300L, 5L))
  # refit one draw: parameters come back
  d2 <- dd
  d2$y <- s[[1L]]
  f <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d2)
  expect_lt(abs(fixef(f)$mu[["x"]] - 0.5), 0.15)
  expect_lt(abs(sqrt(VarCorr(f)[[1]][1, 1]) - 0.8), 0.5)

  # fixed b: identical group structure across draws
  s2 <- frm_simulate(form, dd, nsim = 2, seed = 2,
                     newparams = list(beta = c(1, 0.5), betad = log(0.1),
                                      theta = log(1),
                                      b = seq(-2, 2, length.out = 15)))
  gm1 <- tapply(s2[[1L]] - 1 - 0.5 * dd$x, dd$g, mean)
  expect_lt(max(abs(gm1 - seq(-2, 2, length.out = 15))), 0.2)

  expect_error(frm_simulate(form, dd, newparams = list(beta = 1)),
               "length 2")
  expect_error(frm_simulate(form, dd, newparams = list(zeta = 1)),
               "Unknown newparams")
})
