# emmeans on nonlinear, multivariate and distributional predictors, with
# brms's arguments (dpar, nlpar, resp, epred, re_formula). The numbers
# behind these checks, with a parametric bootstrap, are in
# dev/emm-validate.R and dev/emm-findings.md.

skip_if_not_installed("emmeans")

emm_data <- function(n = 150, seed = 1) {
  set.seed(seed)
  d <- data.frame(f = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = runif(n))
  d$y1 <- 1 + as.numeric(d$f) * 0.5 + 2 * d$x + rnorm(n)
  d$y2 <- 0.5 * as.numeric(d$f) + rnorm(n)
  d$cnt <- rpois(n, exp(0.3 + 0.4 * as.numeric(d$f)))
  d
}

emm_df <- function(e) as.data.frame(summary(e))

# The largest relative difference between two vectors.
emm_rel <- function(a, b) max(abs(a - b) / pmax(abs(b), abs(a)))

test_that("nlpar and whole-mu means match the model fitted linearly", {
  d <- emm_data()
  fnl <- frm(bf(y1 ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d)
  flin <- frm(bf(y1 ~ f + x), data = d)
  # two optimizers stop at slightly different points; every comparison
  # is held to a multiple of that measured floor
  floor <- emm_rel(unname(fixef(fnl)[c("a_Intercept", "a_fB", "a_fC",
                                       "b_Intercept"), "Estimate"]),
                   unname(fixef(flin)[c("Intercept", "fB", "fC", "x"),
                                      "Estimate"]))
  tol <- 50 * max(floor, sqrt(.Machine$double.eps))

  # a is mu at x = 0, and its own predictor is linear in its coefficients
  a_nl <- emm_df(emmeans::emmeans(fnl, "f", nlpar = "a"))
  a_lin <- emm_df(emmeans::emmeans(flin, "f", at = list(x = 0)))
  expect_lt(emm_rel(a_nl$emmean, a_lin$emmean), tol)
  expect_lt(emm_rel(a_nl$SE, a_lin$SE), tol)

  # the whole mu goes through the Jacobian of the body, which is exact
  # for a body that is linear in its parameters
  m_nl <- emmeans::emmeans(fnl, "f")
  m_lin <- emmeans::emmeans(flin, "f")
  expect_lt(emm_rel(emm_df(m_nl)$emmean, emm_df(m_lin)$emmean), tol)
  expect_lt(emm_rel(emm_df(m_nl)$SE, emm_df(m_lin)$SE), tol)
  p_nl <- emm_df(pairs(m_nl))
  p_lin <- emm_df(pairs(m_lin))
  expect_lt(emm_rel(p_nl$SE, p_lin$SE), tol)
  expect_lt(max(abs(p_nl$estimate - p_lin$estimate)), tol * max(p_lin$SE))

  # an identity link makes the expected response the predictor itself
  e_nl <- emm_df(emmeans::emmeans(fnl, "f", epred = TRUE))
  expect_identical(e_nl$emmean, emm_df(m_nl)$emmean)
  expect_lt(emm_rel(e_nl$SE, emm_df(m_nl)$SE), sqrt(.Machine$double.eps))
})

test_that("resp selects one response and no resp stacks them", {
  d <- emm_data()
  fmv <- frm(mvbf(bf(y1 ~ f), bf(y2 ~ f)), data = d)
  fu1 <- frm(bf(y1 ~ f), data = d)
  fu2 <- frm(bf(y2 ~ f), data = d)
  u1 <- emm_df(emmeans::emmeans(fu1, "f"))
  u2 <- emm_df(emmeans::emmeans(fu2, "f"))
  floor <- max(emm_rel(unname(fixef(fmv)[1:3, "Estimate"]),
                       unname(fixef(fu1)[, "Estimate"])),
               emm_rel(unname(fixef(fmv)[4:6, "Estimate"]),
                       unname(fixef(fu2)[, "Estimate"])))
  tol <- 50 * max(floor, sqrt(.Machine$double.eps))

  r1 <- emm_df(emmeans::emmeans(fmv, "f", resp = "y1"))
  expect_lt(emm_rel(r1$emmean, u1$emmean), tol)
  expect_lt(emm_rel(r1$SE, u1$SE), tol)

  # brms's layout: the responses are the levels of rep.meas
  st <- emm_df(emmeans::emmeans(fmv, ~ f | rep.meas))
  expect_identical(levels(st$rep.meas), c("y1", "y2"))
  expect_lt(emm_rel(st$emmean[st$rep.meas == "y2"], u2$emmean), tol)
  # without rescor the two responses' coefficients are independent, so
  # a contrast across responses has the root-sum-of-squares error
  ct <- emm_df(pairs(emmeans::emmeans(fmv, ~ rep.meas | f)))
  expect_lt(emm_rel(ct$SE, sqrt(u1$SE^2 + u2$SE^2)), tol)
})

test_that("the delta-method covariance matches a numerical Jacobian", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  set.seed(11)
  n <- 150
  d <- data.frame(f = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = runif(n, 0, 2))
  d$y <- c(A = 2, B = 3, C = 4)[as.character(d$f)] * exp(-0.8 * d$x) +
    rnorm(n, sd = 0.3)
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ f, b ~ 1, nl = TRUE), data = d)
  rg <- emmeans::ref_grid(fit, at = list(x = c(0.2, 1.6)))
  th <- frmtmb:::interop_coef_vector(fit)
  pred <- function(p) {
    frm_linpred(frmtmb:::set_coef.frmtmb_fit(fit, p), newdata = rg@grid)
  }
  expect_identical(unname(rg@bhat), unname(as.numeric(pred(th))))
  Vc <- frmtmb:::interop_vcov(fit)
  jv <- function(J) J %*% Vc %*% t(J)
  V_rich <- jv(numDeriv::jacobian(pred, th))
  V_fwd <- jv(numDeriv::jacobian(pred, th, method = "simple"))
  # the taped Jacobian sits closer to Richardson's answer than a plain
  # forward difference does, by orders of magnitude
  err_emm <- max(abs(rg@V - V_rich))
  err_fwd <- max(abs(V_fwd - V_rich))
  expect_lt(err_emm, err_fwd / 100)
})

test_that("dpar selects that parameter's predictor and its link", {
  d <- emm_data()
  fit <- frm(bf(y1 ~ f + x, sigma ~ f), data = d)
  e <- emmeans::emmeans(fit, "f", dpar = "sigma")
  s <- emm_df(e)
  nd <- data.frame(f = s$f, x = mean(d$x))
  # this used to be mu's marginal means under a dpar = "sigma" label
  expect_equal(s$emmean,
               unname(frm_linpred(fit, newdata = nd, dpar = "sigma")),
               tolerance = sqrt(.Machine$double.eps))
  r <- as.data.frame(summary(e, type = "response"))
  expect_equal(r$response, exp(s$emmean),
               tolerance = sqrt(.Machine$double.eps))
})

test_that("type = 'response' applies the inverse link once", {
  d <- emm_data()
  fit <- frm(bf(cnt ~ f + x) + poisson(), data = d)
  nd <- data.frame(f = factor(levels(d$f), levels = levels(d$f)),
                   x = mean(d$x))
  mu <- unname(frm_linpred(fit, newdata = nd, type = "response"))
  # the link was not declared to emmeans, so "response" was the log scale
  r <- as.data.frame(summary(emmeans::emmeans(fit, "f"),
                             type = "response"))
  expect_equal(r$rate, mu, tolerance = sqrt(.Machine$double.eps))
  # epred is on the response scale already: no link is applied again
  ep <- emmeans::emmeans(fit, "f", epred = TRUE)
  expect_equal(emm_df(ep)$emmean, mu, tolerance = sqrt(.Machine$double.eps))
  expect_equal(as.data.frame(summary(ep, type = "response"))$emmean,
               emm_df(ep)$emmean)
})

test_that("an addition term on the response is not read as a transform", {
  d <- emm_data()
  d$w <- runif(nrow(d), 0.5, 1.5)
  fit <- frm(bf(y1 | weights(w) ~ f), data = d)
  # emmeans read `y1 | weights(w)` as a transformation "|.weights"
  expect_null(emmeans::emmeans(fit, "f")@misc$tran)
  fl <- frm(bf(log(cnt + 1) ~ f), data = d)
  expect_identical(emmeans::emmeans(fl, "f")@misc$tran, "log")
})

test_that("a smooth or mo() term reaches the marginal means", {
  set.seed(3)
  n <- 200
  d <- data.frame(f = factor(sample(c("A", "B", "C"), n, TRUE)),
                  x = runif(n, 0, 3), e = sample(1:4, n, TRUE))
  d$y <- as.numeric(d$f) + sin(2 * d$x) + 0.5 * d$e + rnorm(n, sd = 0.3)
  fs <- frm(bf(y ~ f + s(x)), data = d)
  s <- emm_df(emmeans::emmeans(fs, "f"))
  nd <- data.frame(f = s$f, x = mean(d$x))
  # the design basis used to leave the smooth out of the means
  expect_equal(s$emmean,
               unname(frm_linpred(fs, newdata = nd, re_formula = NA)),
               tolerance = sqrt(.Machine$double.eps))
  fm <- frm(bf(y ~ f + mo(e)), data = d)
  s <- emm_df(emmeans::emmeans(fm, "f"))
  # the design basis left the mo() contribution out of the means
  lb <- frm_lp_basis(fm, newdata = data.frame(f = s$f, e = mean(d$e)),
                     re_formula = NA)
  expect_equal(s$emmean, unname(lb$eta),
               tolerance = sqrt(.Machine$double.eps))
  expect_equal(s$SE, unname(sqrt(diag(lb$A %*% lb$V %*% t(lb$A)))),
               tolerance = sqrt(.Machine$double.eps))
})

test_that("re_formula = NULL averages over the grouping levels", {
  set.seed(4)
  n <- 180
  d <- data.frame(f = factor(sample(c("A", "B"), n, TRUE)),
                  g = factor(sample(1:6, n, TRUE)))
  d$y <- as.numeric(d$f) + rnorm(6)[d$g] + rnorm(n, sd = 0.5)
  fit <- frm(bf(y ~ f + (1 | g)), data = d)
  e <- emmeans::emmeans(fit, "f", re_formula = NULL)
  nd <- expand.grid(f = levels(d$f), g = levels(d$g))
  nd$f <- factor(nd$f, levels = levels(d$f))
  nd$g <- factor(nd$g, levels = levels(d$g))
  p <- frm_linpred(fit, newdata = nd, re_formula = NULL)
  expect_equal(emm_df(e)$emmean, as.numeric(tapply(p, nd$f, mean)),
               tolerance = sqrt(.Machine$double.eps))
  expect_identical(e@misc$avgd.over, "g")
})

test_that("each refusal reaches the user with its reason", {
  d <- emm_data()
  fnl <- frm(bf(y1 ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d)
  # emmeans used to replace every one of these with "Perhaps a 'data'
  # or 'params' argument is needed"
  expect_error(emmeans::emmeans(fnl, "f", nlpar = "zz"),
               "Non-linear parameter 'zz' is not part of the model")
  expect_error(emmeans::emmeans(fnl, "f", nlpar = "a", dpar = "sigma"),
               "'dpar' and 'nlpar' cannot be specified at the same time")
  expect_error(emmeans::emmeans(fnl, "f", dpar = "nu"),
               "Distributional parameter 'nu' is not part of the model")
  expect_error(emmeans::emmeans(fnl, "f", epred = NA),
               "`epred` must be TRUE or FALSE")
  expect_error(emmeans::emmeans(fnl, "f", re_formula = ~ (1 | zz)),
               "matches no group-level term")

  fmv <- frm(mvbf(bf(y1 ~ f), bf(cnt ~ f) + poisson()), data = d)
  expect_error(emmeans::emmeans(fmv, "f", resp = "y3"),
               "Invalid argument 'resp'. Valid response variables are")
  expect_error(emmeans::emmeans(fmv, "f"), "have different links")
  # the same fit answers on the response scale, where the links do not
  # matter
  ep <- emm_df(emmeans::emmeans(fmv, ~ f | rep.meas, epred = TRUE))
  expect_identical(nrow(ep), 6L)

  d$ord <- cut(d$y1, c(-Inf, 2, 3, Inf), labels = FALSE)
  fo <- frm(bf(ord ~ f) + cumulative(), data = d)
  expect_error(emmeans::emmeans(fo, "f", epred = TRUE),
               "predicts a distribution over the categories")
})

test_that("an exact gp() off its fitted positions is refused by name", {
  skip_on_cran()
  set.seed(5)
  n <- 40
  d <- data.frame(f = factor(rep(c("A", "B"), n / 2)),
                  x = round(runif(n, 0, 3), 1))
  d$y <- as.numeric(d$f) + sin(d$x) + rnorm(n, sd = 0.2)
  fit <- frm(bf(y ~ f + gp(x)), data = d)
  expect_error(emmeans::emmeans(fit, "f"), "exact gp[(][)] term")
  # at an observed position there is no kriging variance to omit
  s <- emm_df(emmeans::emmeans(fit, "f", at = list(x = d$x[1])))
  expect_true(all(is.finite(s$SE)))
})
