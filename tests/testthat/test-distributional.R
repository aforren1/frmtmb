sim_disp_data <- function(seed = 202, n = 400) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  g <- factor(rep(seq_len(20), length.out = n))
  data.frame(
    y = rnorm(n, 1 + 2 * x, exp(0.2 + 0.4 * z)),
    x = x, z = z, g = g
  )
}

test_that("gaussian sigma ~ x matches glmmTMB dispformula", {
  skip_if_not_installed("glmmTMB")
  dd <- sim_disp_data()
  fit <- frm(bf(y ~ x, sigma ~ z) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, dispformula = ~z, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
  # glmmTMB's gaussian dispformula models log(sd): coefficients match 1:1
  expect_vector_equal(fixef(fit)$sigma, unname(glmmTMB::fixef(ref)$disp),
                      tol = 1e-4)
})

test_that("negbinomial shape ~ x matches glmmTMB nbinom2 dispformula", {
  skip_if_not_installed("glmmTMB")
  set.seed(303)
  n <- 500
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(0.5 + 0.4 * x)
  shape <- exp(0.8 + 0.5 * z)
  y <- rnbinom(n, size = shape, mu = mu)
  dd <- data.frame(y, x, z)
  fit <- frm(bf(y ~ x, shape ~ z) + negbinomial(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, dispformula = ~z, data = dd,
                          family = glmmTMB::nbinom2)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
  expect_vector_equal(fixef(fit)$shape, unname(glmmTMB::fixef(ref)$disp),
                      tol = 1e-3)
})

test_that("random effects in sigma match a hand-rolled RTMB objective", {
  dd <- sim_disp_data()
  # add a group effect in the dispersion
  set.seed(9)
  sig_re <- rnorm(nlevels(dd$g), 0, 0.3)
  dd$y <- rnorm(nrow(dd), 1 + 2 * dd$x,
                exp(0.2 + 0.4 * dd$z + sig_re[dd$g]))

  fit <- frm(bf(y ~ x, sigma ~ z + (1 | g)) + gaussian(), data = dd)

  # reference: same model written directly against RTMB
  y <- dd$y; x <- dd$x; z <- dd$z; gi <- as.integer(dd$g)
  nll_ref <- function(p) {
    nll <- -sum(RTMB::dnorm(p$u, 0, exp(p$log_sd_u), log = TRUE))
    mu <- p$bm[1] + p$bm[2] * x
    sig <- exp(p$bs[1] + p$bs[2] * z + p$u[gi])
    nll - sum(RTMB::dnorm(y, mu, sig, log = TRUE))
  }
  obj <- RTMB::MakeADFun(
    nll_ref,
    list(bm = c(0, 0), bs = c(0, 0), u = numeric(nlevels(dd$g)),
         log_sd_u = 0),
    random = "u", silent = TRUE
  )
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))

  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  est <- opt$par
  expect_vector_equal(fixef(fit)$mu, est[names(est) == "bm"], tol = 1e-4)
  expect_vector_equal(fixef(fit)$sigma, est[names(est) == "bs"],
                      tol = 1e-4)
  # dispersion RE standard deviation
  vc <- VarCorr(fit)
  expect_lt(abs(sqrt(vc[["sigma: 1 | g"]][1, 1]) -
                  exp(est[names(est) == "log_sd_u"])), 1e-4)
})

test_that("constant dpars are fixed via map", {
  dd <- sim_disp_data(n = 200)
  fit <- frm(bf(y ~ x, sigma = 2) + gaussian(), data = dd)
  # sigma must not be estimated
  expect_false("sigma_(Intercept)" %in% outer_par_names(fit))
  expect_identical(unname(fixef(fit)$sigma), log(2))
  expect_equal(unique(predict(fit, dpar = "sigma", type = "response")), 2)
  # loglik equals a direct fixed-sigma ML fit
  nll <- function(p) -sum(dnorm(dd$y, p[1] + p[2] * dd$x, 2, log = TRUE))
  o <- nlminb(c(0, 0), nll)
  expect_lt(abs(as.numeric(logLik(fit)) - (-o$objective)), 1e-6)
  # and the constant survives summary
  s <- summary(fit)
  expect_true(length(s$fixed_dpars) == 1)
})

test_that("dpar names are validated against the family", {
  expect_error(frm(bf(y ~ x, shape ~ z) + gaussian(),
                      data = NULL, dry_run = "spec"),
               "not available for family")
  expect_error(frm(bf(y ~ x, sigma ~ z, sigma = 1) + gaussian(),
                      data = NULL, dry_run = "spec"),
               "Duplicated dpar")
  expect_error(bf(y ~ x, 3), "named")
})

test_that("distributional model with REs in both mu and sigma fits", {
  dd <- sim_disp_data()
  set.seed(11)
  mu_re <- rnorm(nlevels(dd$g), 0, 1)
  dd$y <- dd$y + mu_re[dd$g]
  fit <- frm(bf(y ~ x + (1 | g), sigma ~ z + (1 | g)) + gaussian(),
                data = dd)
  expect_length(fit$frame$re_blocks, 2)
  expect_identical(vapply(fit$frame$re_blocks, `[[`, "", "dpar"),
                   c("mu", "sigma"))
  d <- diagnose(fit, quiet = TRUE)
  expect_equal(d$convergence, 0)
  expect_true(d$pdHess)
})
