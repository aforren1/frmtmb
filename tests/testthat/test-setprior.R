# brms-style set_prior() specification.

fit_sp <- local({
  set.seed(401)
  dd <- data.frame(x = rnorm(150), f = factor(rep(c("a", "b"), 75)),
                   g = factor(rep(1:15, 10)))
  dd$y <- rnorm(150, 1 + 0.5 * dd$x + (dd$f == "b") +
                  rnorm(15, 0, 0.7)[dd$g], 1)
  frm(bf(y ~ x + f + (1 | g), sigma ~ x) + gaussian(), data = dd)
})

test_that("set_prior parses and combines", {
  pl <- set_prior("normal(0, 5)", class = "b") +
    set_prior("student_t(3, 0, 2)", class = "Intercept") +
    set_prior("exponential(1)", class = "sd", lb = 0.01)
  expect_s3_class(pl, "frmtmb_priorlist")
  expect_length(unclass(pl), 3)
  expect_error(set_prior("bogus(1)"), "Unsupported prior")
  expect_error(set_prior(""), "distribution, bounds, or both")
  expect_output(print(pl), "normal")
})

test_that("class targeting resolves correctly", {
  ri <- frmtmb:::resolve_priorlist(fit_sp,
    set_prior("normal(0, 5)", class = "b"))
  # b: non-intercept location coefs only (x, fb) - not sigma's
  expect_length(ri$entries, 2)
  expect_true(all(vapply(ri$entries, `[[`, "", "comp") == "beta"))

  ri2 <- frmtmb:::resolve_priorlist(fit_sp,
    set_prior("normal(0, 2)", class = "b", dpar = "sigma"))
  expect_length(ri2$entries, 1)
  expect_identical(ri2$entries[[1]]$comp, "betad")

  ri3 <- frmtmb:::resolve_priorlist(fit_sp,
    set_prior("exponential(1)", class = "sd"))
  expect_length(ri3$entries, 1)
  expect_identical(ri3$entries[[1]]$scale, "sd")

  # coefficient-specific overrides class-wide (later wins)
  ri4 <- frmtmb:::resolve_priorlist(fit_sp,
    set_prior("normal(0, 5)", class = "b") +
      set_prior("normal(0, 0.1)", class = "b", coef = "x"))
  x_entry <- Filter(function(e) {
    e$comp == "beta" &&
      names(fit_sp$frame$par_template$beta)[e$idx] == "x"
  }, ri4$entries)[[1]]
  expect_equal(x_entry$dist$scale, 0.1)

  expect_error(frmtmb:::resolve_priorlist(fit_sp,
    set_prior("normal(0, 1)", class = "b", coef = "zzz")),
    "not found")
  expect_error(frmtmb:::resolve_priorlist(fit_sp,
    set_prior("normal(0, 1)", class = "sd", group = "zzz")),
    "No random-effect SDs")
})

test_that("sd-class priors act on the natural scale with the Jacobian", {
  ri <- frmtmb:::resolve_priorlist(fit_sp,
    set_prior("exponential(1)", class = "sd"))
  e <- ri$entries[[1]]
  th <- 0.4   # internal log-sd value
  # log p_theta(th) = log p_sd(exp(th)) + th
  expect_equal(as.numeric(frmtmb:::prior_logdens(th, e$dist, e$scale)),
               stats::dexp(exp(th), 1, log = TRUE) + th,
               tolerance = 1e-12)
})


test_that("RTMBdist families: beta_binomial matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(402)
  n <- 400
  x <- rnorm(n)
  mu <- plogis(0.2 + 0.6 * x)
  phi <- 6
  p <- rbeta(n, mu * phi, (1 - mu) * phi)
  dd <- data.frame(y = rbinom(n, 10, p), x = x, size = 10)
  fit <- frm(bf(y | trials(size) ~ x) + beta_binomial(), data = dd)
  ref <- glmmTMB::glmmTMB(cbind(y, size - y) ~ x,
                          family = glmmTMB::betabinomial, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-5)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-3)
})

test_that("RTMBdist families: skew_normal, inverse.gaussian, exgaussian", {
  set.seed(403)
  n <- 800
  x <- rnorm(n)

  y_sn <- RTMBdist::rskewnorm2(n, 1 + 0.5 * x, 1.2, 4)
  d1 <- data.frame(y = y_sn, x = x)
  f1 <- suppressWarnings(frm(bf(y ~ x) + skew_normal(), data = d1))
  expect_vector_equal(fixef(f1)$mu, c(1, 0.5), tol = 0.15)
  # alpha is weakly identified in the skew normal (its ML is known to
  # drift); assert the detected skew direction only
  expect_gt(fixef(f1)$alpha[[1]], 1)

  y_ig <- RTMBdist::rinvgauss(n, exp(0.5 + 0.3 * x), 2)
  d2 <- data.frame(y = y_ig, x = x)
  f2 <- frm(bf(y ~ x) + inverse.gaussian(link = "log"), data = d2)
  expect_vector_equal(fixef(f2)$mu, c(0.5, 0.3), tol = 0.1)
  expect_lt(abs(exp(fixef(f2)$shape[[1]]) - 2), 0.4)

  y_exg <- RTMBdist::rexgauss(n, 2 - 1, 0.5, 1)   # mean 2, beta 1
  d3 <- data.frame(y = y_exg, x = rnorm(n))
  f3 <- frm(bf(y ~ 1) + exgaussian(), data = d3)
  expect_lt(abs(fixef(f3)$mu[[1]] - 2), 0.15)
  expect_lt(abs(exp(fixef(f3)$beta[[1]]) - 1), 0.3)
  # simulators round-trip
  expect_length(simulate(f3, nsim = 1, seed = 1)$sim_1, n)
})