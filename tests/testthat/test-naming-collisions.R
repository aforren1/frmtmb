# Reserved-name shadowing in hypothesis(), and the bare-nlpar spelling
# of bounds and confint(parm =). Both are about one name meaning two
# things; see the "Naming collisions" section of vignette("inputs").

sim_shadow <- function(seed = 3, nm = "sigma") {
  set.seed(seed)
  n <- 200
  dd <- data.frame(v = rnorm(n), g = factor(rep(1:10, length.out = n)))
  dd$y <- rnorm(n, 1 + 0.7 * dd$v + rnorm(10, 0, 0.5)[dd$g], 1)
  names(dd)[1L] <- nm
  dd
}

test_that("a covariate named sigma shadows the residual SD, once, by name", {
  dd <- sim_shadow()
  fit <- frm(bf(y ~ sigma + (1 | g)) + gaussian(), data = dd)

  # one message per call, naming both meanings and the winner
  msgs <- capture_messages(h <- hypothesis(fit, c("sigma = 0", "sigma > 0")))
  expect_length(msgs, 1L)
  expect_match(msgs, "'sigma'", fixed = TRUE)
  expect_match(msgs, "residual standard deviation")
  expect_match(msgs, "'.sigma'", fixed = TRUE)

  # the coefficient is what was tested, both rows of it
  expect_equal(h$estimate, rep(unname(fixef(fit)$mu[["sigma"]]), 2L),
               tolerance = 1e-10)

  # ... and the shadowed quantity is reachable under the dot spelling
  expect_equal(suppressMessages(hypothesis(fit, ".sigma"))$estimate,
               unname(sigma(fit)), tolerance = 1e-8)

  # both names are listed, so the escape hatch is discoverable
  vv <- variables(fit)
  expect_true(all(c("sigma", ".sigma") %in% vv))

  # the note is per call, not per session
  expect_length(capture_messages(hypothesis(fit, "sigma = 0")), 1L)
})

test_that("a coefficient shadowing an sd_ name gets the same treatment", {
  dd <- sim_shadow(nm = "sd_g__Intercept")
  fit <- frm(bf(y ~ sd_g__Intercept + (1 | g)) + gaussian(), data = dd)

  msgs <- capture_messages(hypothesis(fit, "sd_g__Intercept = 0"))
  expect_length(msgs, 1L)
  expect_match(msgs, "random-effect standard deviation")
  expect_match(msgs, "'.sd_g__Intercept'", fixed = TRUE)

  # the dot name is the natural-scale standard deviation itself
  vc <- confint_varcorr(fit)
  expect_equal(suppressMessages(hypothesis(fit, ".sd_g__Intercept"))$estimate,
               vc$estimate[vc$type == "sd"][1L], tolerance = 1e-8)
})

test_that("a clean model says nothing and grows no dot names", {
  dd <- sim_shadow(nm = "v")
  fit <- frm(bf(y ~ v + (1 | g)) + gaussian(), data = dd)
  expect_no_message(hypothesis(fit, c("v = 0", "sigma > 0")))
  expect_false(any(startsWith(variables(fit), ".")))
  # sigma still means the residual SD when nothing has taken the name
  expect_equal(hypothesis(fit, "sigma")$estimate, unname(sigma(fit)),
               tolerance = 1e-8)
})

# --- bare nonlinear-parameter names in bounds and confint(parm =) -----

sim_nl <- function(seed = 1, n = 120) {
  set.seed(seed)
  dd <- data.frame(x = runif(n, 0, 3))
  dd$y <- rnorm(n, 2.5 * exp(-0.8 * dd$x), 0.2)
  dd$g <- factor(rep(1:4, length.out = n))
  dd
}

test_that("an intercept-only nlpar may be bounded by its bare name", {
  dd <- sim_nl()
  form <- bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)
  fit <- frm(form + gaussian(), data = dd)
  nm <- outer_par_names(fit)

  bare <- resolve_bounds(fit, c(b = 0.1), c(b = 2))
  expect_identical(bare$lower[nm == "b_(Intercept)"], 0.1)
  expect_identical(bare$upper[nm == "b_(Intercept)"], 2)
  expect_true(all(is.infinite(bare$lower[nm != "b_(Intercept)"])))

  # the full and parenthesis-free spellings are unchanged, and the bare
  # name is exactly equivalent to them
  full <- resolve_bounds(fit, c(`b_(Intercept)` = 0.1), c(`b_(Intercept)` = 2))
  bareless <- resolve_bounds(fit, c(b_Intercept = 0.1), c(b_Intercept = 2))
  expect_identical(bare, full)
  expect_identical(bare, bareless)

  # and it reaches frm() itself, which is where a user writes it
  bounded <- frm(form + gaussian(), data = dd,
                 prior = set_prior("", nlpar = "b", lb = 0.1, ub = 2))
  expect_equal(unname(fixef(bounded)$b), unname(fixef(fit)$b),
               tolerance = 1e-6)

  # confint(parm =) takes it too, silently: it is a spelling of one
  # internal parameter, not a change of scale
  expect_no_message(ci <- confint(fit, parm = "b", method = "wald"))
  expect_identical(rownames(ci), "b_(Intercept)")
  expect_equal(ci, confint(fit, parm = "b_(Intercept)", method = "wald"))
})

test_that("a bare nlpar with several coefficients is refused by name", {
  dd <- sim_nl()
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1 + g, b ~ 1, nl = TRUE) +
               gaussian(), data = dd)
  expect_error(resolve_bounds(fit, c(a = 0), NULL),
               "more than one coefficient")
  expect_error(resolve_bounds(fit, c(a = 0), NULL), "a_g2")
  expect_error(confint(fit, parm = "a"), "more than one coefficient")
  # the unambiguous one in the same model still resolves
  expect_no_error(resolve_bounds(fit, c(b = 0), NULL))
})

test_that("the unknown-parameter errors advertise the bare spelling", {
  dd <- sim_nl()
  fit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) + gaussian(),
             data = dd)
  expect_error(resolve_bounds(fit, c(zzz = 0), NULL),
               "intercept-only nonlinear parameters may be named bare")
  expect_error(confint(fit, parm = "zzz"),
               "intercept-only nonlinear parameters may be named bare")
  # a dpar is deliberately NOT aliased: bare `sigma` is the natural-scale
  # summary elsewhere, and bounds are on the internal scale
  expect_error(resolve_bounds(fit, c(sigma = 0), NULL),
               "Unknown parameter\\(s\\) in bounds: sigma")
})
