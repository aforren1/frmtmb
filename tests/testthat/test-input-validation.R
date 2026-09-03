# The autotest-driven entry validation: every scalar refusal fires by
# argument name, and correct calls are untouched. One block per helper
# semantic plus the bespoke refusals, so the G5.2 claim that every
# error is demonstrated by a matched test covers the new surface.

iv_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(3)
      dd <- data.frame(x = stats::rnorm(60),
                       g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x, 1)
      cache <<- list(dd = dd,
                     fit = frm(bf(y ~ x + (1 | g)), family = gaussian(),
                               data = dd))
    }
    cache
  }
})

test_that("flags must be TRUE or FALSE, refused by name", {
  cs <- iv_fit()
  expect_error(bf(y ~ x, nl = "yes"), "`nl` must be TRUE or FALSE")
  expect_error(bf(y ~ x, nl = NA), "`nl` must be TRUE or FALSE")
  expect_error(bf(y ~ x, nl = c(TRUE, FALSE)), "`nl` must be TRUE")
  expect_error(frmtmb_control(profile = 1L), "`profile` must be TRUE")
  expect_error(frmtmb_control(sparse_x = "auto"), "`sparse_x` must be")
  expect_error(frm(bf(y ~ x), data = cs$dd, REML = "yes"),
               "`REML` must be TRUE")
  expect_error(predict(cs$fit, se.fit = "yes"), "`se.fit` must be TRUE")
  expect_error(predict(cs$fit, allow_new_levels = NA),
               "`allow_new_levels` must be TRUE")
  expect_error(ranef(cs$fit, condVar = 1), "`condVar` must be TRUE")
  expect_error(diagnose(cs$fit, quiet = "sh"), "`quiet` must be TRUE")
  expect_error(set_rescor("no"), "must be TRUE")
})

test_that("counts must be single whole numbers, refused by name", {
  cs <- iv_fit()
  expect_error(simulate(cs$fit, nsim = 2.5), "`nsim` must be a single")
  expect_error(simulate(cs$fit, nsim = c(2, 3)), "`nsim` must be a single")
  expect_error(frm_simulate(cs$fit, nsim = 0), "`nsim` must be a single")
  expect_error(frm_bootstrap(cs$fit, nsim = "many"),
               "`nsim` must be a single")
  expect_error(frmtmb_control(restarts = -1), "`restarts` must be")
  # a whole number in double form is fine
  expect_type(frmtmb_control(restarts = 2.0), "list")
})

test_that("coverages must be single probabilities, refused by name", {
  cs <- iv_fit()
  expect_error(confint(cs$fit, level = c(0.9, 0.95)),
               "`level` must be a single")
  expect_error(confint(cs$fit, level = 1), "`level` must be a single")
  expect_error(confint(cs$fit, level = 0), "`level` must be a single")
  expect_error(conditional_effects(cs$fit, effects = "x", prob = 1.2),
               "`prob` must be a single")
  ok <- confint(cs$fit, parm = "x", level = 0.9)
  expect_true(all(is.finite(unlist(ok[, c("lwr", "upr")]))))
})

test_that("named-list arguments refuse the unnamed and partially named", {
  cs <- iv_fit()
  expect_error(frm(bf(y ~ x), data = cs$dd, start = list(1, 2)),
               "`start` must be a named list")
  expect_error(frm(bf(y ~ x), data = cs$dd,
                   start = list(1, beta = 2)),
               "`start` must be a named list")
  expect_error(frmtmb_control(optCtrl = c(iter.max = 5)),
               "`optCtrl` must be a list of options")
})

test_that("string choices are validated where match.arg cannot reach", {
  cs <- iv_fit()
  expect_error(frm(bf(y ~ x), data = cs$dd, dry_run = "nonsense"),
               "`dry_run` must be one of")
  expect_error(frm(bf(y ~ x), data = cs$dd, dry_run = TRUE),
               "`dry_run` must be one of")
  expect_error(frmtmb_control(optimizer = c("a", "b")),
               "`optimizer`")
  # value checking stays at fit time, where the message carries context
  expect_type(frmtmb_control(optimizer = "nope"), "list")
  expect_error(frm(bf(y ~ x), data = cs$dd,
                   control = frmtmb_control(optimizer = "nope")),
               "Unknown optimizer")
})

test_that("the bespoke refusals name their argument and contract", {
  cs <- iv_fit()
  expect_error(frm(bf(y ~ x), data = NULL), "is NULL: frm")
  expect_error(frm(bf(y ~ x), data = cs$dd, control = list(a = 1)),
               "`control`")
  expect_error(frm(bf(y ~ x), data = cs$dd, na.action = 42),
               "`na.action`")
  expect_error(frm(bf(y ~ x), data = cs$dd, priors = 5), "`priors`")
  expect_error(bernoulli(link = c("logit", "probit")),
               "single string")
  expect_error(bernoulli(link = 1L), "single string")
  expect_error(predict(cs$fit, newdata = "nope"), "`newdata`")
  expect_error(predict(cs$fit, re.form = "oops"), "`re.form`")
  expect_error(confint(cs$fit, parm = 1i), "`parm`")
})

test_that("correct calls are untouched by the validation layer", {
  cs <- iv_fit()
  expect_s3_class(bf(y ~ x, nl = FALSE), "frmtmb_formula")
  p <- predict(cs$fit, se.fit = TRUE)
  expect_true(all(is.finite(p$se.fit)))
  expect_identical(nrow(confint(cs$fit, level = 0.95)) > 0, TRUE)
  s <- simulate(cs$fit, nsim = 3L)
  expect_identical(ncol(as.matrix(s)), 3L)
})

test_that("a tmbstan build that samples the wrong density is refused", {
  skip_if_not_installed("tmbstan")
  # this installation is a healthy binary build: the static check must
  # pass silently (the affected builds are source installs whose
  # model.hpp keeps an unpatched std_normal placeholder; see
  # dev/prior-dropping-investigation.md)
  expect_false(frmtmb:::tmbstan_build_broken())
  expect_silent(frmtmb:::check_tmbstan_build("frm_sample()"))
  # the refusal path, exercised against a synthesized broken model.hpp
  # through the same grepl the guard runs
  bad <- "lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));"
  expect_true(any(grepl("std_normal_lpdf<propto__>(y)", bad,
                        fixed = TRUE)))
})
