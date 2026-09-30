# make_conditions() and update_adterms(): brms's helpers, whose output
# is compared with brms's own where brms is installed and with literal
# values where it is not.

test_that("make_conditions() crosses levels and mean +/- sd", {
  d <- data.frame(a = c(0, 2, 4), f = factor(c("u", "v", "u")),
                  s = c("p", "q", "q"), stringsAsFactors = FALSE)
  mc <- make_conditions(d, c("a", "f"))
  expect_identical(names(mc), c("a", "f", "cond__"))
  expect_equal(dim(mc), c(6, 3))
  expect_equal(mc$a, rep(c(0, 2, 4), each = 2))
  expect_identical(as.character(mc$f), rep(c("u", "v"), 3))
  expect_identical(mc$cond__[1:2], c("a = 0 & f = u", "a = 0 & f = v"))
  # one variable without its name is the variable itself, as in brms
  expect_identical(make_conditions(d, "s", incl_vars = FALSE)$cond__,
                   factor(c("p", "q")))
  fit <- frm(bf(a ~ f), family = gaussian(), data = d)
  expect_identical(make_conditions(fit, "f"), make_conditions(d, "f"))
  expect_error(make_conditions(d, "nope"), "`nope` is not a column")
  skip_if_not_installed("brms")
  for (v in list(c("a", "f"), "s", c("f", "a"))) {
    expect_identical(make_conditions(d, v), brms::make_conditions(d, v))
  }
  expect_identical(make_conditions(d, c("a", "f"), digits = 0, sep = "|"),
                   brms::make_conditions(d, c("a", "f"), digits = 0,
                                         sep = "|"))
})

test_that("update_adterms() replaces, adds and removes addition terms", {
  form <- y | trials(size) ~ x
  expect_equal(update_adterms(form, ~ trials(10)), y | trials(10) ~ x)
  expect_equal(update_adterms(form, ~ weights(w)),
               y | trials(size) + weights(w) ~ x)
  expect_equal(update_adterms(form, ~ weights(w), action = "replace"),
               y | weights(w) ~ x)
  expect_equal(update_adterms(y ~ x, ~ trials(10)), y | trials(10) ~ x)
  expect_equal(update_adterms(form, ~ 1, action = "replace"), y ~ x)
  # the environment is the formula's own
  e <- new.env()
  f2 <- stats::as.formula("y | se(s) ~ x", env = e)
  expect_identical(environment(update_adterms(f2, ~ se(t))), e)
  # a bf() formula keeps its other parts
  b <- bf(y | trials(size) ~ x, phi ~ x)
  ub <- update_adterms(b, ~ trials(10))
  expect_s3_class(ub, "frmtmb_formula")
  expect_equal(ub$formula, y | trials(10) ~ x)
  expect_identical(ub$pforms, b$pforms)
  expect_error(update_adterms(~ x, ~ weights(w)), "needs a formula with a")
  expect_error(update_adterms(mvbf(bf(y ~ x), bf(z ~ x)), ~ weights(w)),
               "updates one response formula")
  expect_error(update_adterms(y ~ x, ~ weights(w), action = "add"),
               "action")
  skip_if_not_installed("brms")
  cases <- list(
    list(form, ~ trials(10)),
    list(y | se(s, sigma = TRUE) + weights(w) ~ x + (1 | g), ~ se(s2)),
    list(y | resp_se(s) ~ x, ~ se(s3) + cens(c)),
    list(log(y) | trunc(lb = 0) ~ s(x), ~ trunc(ub = 5))
  )
  for (cs in cases) {
    expect_identical(do.call(update_adterms, cs),
                     do.call(brms::update_adterms, cs))
  }
})
