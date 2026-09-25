# A prior with no `resp` in a multivariate model, and one with no
# `nlpar` on a nonlinear location, against brms 2.23.0's stancode() on
# the same models (dev/mvprior-log/brms-probe.txt, script
# dev/mvprior-brms-probe.R, data dev/mvprior-cases.R). brms keys these
# rows by response and by nonlinear parameter, and refuses a
# specification whose key is empty where the model has none empty:
#
#   bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)):
#     class b, b coef x, Intercept, sigma, bounds-only b, no resp: refused
#     class b, resp = "y1":  lprior += normal_lpdf(b_y1 | 0, 5)
#   the same refusal for nu, b and Intercept with dpar = "sigma", ar, ma,
#   cosy and cortime, and in a mi() model
#   bf(ynl ~ a * exp(-b * ax), a ~ 1 + (1 | g), b ~ 1, nl = TRUE):
#     class b, b coef Intercept, Intercept, bounds-only b, no nlpar:
#     refused
#   rescor and cor with no resp: accepted
#
# Through 0.62.0 frmtmb applied the multivariate ones to every response
# that had the slot, class b on the nonlinear model to nothing, and
# class Intercept there to every nonlinear parameter's intercept.

mvr_data <- function() {
  set.seed(2309)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(1:8, length.out = n)),
                  t = rep(1:10, each = 8))
  d$y1 <- 1 + d$x + stats::rnorm(8)[d$g] + stats::rnorm(n)
  d$y2 <- -1 + 0.5 * d$x + stats::rnorm(8)[d$g] + stats::rnorm(n)
  d$cnt <- stats::rpois(n, exp(0.5 + 0.3 * d$x))
  d$ax <- abs(d$x)
  d$ynl <- (2 + stats::rnorm(8, 0, 0.3)[d$g]) * exp(-0.5 * d$ax) +
    stats::rnorm(n, 0, 0.1)
  d
}

# the internal parameters a prior puts a density or a bound on
mvr_reached <- function(formula, data, pl) {
  des <- frmtmb:::prior_design(formula, data, NULL, list())
  r <- frmtmb:::resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  nm <- unlist(lapply(r$entries, function(e) {
    frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
  }))
  sort(c(nm, names(r$lower), names(r$upper)))
}

mvr_msg <- "names no parameter of a multivariate model"

test_that("multivariate: class b and Intercept need resp", {
  d <- mvr_data()
  f <- bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)) + set_rescor(FALSE)
  for (p in list(set_prior("normal(0, 5)", class = "b"),
                 set_prior("normal(0, 5)", class = "b", coef = "x"),
                 set_prior("normal(0, 5)", class = "Intercept"),
                 set_prior("", class = "b", lb = 0))) {
    expect_error(validate_prior(p, f, data = d), mvr_msg,
                 class = "frmtmb_error")
  }
})

test_that("a distributional parameter's class with no resp is refused", {
  d <- mvr_data()
  f <- bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(FALSE)
  expect_error(validate_prior(set_prior("student_t(3, 0, 2.5)",
                                        class = "sigma"), f, data = d),
               mvr_msg, class = "frmtmb_error")
  expect_error(validate_prior(set_prior("gamma(2, 0.1)", class = "nu"),
                              f + student(), data = d),
               mvr_msg, class = "frmtmb_error")
})

# One per block: 0.62.0 refused these two with the shape message, judged
# against y2's sigma, so each errors there rather than failing and would
# hide any assertion after it.
test_that("class b with dpar and no resp is refused, not shape-judged", {
  d <- mvr_data()
  fs <- bf(y1 ~ x, sigma ~ z) + bf(y2 ~ x) + set_rescor(FALSE)
  expect_error(validate_prior(set_prior("normal(0, 5)", class = "b",
                                        dpar = "sigma"), fs, data = d),
               mvr_msg, class = "frmtmb_error")
})

test_that("class Intercept with dpar and no resp: refused, not shape-judged", {
  d <- mvr_data()
  fs <- bf(y1 ~ x, sigma ~ z) + bf(y2 ~ x) + set_rescor(FALSE)
  expect_error(validate_prior(set_prior("normal(0, 5)", class = "Intercept",
                                        dpar = "sigma"), fs, data = d),
               mvr_msg, class = "frmtmb_error")
})

test_that("a residual-correlation class with no resp is refused", {
  d <- mvr_data()
  d$t4 <- rep(1:4, length.out = nrow(d))
  d$g4 <- factor(rep(seq_len(nrow(d) / 4), each = 4))
  far <- bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)) + bf(y2 ~ x) +
    set_rescor(FALSE)
  expect_error(validate_prior(set_prior("normal(0, 0.5)", class = "ar"),
                              far, data = d),
               mvr_msg, class = "frmtmb_error")
  fma <- bf(y1 ~ x + ma(time = t, gr = g, cov = TRUE)) +
    bf(y2 ~ x + cosy(time = t, gr = g)) + set_rescor(FALSE)
  expect_error(validate_prior(set_prior("normal(0, 0.5)", class = "ma"),
                              fma, data = d),
               mvr_msg, class = "frmtmb_error")
  expect_error(validate_prior(set_prior("beta(2, 2)", class = "cosy"),
                              fma, data = d),
               mvr_msg, class = "frmtmb_error")
  fun <- bf(y1 ~ x + unstr(time = t4, gr = g4)) + bf(y2 ~ x) +
    set_rescor(FALSE)
  expect_error(validate_prior(set_prior("lkj(2)", class = "cortime"),
                              fun, data = d),
               mvr_msg, class = "frmtmb_error")
})

test_that("frm() and a brms prior table meet the same refusal", {
  d <- mvr_data()
  f <- bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(FALSE) + gaussian()
  expect_error(frm(f, data = d, prior = set_prior("normal(0, 5)",
                                                  class = "b")),
               mvr_msg, class = "frmtmb_error")
  # one bad specification among good ones is enough
  expect_error(frm(f, data = d,
                   prior = set_prior("normal(0, 5)", class = "Intercept") +
                     set_prior("normal(0, 5)", class = "b", resp = "y1")),
               mvr_msg, class = "frmtmb_error")
  skip_if_not_installed("brms")
  bp <- brms::set_prior("normal(0, 5)", class = "b")
  expect_error(frm(f, data = d, prior = bp), mvr_msg,
               class = "frmtmb_error")
})

test_that("the refusal names the responses that have the slot", {
  d <- mvr_data()
  fm <- bf(y1 ~ x, family = gaussian()) + bf(cnt ~ x, family = poisson()) +
    set_rescor(FALSE)
  e <- tryCatch(validate_prior(set_prior("student_t(3, 0, 2.5)",
                                         class = "sigma"), fm, data = d),
                error = identity)
  expect_s3_class(e, "frmtmb_error")
  # read without conditionMessage(), which errors on the table 0.62.0
  # returns and would stop the block before the assertions below
  msg <- if (inherits(e, "error")) conditionMessage(e) else ""
  # a poisson response has no sigma, so only y1 is offered
  expect_match(msg, "resp = \"y1\")", fixed = TRUE)
  expect_no_match(msg, "resp = \"cnt\"", fixed = TRUE)
  # the advice resolves as written
  expect_identical(mvr_reached(fm, d, set_prior("student_t(3, 0, 2.5)",
                                                class = "sigma",
                                                resp = "y1")),
                   "y1_sigma_(Intercept)")
})

test_that("nonlinear location: class b and Intercept need nlpar", {
  d <- mvr_data()
  f <- bf(ynl ~ a * exp(-b * ax), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)
  for (p in list(set_prior("normal(0, 5)", class = "b"),
                 set_prior("normal(0, 5)", class = "b", coef = "Intercept"),
                 set_prior("normal(0, 5)", class = "Intercept"),
                 set_prior("", class = "b", lb = 0))) {
    expect_error(validate_prior(p, f, data = d), "location is nonlinear",
                 class = "frmtmb_error")
  }
})

# The accepted spellings. These pass on 0.62.0 too; they guard against
# the refusals above reaching further than brms's do.
test_that("a prior with resp reaches that response only", {
  d <- mvr_data()
  f <- bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)) + set_rescor(FALSE)
  expect_identical(mvr_reached(f, d, set_prior("normal(0, 5)", class = "b",
                                               resp = "y1")), "y1_x")
  expect_identical(mvr_reached(f, d, set_prior("normal(0, 5)",
                                               class = "Intercept",
                                               resp = "y2")),
                   "y2_(Intercept)")
  expect_identical(mvr_reached(f, d, set_prior("", class = "b", resp = "y1",
                                               lb = 0)), "y1_x")
  expect_identical(mvr_reached(f, d, set_prior("student_t(3, 0, 2.5)",
                                               class = "sigma",
                                               resp = "y2")),
                   "y2_sigma_(Intercept)")
})

test_that("rescor, cor and theta need no resp", {
  d <- mvr_data()
  f <- bf(mvbind(y1, y2) ~ x + (1 + x | g)) + set_rescor(TRUE)
  expect_identical(mvr_reached(f, d, set_prior("lkj(2)", class = "rescor")),
                   "thetar_1")
  expect_length(mvr_reached(f, d, set_prior("lkj(2)", class = "cor")), 2L)
  expect_identical(mvr_reached(f, d, set_prior("normal(0, 1)",
                                               class = "theta",
                                               coef = "theta_1")),
                   "theta_1")
})

test_that("a univariate or nonlinear model keeps its prefix-less spellings", {
  d <- mvr_data()
  d$y <- d$y1
  fu <- bf(y ~ x + (1 | g), sigma ~ z)
  expect_identical(mvr_reached(fu, d, set_prior("normal(0, 5)",
                                                class = "b")), "x")
  expect_identical(mvr_reached(fu, d, set_prior("normal(0, 5)",
                                                class = "Intercept")),
                   "(Intercept)")
  expect_identical(mvr_reached(fu, d, set_prior("normal(0, 5)", class = "b",
                                                dpar = "sigma")), "sigma_z")
  expect_identical(mvr_reached(bf(y ~ x), d,
                               set_prior("student_t(3, 0, 2.5)",
                                         class = "sigma")),
                   "sigma_(Intercept)")
  fnl <- bf(ynl ~ a * exp(-b * ax), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)
  expect_identical(mvr_reached(fnl, d, set_prior("normal(0, 5)",
                                                 nlpar = "a")),
                   "a_(Intercept)")
  expect_identical(mvr_reached(fnl, d, set_prior("student_t(3, 0, 2.5)",
                                                 class = "sigma")),
                   "sigma_(Intercept)")
})
