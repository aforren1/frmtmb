# par_template(): the parameter vocabulary, before and after fitting,
# and the named form of frm(start =) that it feeds.
#
# The positional contract is characterized here rather than changed;
# test-methods.R ("start values are validated") and test-nl.R already
# depend on it, and both still pass.

# This file tests the parameter VOCABULARY, not a likelihood, so the
# response is drawn from the model rather than written out beside it.
pt_data <- function() {
  set.seed(11)
  dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)), y = 0)
  # The draw takes its own seed, away from the fixture's: reusing
  # the fixture seed restarts the same random stream that made the
  # covariates, and the residuals come out equal to x.
  dd$y <- frm_simulate(bf(y ~ x) + gaussian(), dd,
                       newparams = list(Intercept = 1, x = 2, sigma = 1),
                       nsim = 1, seed = 1011)[[1]]
  dd
}

pt_frame <- function(dd) {
  frm(bf(y ~ x) + gaussian(), dd, dry_run = "frame")
}

test_that("par_template() on a formula returns the layout and the starts", {
  dd <- pt_data()
  tp <- par_template(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_s3_class(tp, "frmtmb_par_template")
  expect_true(is.list(tp))
  # every component is a NAMED numeric vector, whatever the frame stored
  for (nm in names(tp)) {
    expect_type(tp[[nm]], "double")
    expect_true(all(nzchar(names(tp[[nm]]))))
  }
  expect_named(tp[["beta"]], c("(Intercept)", "x"))
  # the covariance components take the confint() spelling
  expect_named(tp[["theta"]], "theta_1")

  # the values ARE the ones frm() would start from
  frame <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd, dry_run = "frame")
  expect_equal(unclass(unname(tp)), unname(make_start(frame, NULL)),
               ignore_attr = TRUE)
})

test_that("par_template() on a fit returns the estimates", {
  dd <- pt_data()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), dd)
  tp <- par_template(fit)
  expect_s3_class(tp, "frmtmb_par_template")
  expect_equal(unname(tp[["beta"]]), unname(fit$estimates$beta))
  expect_equal(tp[["beta"]][["x"]], unname(fixef(fit)$mu[["x"]]))
  expect_true(attr(tp, "fitted"))
  expect_false(attr(par_template(bf(y ~ x) + gaussian(), data = dd),
                    "fitted"))
})

test_that("par_template() prints its components and refuses a bare formula", {
  dd <- pt_data()
  tp <- par_template(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  out <- utils::capture.output(print(tp))
  expect_true(any(grepl("starting values", out, fixed = TRUE)))
  expect_true(any(grepl("$beta", out, fixed = TRUE)))
  expect_true(any(grepl("(Intercept)", out, fixed = TRUE)))
  # long components are truncated rather than dumped
  short <- utils::capture.output(print(tp, n = 2L))
  expect_true(any(grepl("more", short, fixed = TRUE)))
  expect_output(print(par_template(frm(bf(y ~ x) + gaussian(), dd))),
                "estimates")
  expect_error(par_template(bf(y ~ x) + gaussian()), "needs `data`")
})

test_that("a discovered template round-trips as start=", {
  dd <- pt_data()
  frame <- pt_frame(dd)
  st <- par_template(bf(y ~ x) + gaussian(), data = dd)
  st$beta["x"] <- 2
  # make_start() reproduces the edited template exactly: the fit begins
  # where the user put it
  expect_equal(make_start(frame, st), unclass(st), ignore_attr = TRUE)
  fit <- frm(bf(y ~ x) + gaussian(), dd, start = st)
  expect_s3_class(fit, "frmtmb_fit")
  # and a fit's own template feeds straight back
  expect_silent(frm(bf(y ~ x) + gaussian(), dd,
                    start = par_template(fit)))
})

test_that("a named start component matches by name, both spellings", {
  dd <- pt_data()
  frame <- pt_frame(dd)
  paren <- make_start(frame, list(beta = c("(Intercept)" = 5)))
  bare <- make_start(frame, list(beta = c(Intercept = 5)))
  expect_equal(paren, bare)
  expect_equal(unname(paren$beta[["(Intercept)"]]), 5)
  # a partial override leaves everything else at its default
  cold <- make_start(frame, NULL)
  expect_equal(paren$beta[["x"]], cold$beta[["x"]])
  expect_equal(paren$betad, cold$betad)
})

test_that("named matching reaches every component, not just beta", {
  dd <- pt_data()
  frame <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), dd,
               dry_run = "frame")
  st <- make_start(frame, list(betad = c("sigma_x" = 0.4),
                               theta = c(theta_1 = -1)))
  cold <- make_start(frame, NULL)
  expect_equal(unname(st$betad[["sigma_x"]]), 0.4)
  expect_equal(st$betad[["sigma_(Intercept)"]],
               cold$betad[["sigma_(Intercept)"]])
  expect_equal(unname(st$theta[[1L]]), -1)
})

test_that("a named start refuses what it cannot resolve", {
  dd <- pt_data()
  frame <- pt_frame(dd)
  expect_error(make_start(frame, list(beta = c(zzz = 1))),
               "names no parameter of this model")
  expect_error(make_start(frame, list(beta = c(Intercept = 1, 2))),
               "mixes named and unnamed entries")
  expect_error(make_start(frame,
                          list(beta = c(Intercept = 1,
                                        "(Intercept)" = 2))),
               "addresses one parameter more than once")
  expect_error(make_start(frame, list(beta = "a")),
               "must be a numeric vector")
  # the component check is unchanged
  expect_error(make_start(frame, list(bogus = 1)),
               "Unknown start component")
})

test_that("an unnamed start keeps the positional contract exactly", {
  dd <- pt_data()
  frame <- pt_frame(dd)
  st <- make_start(frame, list(beta = c(9, 8)))
  expect_equal(unname(st$beta), c(9, 8))
  expect_error(make_start(frame, list(beta = 1)),
               "must have length 2")
  # and through frm(), as test-methods.R asserts
  expect_error(frm(bf(y ~ x) + gaussian(), dd, start = list(beta = 1)),
               "must have length")
})

test_that("a located prior places a nonlinear start", {
  set.seed(4)
  nd <- data.frame(t = rep(seq(0, 10, length.out = 25), 4))
  nd$y <- 8 * (1 - exp(-nd$t / 3)) + stats::rnorm(nrow(nd), 0, 0.3)
  nf <- bf(y ~ asym * (1 - exp(-t / lrc)), asym ~ 1, lrc ~ 1,
           nl = TRUE) + gaussian()
  pr <- prior(normal(5, 5), nlpar = "asym") +
    prior(normal(1, 5), nlpar = "lrc")

  cold <- par_template(nf, data = nd)
  expect_equal(unname(cold$beta), c(0, 0))

  placed <- par_template(nf, data = nd, prior = pr)
  expect_equal(unname(placed$beta), c(5, 1))

  # the fit that would die at zero now runs, and says so
  expect_message(fit <- frm(nf, data = nd, prior = pr),
                 "placed at the prior locations")
  expect_equal(unname(fixef(fit)$asym[["(Intercept)"]]), 8,
               tolerance = 0.1)

  # scoped to nonlinear parameters: a prior on an ordinary coefficient
  # moves no starting value
  dd <- pt_data()
  lin <- frm(bf(y ~ x) + gaussian(), dd, dry_run = "frame")
  ent <- resolve_prior_input(
    list(frame = lin, spec = lin$spec),
    set_prior("normal(100, 1)", class = "b"))$entries
  expect_equal(make_start(lin, NULL, ent), make_start(lin, NULL))
})

test_that("start= overrides a prior-placed nonlinear start silently", {
  set.seed(5)
  nd <- data.frame(t = rep(seq(0, 10, length.out = 25), 4))
  nd$y <- 8 * (1 - exp(-nd$t / 3)) + stats::rnorm(nrow(nd), 0, 0.3)
  nf <- bf(y ~ asym * (1 - exp(-t / lrc)), asym ~ 1, lrc ~ 1,
           nl = TRUE) + gaussian()
  pr <- prior(normal(5, 5), nlpar = "asym") +
    prior(normal(1, 5), nlpar = "lrc")
  st <- list(beta = c("asym_(Intercept)" = 6, "lrc_(Intercept)" = 2))

  frame <- frm(nf, data = nd, dry_run = "frame")
  ent <- resolve_prior_input(list(frame = frame, spec = frame$spec),
                             pr)$entries
  expect_equal(unname(make_start(frame, st, ent)$beta), c(6, 2))
  expect_no_message(frm(nf, data = nd, prior = pr, start = st))
})

test_that("a discovered template is accepted as newparams", {
  dd <- data.frame(x = stats::rnorm(40), g = factor(rep(1:4, 10)), y = 0)
  tp <- par_template(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  tp$beta[] <- c(1, 0.5)
  tp$betad[] <- log(0.7)
  tp$theta[] <- log(0.5)
  sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                       newparams = tp[c("beta", "betad", "theta")],
                       nsim = 2, seed = 1)
  expect_equal(ncol(sims), 2L)
  expect_true(all(is.finite(as.matrix(sims))))
})
