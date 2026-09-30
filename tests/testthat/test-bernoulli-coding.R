# bernoulli() on a two-valued response that is not 0/1. brms's
# data_response() codes the two values 0 and 1 by level order,
# levels(as.factor(y)), and stores the coding with the fit
# (frame$basis$resp_levels) so that newdata is coded the same way.
# Before, frmtmb refused anything but 0/1 ("bernoulli: response must be
# 0/1", ledger row standata:75). brms's codes for each input type are in
# dev/formrobust-log/brms-binary.txt.

bc_data <- local({
  set.seed(5)
  n <- 120
  d <- data.frame(x = stats::rnorm(n))
  yy <- stats::rbinom(n, 1, stats::plogis(-0.3 + 1.2 * d$x))
  d$y01 <- yy
  d$ym <- ifelse(yy == 1, -1, -2)         # -2 codes 0, -1 codes 1
  d$yf <- factor(ifelse(yy == 1, "yes", "no"))
  d$yl <- yy == 1
  d$y12 <- yy + 1
  d$yc <- ifelse(yy == 1, "yes", "no")
  d
})

bc_fit <- function(v) {
  frm(bf(stats::as.formula(paste(v, "~ x"))), data = bc_data,
      family = bernoulli())
}

test_that("every two-valued type fits the 0/1 model, bitwise", {
  ref <- bc_fit("y01")
  for (v in c("ym", "yf", "yl", "y12", "yc")) {
    f <- bc_fit(v)
    # an identity: the same codes reach the same tape
    expect_identical(logLik(f), logLik(ref), info = v)
    expect_identical(fixef(f), fixef(ref), info = v)
  }
})

test_that("the codes are brms's level order", {
  cd <- function(y) {
    fr <- frm(bf(y ~ 1), data = data.frame(y = y), family = bernoulli(),
              dry_run = "frame")
    fr$y[[1L]]
  }
  expect_identical(cd(rep(-c(1, 2), 5)), rep(c(1, 0), 5))
  expect_identical(cd(c(1, 2, 2, 1)), c(0, 1, 1, 0))
  expect_identical(cd(c(5, 5, 5)), c(1, 1, 1))
  expect_identical(cd(c(0, 0, 0)), c(0, 0, 0))
  expect_identical(cd(factor(c("a", "z", "a"), levels = c("z", "a"))),
                   c(1, 0, 1))
  # brms codes an all-TRUE logical 0, since as.factor(TRUE) has one
  # level; frmtmb reads a logical as its number (a recorded divergence)
  expect_identical(cd(c(TRUE, TRUE)), c(1, 1))
  expect_error(cd(c(0, 1, 2)), "requires responses to contain only two ",
               fixed = TRUE)
  # two fractions are coded as brms codes them, with a warning, since
  # they are usually a proportion (lme4#682)
  pw <- "lie strictly between 0 and 1"
  allow_warnings(expect_identical(cd(c(0.1, 0.9, 0.1)), c(0, 1, 0)), pw,
                 require = pw)
  # the guard-absent cases: effect coding, and a value at 0 or 1, are
  # not a proportion and do not warn
  expect_no_warning(expect_identical(cd(c(-0.5, 0.5, 0.5)), c(0, 1, 1)))
  expect_no_warning(expect_identical(cd(c(0, 0.5, 0)), c(0, 1, 0)))
})

test_that("brms's standata agrees", {
  skip_unless_brms()
  for (y in list(rep(-c(1, 2), 5), c(1, 2, 2, 1), c(5, 5, 5),
                 c(-0.5, 0.5, 0.5), c(0.1, 0.9, 0.9),
                 factor(c("b", "a", "b")), c("yes", "no", "yes"),
                 c(TRUE, FALSE, TRUE))) {
    d <- data.frame(y = y)
    sd <- brms_standata(y ~ 1, data = d, family = brms::bernoulli())
    fr <- NULL
    allow_warnings(fr <- frm(bf(y ~ 1), data = d, family = bernoulli(),
                             dry_run = "frame"),
                   "lie strictly between 0 and 1")
    expect_identical(as.numeric(sd$Y), fr$y[[1L]])
  }
})

test_that("newdata is coded as the fit coded its response", {
  fm <- bc_fit("ym")
  ref <- bc_fit("y01")
  # a newdata holding only -2: coded afresh it would read as the 1 of
  # c(0, -2), so this is the case the stored coding exists for
  nd <- bc_data[bc_data$ym == -2, ][1:6, ]
  expect_identical(residuals(fm, newdata = nd),
                   residuals(ref, newdata = nd))
  nd1 <- bc_data[bc_data$ym == -1, ][1:6, ]
  expect_identical(residuals(fm, newdata = nd1),
                   residuals(ref, newdata = nd1))
  # pp_check(newdata = ) plots newdata's response on the draws' codes
  expect_identical(pp_check_newdata_y(fm, fm$spec$responses[[1L]], nd),
                   rep(0, 6))
  bad <- nd
  bad$ym[1] <- 3
  expect_error(residuals(fm, newdata = bad), "'3' are neither", fixed = TRUE)
})

test_that("post-fit methods agree with the 0/1 fit", {
  fm <- bc_fit("ym")
  ref <- bc_fit("y01")
  expect_identical(fitted(fm), fitted(ref))
  expect_identical(residuals(fm), residuals(ref))
  expect_identical(simulate(fm, nsim = 2, seed = 1)[[1L]],
                   simulate(ref, nsim = 2, seed = 1)[[1L]])
  set.seed(4)
  p1 <- predict(fm)
  set.seed(4)
  p0 <- predict(ref)
  expect_identical(p1, p0)
  expect_identical(AIC(fm), AIC(ref))
  i1 <- influence(fm)
  i0 <- influence(ref)
  expect_identical(i1$fixef, i0$fixef)
  b1 <- frm_bootstrap(fm, FUN = function(f) fixef(f)[, 1], nsim = 3,
                      seed = 2)
  b0 <- frm_bootstrap(ref, FUN = function(f) fixef(f)[, 1], nsim = 3,
                      seed = 2)
  expect_identical(b1$t, b0$t)
  ce1 <- conditional_effects(fm, effects = "x")[[1L]]
  ce0 <- conditional_effects(ref, effects = "x")[[1L]]
  expect_identical(ce1$estimate__, ce0$estimate__)
  # refit() takes the 0/1 codes simulate() returns; the response's own
  # values are refused by name, where they fitted a logLik of 3e304
  expect_error(refit(fm, replace(bc_data$y01, 1, NA)), "1 NA value",
               fixed = TRUE)
  expect_error(refit(fm, bc_data$ym), "refitted on its 0/1 codes",
               fixed = TRUE)
  expect_identical(logLik(refit(fm, bc_data$y01)),
                   logLik(refit(ref, bc_data$y01)))
  # the coding travels with the fit's spec, which a refit reads
  expect_identical(fm$spec$responses[[1L]]$family[["bin_levels"]],
                   c("-2", "-1"))
})
