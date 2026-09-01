# Two open-tracker defects that share the "looks fine, is wrong" shape:
# prediction at a cell whose coefficient was aliased away (lme4#303) and
# grouping factors written as calls (lme4#464, #156, via reformulas).

# A rank-deficient design where h is a coarsening of f, so hB = fb + fc
# and one cell of the f x h table is unobservable.
aliased_data <- function() {
  set.seed(303)
  d <- data.frame(g = factor(rep(1:10, each = 6)))
  d$f <- factor(rep(c("a", "b", "c"), 20))
  d$h <- factor(ifelse(d$f == "a", "A", "B"))
  d$y <- rnorm(60) + as.integer(d$f)
  d
}

test_that("prediction at an aliased cell returns NA, not a partial sum (lme4#303)", {
  d <- aliased_data()
  expect_message(
    fit <- frm(bf(y ~ f + h + (1 | g)) + gaussian(), data = d),
    "dropping column\\(s\\): hB"
  )
  # the full-rank reparameterization of the same model
  ref <- frm(bf(y ~ f + (1 | g)) + gaussian(), data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_identical(fit$frame$linpreds[["y.mu"]]$dropped_colnames, "hB")

  lv_f <- levels(d$f)
  lv_h <- levels(d$h)
  lv_g <- levels(d$g)
  ok <- data.frame(f = factor(c("a", "b", "c"), levels = lv_f),
                   h = factor(c("A", "B", "B"), levels = lv_h),
                   g = factor(rep(1, 3), levels = lv_g))
  # estimable rows are exact, not merely close: same numbers as the
  # reparameterized fit, and no warning
  expect_silent(p_ok <- predict(fit, newdata = ok))
  expect_equal(unname(p_ok),
               unname(predict(ref, newdata = ok[, c("f", "g")])),
               tolerance = 1e-8)

  # (f = a, h = B) is the cell the data never saw
  mixed <- rbind(ok, data.frame(f = factor("a", levels = lv_f),
                                h = factor("B", levels = lv_h),
                                g = factor(1, levels = lv_g)))
  expect_warning(p_mix <- predict(fit, newdata = mixed),
                 "not estimable.*hB")
  expect_equal(unname(p_mix[1:3]), unname(p_ok), tolerance = 1e-12)
  expect_true(is.na(p_mix[4]))

  # one warning per call, however many rows are non-estimable
  two_bad <- mixed[c(4, 4), , drop = FALSE]
  w <- character(0)
  p_two <- withCallingHandlers(
    predict(fit, newdata = two_bad),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 1L)
  expect_true(all(is.na(p_two)))

  # se.fit follows the point prediction: NA there, unchanged elsewhere
  s_ok <- predict(fit, newdata = ok, se.fit = TRUE)
  s_mix <- suppressWarnings(predict(fit, newdata = mixed, se.fit = TRUE))
  expect_equal(s_mix$se.fit[1:3], s_ok$se.fit, tolerance = 1e-12)
  expect_true(is.na(s_mix$fit[4]))
  expect_true(is.na(s_mix$se.fit[4]))

  # in-sample paths use the fitted design, where the dropped columns
  # were consistently absent, so nothing changes there
  expect_silent(p_in <- predict(fit))
  expect_equal(unname(p_in), unname(predict(ref)), tolerance = 1e-8)
  expect_equal(unname(fitted(fit)), unname(fitted(ref)), tolerance = 1e-8)
  expect_silent(predict(fit, newdata = d))
  expect_false(anyNA(predict(fit, newdata = d)))
})

test_that("estimability tests the null space, not the dropped columns (lme4#303)", {
  set.seed(304)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$x2 <- 2 * dd$x                       # perfectly collinear
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  expect_message(m <- frm(bf(y ~ x + x2 + (1 | g)) + gaussian(), data = dd),
                 "rank deficient")
  m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  # x2 is nonzero on every row, yet every row restates the kept column,
  # so a "dropped column is nonzero" rule would blank the whole frame
  expect_silent(p <- predict(m, newdata = dd))
  expect_equal(p, predict(m0, newdata = dd), tolerance = 1e-6)
  # break the collinearity in newdata and the rows become non-estimable
  off <- transform(dd, x2 = x2 + 1)
  expect_warning(p_off <- predict(m, newdata = off), "not estimable")
  expect_true(all(is.na(p_off)))
})

test_that("grouping factors written as calls fit (lme4#464, #156)", {
  data(sleepstudy, package = "lme4")
  ss <- sleepstudy
  ss$xn <- as.integer(ss$Subject)
  ss$xf <- factor(ss$xn)
  ss$a <- factor(ss$Days %% 3)
  ss$b <- factor(ss$Days %% 2)
  ss$ab <- interaction(ss$a, ss$b)

  # reformulas re-evaluates the bar RHS inside the model frame, where
  # the call's arguments are not columns; the failure surfaced as an
  # error raised several frames down ("unique() applies only to vectors")
  f_call <- frm(bf(Reaction ~ Days + (1 | factor(xn))) + gaussian(),
                data = ss)
  f_col <- frm(bf(Reaction ~ Days + (1 | xf)) + gaussian(), data = ss)
  expect_equal(as.numeric(logLik(f_call)), as.numeric(logLik(f_col)),
               tolerance = 1e-10)
  expect_equal(unname(predict(f_call)), unname(predict(f_col)),
               tolerance = 1e-10)
  # the expression, not a synthetic column name, labels the term
  expect_identical(names(ranef(f_call)), "1 | factor(xn)")
  # prediction re-evaluates the expression against newdata
  expect_equal(unname(predict(f_call, newdata = ss)),
               unname(predict(f_col, newdata = ss)), tolerance = 1e-10)

  g_call <- frm(bf(Reaction ~ Days + (1 | interaction(a, b))) + gaussian(),
                data = ss)
  g_col <- frm(bf(Reaction ~ Days + (1 | ab)) + gaussian(), data = ss)
  expect_equal(as.numeric(logLik(g_call)), as.numeric(logLik(g_col)),
               tolerance = 1e-10)
  expect_equal(unname(predict(g_call, newdata = ss)),
               unname(predict(g_col, newdata = ss)), tolerance = 1e-10)

  # ':' and '/' groupings keep their reformulas expansion untouched
  h_colon <- frm(bf(Reaction ~ Days + (1 | Subject:a)) + gaussian(),
                 data = ss)
  expect_s3_class(h_colon, "frmtmb_fit")
  h_slash <- frm(bf(Reaction ~ Days + (1 | Subject/a)) + gaussian(),
                 data = ss)
  expect_length(ranef(h_slash), 2L)
  # a call nested inside ':' is resolved too
  h_mixed <- frm(bf(Reaction ~ Days + (1 | a:factor(b))) + gaussian(),
                 data = ss)
  expect_identical(h_mixed$frame$re_blocks[[1]]$levels,
                   levels(droplevels(ss$a:factor(ss$b))))
})
