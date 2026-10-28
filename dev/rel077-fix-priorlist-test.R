test_that("a list of specifications points to + and c()", {
  # the 0.77.0 consolidation: list(sp, sp2) got the generic "not a list
  # of length 2"; brms combines specifications with + or c()
  d <- ros_data()
  sp <- set_prior("normal(0, 1)", class = "b")
  sp2 <- set_prior("normal(0, 2)", class = "Intercept")
  expect_error(frm(bf(y ~ x), data = d, prior = list(sp, sp2)),
               "Combine them with + or c(), as in brms: sp1 + sp2",
               fixed = TRUE)
  # the guard absent: both spellings fit, and give the same fit
  f1 <- frm(bf(y ~ x), data = d, prior = sp + sp2)
  f2 <- frm(bf(y ~ x), data = d, prior = c(sp, sp2))
  expect_identical(f1$opt$par, f2$opt$par)
})
