test_that("c() reads a prior table as brms's c() does, keeping its route", {
  # the 0.77.0 release review, finding 3: c(set_prior(), default_prior())
  # was refused as "a plain data frame"; 0.76.0 combined them, and brms
  # 2.23.0's c() combines a brmsprior and a prior table in either order
  # (dev/rel077-fix2-brms-c.R)
  d <- ros_data()
  f <- y ~ x + (1 + x | g)
  pf <- default_prior(f, data = d)
  sp <- set_prior("normal(0, 1)", class = "b", coef = "x")
  for (pl in list(c(sp, pf), c(pf, sp))) {
    expect_s3_class(pl, "frmtmb_priorlist")
    expect_identical(attr(pl, "route"), "fit")
    expect_length(pl, length(as.brmsprior(pf)) + 1L)
  }
  # the rows are the specification and the table read as as.brmsprior()
  # reads it
  expect_equal(as.data.frame(c(sp, pf))[-1L, ],
               as.data.frame(as.brmsprior(pf)), ignore_attr = TRUE)
  # the guard absent: what is not a specification or a table is refused
  expect_error(c(sp, prior_normal(0, 1)),
               "a density object alone names no parameter", fixed = TRUE)
})
