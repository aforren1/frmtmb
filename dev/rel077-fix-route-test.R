test_that("as.data.frame() drops the route, and the (flat) rows are refused", {
  # the 0.77.0 consolidation: as.data.frame() of a default_prior() table
  # drops its route, and through 0.76.0 plus the lanes,
  # as.brmsprior(as.data.frame(pe)) then put brms's (flat) on the
  # intercept, sigma, sd and cor of y ~ x + (1 + x | g) in silence. Its
  # rows carry frmtmb's printed "(flat)" with source "default", so it is
  # marked, and refused where such a row replaces a proper default
  d <- fp_data()
  f <- y ~ x + (1 + x | g)
  fit <- frm(f, data = d, dry_run = "objective")
  p <- default_prior(f, data = d)
  p$prior[p$class == "b" & p$coef == "x"] <- "normal(0, 1)"
  pl <- as.brmsprior(as.data.frame(p))
  expect_identical(attr(pl, "route"), "table")
  sp <- set_prior("normal(0, 1)", class = "b")
  msg <- "that as.data.frame() turned into a plain data frame"
  for (pr in list(pl, sp + pl, c(pl, sp))) {
    err <- tryCatch(fp_rp(fit, pr), error = conditionMessage)
    expect_true(is.character(err) && grepl(msg, err, fixed = TRUE))
    # it names what the (flat) rows would have replaced
    expect_true(all(vapply(c("Intercept", "sigma", "sd", "cor"), grepl,
                           NA, x = err, fixed = TRUE)))
  }
  expect_error(suppressMessages(frm_sample(fit, prior = pl)), msg,
               fixed = TRUE)
  # the guard absent: the sample route's table, coerced and edited, is
  # brms's defaults with the edited row, as passed whole
  ps <- default_prior(f, data = d, route = "sample")
  ps$prior[ps$class == "b" & ps$coef == "x"] <- "normal(0, 1)"
  pls <- as.brmsprior(as.data.frame(ps))
  key <- function(es) sort(vapply(es, function(e) {
    paste(e$comp, paste(e$idx, collapse = ","), e$dist$kind)
  }, ""))
  expect_identical(key(fp_rp(fit, pls)$ri$entries),
                   key(fp_rp(fit, ps)$ri$entries))
  # and the plain-data-frame refusal says to pass the table as it is
  expect_error(fp_rp(fit, as.data.frame(p)),
               "A default_prior() table is passed as it is", fixed = TRUE)
})
