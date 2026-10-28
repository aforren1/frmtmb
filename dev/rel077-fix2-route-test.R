test_that("subset() and a column selection keep the route, and are refused", {
  # the 0.77.0 release review, finding 1 (seed 11, y ~ x + (1 + x | g)):
  # `[.data.frame` dropped the route on a column index, so subset(pf, ...)
  # and pf[, cols] of frm()'s table were read as written, and (flat)
  # replaced brms's defaults on the intercept, sigma and the sd in silence
  set.seed(11)
  d <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6)[d$g] + stats::rnorm(60)
  f <- y ~ x + (1 + x | g)
  fit <- frm(f, data = d, dry_run = "objective")
  pf <- default_prior(f, data = d)
  pf$prior[pf$class == "b" & pf$coef == "x"] <- "normal(0, 1)"
  msg <- "default_prior() or validate_prior() table of route = \"fit\""
  cols <- c("prior", "class", "coef", "group")
  for (pr in list(subset(pf, class != "cor"), pf[, cols], pf[cols])) {
    expect_identical(attr(pr, "route"), "fit")
    expect_error(fp_rp(fit, pr), msg, fixed = TRUE)
  }
  expect_error(suppressMessages(frm_sample(fit, prior = subset(pf,
                                                             class != "cor"))),
               msg, fixed = TRUE)
  # the review's two smaller gaps: the coerced table with its `source`
  # column dropped, or set to "user", is marked by its printed (flat)
  cm <- "that as.data.frame() turned into a plain data frame"
  pd <- as.data.frame(pf)
  pu <- pd
  pu$source <- "user"
  for (pr in list(pd[, setdiff(names(pd), "source")], pu)) {
    expect_error(fp_rp(fit, as.brmsprior(pr)), cm, fixed = TRUE)
  }
  # the guard absent: subset() of the sample route's table is brms's
  # defaults, as the table passed whole is
  key <- function(es) sort(vapply(es, function(e) {
    paste(e$comp, paste(e$idx, collapse = ","), e$dist$kind)
  }, ""))
  ps <- default_prior(f, data = d, route = "sample")
  sub <- subset(ps, class != "cor")
  expect_identical(attr(sub, "route"), "sample")
  expect_identical(key(fp_rp(fit, sub)$ri$entries),
                   key(fp_rp(fit, NULL)$ri$entries))
})
