# add_criterion() on draws, the pp_mixture() names, and the generics
# that moved to frmtmb (lane surface, 2026-10-06).
#
# Seen to fail on frmtmb.sample 0.16.0 (rellib-r6): add_criterion()
# did not exist ("could not find function"), and pp_mixture()'s summary
# named its components "class1", "class2" where brms's say
# "P(K = 1 | Y)".

test_that("stancode, standata and pp_mixture are frmtmb's generics", {
  for (g in c("stancode", "standata", "pp_mixture", "add_criterion")) {
    expect_true(g %in% getNamespaceExports("frmtmb.sample"), info = g)
    expect_false(g %in% names(frmtmb.sample:::sample_generic_owners),
                 info = g)
    # the binding a call reaches is frmtmb's (or, with brms loaded, the
    # owner's), never a rival defined here
    ns <- asNamespace("frmtmb.sample")
    expect_false(identical(environment(get(g, envir = ns)), ns), info = g)
  }
})

test_that("add_criterion() stores loo and waic, and loo() reads them", {
  skip_sampler()
  skip_if_not_installed("loo")
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x, 1)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x), family = gaussian(), data = dd, chains = 1,
               iter = 400, refresh = 0, seed = 3)))
  # a 60-row fit can raise loo's and waic's own diagnostics, which are
  # theirs and not this method's
  allow_warnings({
  ds2 <- add_criterion(ds, c("loo", "waic"), model_name = "m1")
  expect_s3_class(ds2, "frmtmb_draws")
  expect_setequal(names(ds2$criteria), c("loo", "waic"))
  expect_s3_class(ds2$criteria$loo, "psis_loo")
  expect_identical(attr(ds2$criteria$loo, "model_name"), "m1")
  # brms's use_stored: with nothing else asked, the stored object
  expect_identical(loo(ds2), ds2$criteria$loo)
  expect_identical(waic(ds2), ds2$criteria$waic)
  # a further argument recomputes, as brms's does
  expect_false(identical(attr(loo(ds2, ndraws = 100), "model_name"), "m1"))
  # the stored value is what loo() computes
  expect_equal(ds2$criteria$loo$estimates, loo(ds)$estimates)
  # kept unless overwrite = TRUE
  ds3 <- add_criterion(ds2, "loo", model_name = "m2")
  expect_identical(attr(ds3$criteria$loo, "model_name"), "m1")
  ds4 <- add_criterion(ds2, "loo", model_name = "m2", overwrite = TRUE)
  expect_identical(attr(ds4$criteria$loo, "model_name"), "m2")
  # overwrite = TRUE RECOMPUTES, as brms's clears the stored criterion
  # before it computes: a loo stored from 100 draws is replaced by one
  # from all of them (the review of 2026-10-07, B1: the stored object
  # came back, because loo() returned it)
  dsa <- add_criterion(ds, "loo", ndraws = 100)
  expect_equal(attr(dsa$criteria$loo, "dims")[1L], 100L)
  dsb <- add_criterion(dsa, "loo", overwrite = TRUE)
  expect_equal(attr(dsb$criteria$loo, "dims")[1L], ndraws(ds))
  expect_equal(dsb$criteria$loo$estimates, loo(ds)$estimates)
  # bayes_R2 stored as brms stores it, the draws
  ds5 <- add_criterion(ds, "bayes_R2")
  expect_equal(length(ds5$criteria$bayes_R2), ndraws(ds))
  # the refusals of the criteria frmtmb.sample does not compute
  expect_error(add_criterion(ds, "kfold"), "kfold() is not implemented",
               fixed = TRUE)
  expect_error(add_criterion(ds, "marglik"), "bridge_sampler()",
               fixed = TRUE)
  expect_error(add_criterion(ds, "loo_R2"), "loo_R2", fixed = TRUE)
  expect_error(add_criterion(ds, "aic"), "should be a subset of")
  }, c("Pareto k", "p_waic"))
})

test_that("add_criterion() on a fit refuses and names the draws route", {
  set.seed(1)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- dd$x + stats::rnorm(40)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_error(add_criterion(fit, "loo"),
               "add_criterion(frmtmb.sample::frm_sample(fit), \"loo\")",
               fixed = TRUE)
})

test_that("pp_mixture() on draws names its components as brms does", {
  skip_sampler()
  set.seed(4)
  dd <- data.frame(y = c(stats::rnorm(50, -2), stats::rnorm(50, 3)))
  fit <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()),
             data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 2)))
  st <- pp_mixture(ds)
  expect_equal(dimnames(st)[[3L]], c("P(K = 1 | Y)", "P(K = 2 | Y)"))
  expect_equal(dimnames(st)[[1L]], as.character(seq_len(100L)))
  # the fit's method has the same names
  expect_identical(dimnames(frmtmb::pp_mixture(fit))[2:3],
                   dimnames(st)[2:3])
})

test_that("plot() of draws is brms's display, with brms's arguments", {
  # Seen to fail on frmtmb.sample 0.16.0: plot() of draws refused every
  # call and named mcmc_plot().
  skip_sampler()
  skip_if_not_installed("bayesplot")
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x, 1)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + z), family = gaussian(), data = dd, chains = 2,
               iter = 300, refresh = 0, seed = 3)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # the outer parameters, five to a page by default
  p <- plot(ds, ask = FALSE)
  expect_length(p, ceiling(length(frmtmb.sample:::draws_outer_cols(ds)) / 5))
  # brms's variable and regex, and nvariables per page
  p2 <- plot(ds, variable = "^b_", regex = TRUE, nvariables = 2,
             ask = FALSE)
  expect_length(p2, 2L)
  # brms's deprecated N, honored with brms's warning
  p3 <- allow_warnings(plot(ds, N = 1, variable = "b_x", plot = FALSE),
                       "Argument 'N' is deprecated",
                       require = "Argument 'N' is deprecated")
  expect_length(p3, 1L)
  expect_error(plot(ds, variable = "^nosuch", regex = TRUE),
               "No valid variables selected.", fixed = TRUE)
  expect_error(plot(ds, nvariables = 0),
               "must be a positive integer", fixed = TRUE)
})

test_that("frm_sample() gives a frame-build warning once, as frm() does", {
  # The 0.68.0 consolidation filed this as "6 times"; the count of 6 was
  # every warning of the call, the sampler's five beside this one
  # (dev/surface-findings.md, item 5). Kept as a guard.
  skip_sampler()
  set.seed(20261006)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  u <- stats::rlogis(n) + 0.5 * d$x + 0.4 * d$z
  d$y <- 1L + (u > -1 + 0.15 * d$x) + (u > 0.3) + (u > 1.5 - 0.15 * d$x)
  w <- character()
  suppressMessages(withCallingHandlers(
    frm_sample(y ~ z + cs(x), family = cumulative(), data = d, chains = 1,
               iter = 300, refresh = 0, seed = 3,
               prior = set_prior("normal(0, 2)", class = "b")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }))
  expect_equal(sum(grepl("Category specific effects", w, fixed = TRUE)), 1L)
})
