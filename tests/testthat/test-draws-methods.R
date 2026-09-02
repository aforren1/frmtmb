# The brmsfit-shaped method surface on frm_sample() draws: shape
# conversions, posterior summaries, structural delegations, sampler
# diagnostics, and the refusals that replace "could not find function"
# for a ported brms script.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

skip_sampler <- function() {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
}

dm_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(9)
      dd <- data.frame(x = stats::rnorm(60),
                       g = factor(rep(1:6, 10)))
      dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x +
                             stats::rnorm(6, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 400, refresh = 0, seed = 1)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

fake_draws <- function(fit, n = 4L) {
  lab <- c(frmtmb:::all_par_labels(fit), "lp__")
  structure(list(stanfit = NULL,
                 draws = matrix(0, n, length(lab),
                                dimnames = list(NULL, lab)),
                 fit = fit),
            class = "frmtmb_draws")
}

## ---- shapes ---------------------------------------------------------

test_that("the dimension accessors describe the draws matrix", {
  cs <- dm_case()
  expect_equal(ndraws(cs$ds), nrow(cs$ds$draws))
  expect_equal(nchains(cs$ds), 2L)
  expect_equal(niterations(cs$ds), ndraws(cs$ds) / 2L)
  expect_equal(nvariables(cs$ds), ncol(cs$ds$draws))
  expect_equal(nvariables(cs$ds), length(variables(cs$ds)))
})

test_that("as.array() keeps the chains apart, in draws order", {
  cs <- dm_case()
  a <- as.array(cs$ds)
  expect_equal(dim(a), c(niterations(cs$ds), 2L, nvariables(cs$ds)))
  expect_equal(dimnames(a)[[3L]], colnames(cs$ds$draws))
  # frm_sample() rbinds the chains in order, so chain 2 is the second
  # half of the matrix; a wrong reshape here would silently make every
  # convergence diagnostic meaningless
  k <- niterations(cs$ds)
  expect_equal(as.numeric(a[, 1L, "x"]), unname(cs$ds$draws[seq_len(k), "x"]))
  expect_equal(as.numeric(a[, 2L, "x"]),
               unname(cs$ds$draws[k + seq_len(k), "x"]))
})

test_that("the posterior converters all round-trip the same draws", {
  cs <- dm_case()
  skip_if_not_installed("posterior")
  expect_s3_class(as_draws(cs$ds), "draws_matrix")
  expect_s3_class(as_draws_matrix(cs$ds), "draws_matrix")

  arr <- as_draws_array(cs$ds)
  expect_s3_class(arr, "draws_array")
  expect_equal(posterior::nchains(arr), 2L)
  expect_equal(posterior::ndraws(arr), ndraws(cs$ds))

  df <- as_draws_df(cs$ds)
  expect_s3_class(df, "draws_df")
  expect_equal(nrow(df), ndraws(cs$ds))
  expect_equal(df$x, unname(as.numeric(as_draws_matrix(cs$ds)[, "x"])))

  expect_s3_class(as_draws_list(cs$ds), "draws_list")
  rv <- as_draws_rvars(cs$ds)
  expect_s3_class(rv, "draws_rvars")
  expect_true("x" %in% names(rv))

  expect_equal(as.matrix(cs$ds), cs$ds$draws)
  expect_equal(as.data.frame(cs$ds)$x, unname(cs$ds$draws[, "x"]))
})

test_that("as.mcmc() gives coda one component per chain", {
  cs <- dm_case()
  skip_if_not_installed("coda")
  m <- as.mcmc(cs$ds)
  expect_s3_class(m, "mcmc.list")
  expect_length(m, 2L)
  expect_equal(colnames(m[[1L]]), colnames(cs$ds$draws))
  expect_equal(nrow(m[[1L]]), niterations(cs$ds))
  expect_true(is.finite(coda::gelman.diag(m[, "x"])$psrf[1L, 1L]))

  one <- as.mcmc(cs$ds, combine_chains = TRUE)
  expect_s3_class(one, "mcmc")
  expect_equal(nrow(one), ndraws(cs$ds))
})

## ---- summaries and intervals ----------------------------------------

test_that("posterior_summary() and posterior_interval() summarize draws", {
  cs <- dm_case()
  s <- posterior_summary(cs$ds)
  expect_equal(colnames(s), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  # the group-level modes are not what posterior_summary(ds) is asking
  # for, the same columns summary() and print() leave out
  expect_false(any(grepl("^b\\[", rownames(s))))
  expect_false("lp__" %in% rownames(s))
  expect_equal(unname(s["x", "Estimate"]), mean(cs$ds$draws[, "x"]))
  expect_equal(unname(s["x", "Q2.5"]),
               unname(stats::quantile(cs$ds$draws[, "x"], 0.025)))

  # it works on any matrix of draws, which is what makes the brms idiom
  # posterior_summary(bayes_R2(ds, summary = FALSE)) work
  R2 <- bayes_R2(cs$ds, summary = FALSE)
  expect_equal(posterior_summary(R2)[1, "Estimate"], mean(R2))

  rb <- posterior_summary(cs$ds, robust = TRUE, variable = "x")
  expect_equal(unname(rb[1, "Estimate"]),
               unname(stats::median(cs$ds$draws[, "x"])))

  pi <- posterior_interval(cs$ds, prob = 0.9)
  expect_equal(colnames(pi), c("5%", "95%"))
  expect_true(all(pi[, 1] < pi[, 2]))
  expect_error(posterior_interval(cs$ds, variable = "nope"),
               "which the draws do not contain")
})

test_that("predictive_interval() and predictive_error() use the predictive draws", {
  cs <- dm_case()
  set.seed(1)
  pit <- predictive_interval(cs$ds, prob = 0.8, ndraws = 100)
  expect_equal(dim(pit), c(nrow(cs$dd), 2L))
  expect_equal(colnames(pit), c("10%", "90%"))
  # predictive intervals cover the observed data: they carry the
  # family's own noise on top of the parameter uncertainty
  expect_gt(mean(cs$dd$y > pit[, 1] & cs$dd$y < pit[, 2]), 0.6)

  pe <- predictive_error(cs$ds, ndraws = 20)
  expect_equal(dim(pe), c(20L, nrow(cs$dd)))
  # brms's sign convention: y - yrep, so the errors centre on zero
  expect_lt(abs(mean(pe)), 0.5)
})

## ---- structural delegation ------------------------------------------

test_that("nobs/formula/family/getCall/ngrps report the sampled model", {
  cs <- dm_case()
  expect_equal(stats::nobs(cs$ds), nrow(cs$dd))
  expect_equal(stats::formula(cs$ds), stats::formula(cs$fit))
  expect_equal(stats::family(cs$ds)$family, "gaussian")
  expect_equal(ngrps(cs$ds), ngrps(cs$fit))
  expect_true(is.call(getCall(cs$ds)))
})

test_that("the structural accessors work on formula-route draws", {
  skip_sampler()
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(50), g = factor(rep(1:5, 10)))
  dd$y <- stats::rnorm(50, 1 + 0.5 * dd$x, 1)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
               chains = 1, iter = 200, refresh = 0, seed = 2)))
  # these read structure only, so the "no maximum-likelihood estimate"
  # refusal must not fire on them
  expect_equal(stats::nobs(ds), 50L)
  expect_equal(ngrps(ds), c(g = 5L))
  expect_equal(stats::family(ds)$family, "gaussian")
})

test_that("coef() is fixef broadcast plus each group's own draws", {
  cs <- dm_case()
  cf <- coef(cs$ds)
  expect_named(cf, "g")
  expect_equal(dim(cf$g), c(6L, 4L, 2L))
  expect_equal(dimnames(cf$g)[[2L]],
               c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(dimnames(cf$g)[[3L]], c("(Intercept)", "x"))

  # the estimate is the posterior mean of (fixed + that group's random
  # intercept), which is the fit-side coef() computed per draw
  idx <- frmtmb:::draws_par_index(cs$ds$fit)
  per <- vapply(seq_len(ndraws(cs$ds)), function(i) {
    coef(frmtmb:::draws_fit_at(cs$ds, i, idx))$g[["(Intercept)"]]
  }, numeric(6L))
  expect_equal(unname(cf$g[, "Estimate", "(Intercept)"]),
               unname(rowMeans(per)), tolerance = 1e-12)

  # a slope with no group-level term is the same in every group
  expect_equal(length(unique(round(cf$g[, "Estimate", "x"], 12))), 1L)
})

## ---- sampler diagnostics and plots ----------------------------------

test_that("the bayesplot accessors read the stanfit", {
  cs <- dm_case()
  skip_if_not_installed("bayesplot")
  expect_equal(unname(rhat(cs$ds)),
               unname(bayesplot::rhat(cs$ds$stanfit)))
  expect_equal(unname(neff_ratio(cs$ds)),
               unname(bayesplot::neff_ratio(cs$ds$stanfit)))
  np <- nuts_params(cs$ds)
  expect_true(all(c("Chain", "Iteration", "Parameter", "Value") %in%
                    names(np)))
  expect_true("divergent__" %in% levels(np$Parameter))
  lp <- log_posterior(cs$ds)
  expect_equal(nrow(lp), ndraws(cs$ds))
})

test_that("mcmc_plot() and pairs() call bayesplot on the draws array", {
  cs <- dm_case()
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("ggplot2")
  expect_s3_class(mcmc_plot(cs$ds), "ggplot")
  expect_s3_class(mcmc_plot(cs$ds, type = "trace", variable = "x"),
                  "ggplot")
  expect_s3_class(mcmc_plot(cs$ds, type = "hist"), "ggplot")
  expect_error(mcmc_plot(cs$ds, type = "not_a_plot"),
               "which does not exist")
  expect_s3_class(pairs(cs$ds, variable = c("Intercept", "x")), "bayesplot_grid")
})

## ---- mixture membership ---------------------------------------------

test_that("pp_mixture() propagates parameter uncertainty into the probabilities", {
  skip_sampler()
  set.seed(4)
  dd <- data.frame(y = c(stats::rnorm(50, -2), stats::rnorm(50, 3)))
  fit <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()),
             data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 2)))

  raw <- pp_mixture(ds, summary = FALSE)
  expect_equal(dim(raw), c(ndraws(ds), 100L, 2L))
  expect_equal(unname(apply(raw, c(1, 2), sum)),
               matrix(1, ndraws(ds), 100L), tolerance = 1e-10)

  # each slice is the fit-side computation at that draw
  idx <- frmtmb:::draws_par_index(ds$fit)
  expect_equal(unname(raw[3L, , ]),
               unname(mixture_probs(frmtmb:::draws_fit_at(ds, 3L, idx))),
               tolerance = 1e-12)

  st <- pp_mixture(ds)
  expect_equal(dim(st), c(100L, 4L, 2L))
  expect_equal(unname(st[, "Estimate", 1L]),
               unname(apply(raw[, , 1L], 2L, mean)), tolerance = 1e-12)
  # the data are well separated, so the assignment is nearly certain
  expect_gt(mean(apply(st[, "Estimate", ], 1, max)), 0.95)
})

## ---- refusals and renamed spellings ---------------------------------

test_that("the brms-only methods refuse with a reason and a replacement", {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(40), g = factor(rep(1:4, 10)))
  dd$y <- stats::rnorm(40)
  uf <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            dry_run = "objective")
  fd <- fake_draws(uf)
  expect_error(stancode(fd), "no Stan program")
  expect_error(standata(fd), "Stan data list")
  expect_error(expose_functions(fd), "already plain R")
  expect_error(restructure(fd), "upgrade path")
  expect_error(posterior_samples(fd), "as_draws\\(x\\)")
  expect_error(nsamples(fd), "ndraws\\(x\\)")
  expect_error(parnames(fd), "variables\\(x\\)")
})

test_that("the matrix-response guards name the function that hit them", {
  skip_sampler()
  set.seed(6)
  K <- 3L
  P <- matrix(c(0.5, 0.3, 0.2), nrow = 1)
  Y <- t(vapply(seq_len(40), function(i) {
    stats::rmultinom(1, 10, P)[, 1]
  }, numeric(K)))
  colnames(Y) <- c("a", "b", "c")
  dd <- data.frame(x = stats::rnorm(40))
  dd$Y <- Y
  fit <- frm(bf(Y | trials(10) ~ 1), family = multinomial(K = 3),
             data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 2)))
  expect_error(predictive_interval(ds, ndraws = 5),
               "one predicted number per")
  expect_error(predictive_error(ds, ndraws = 5), "vector response")
})

## ---- conditional_effects --------------------------------------------

test_that("conditional_effects() bands the drawn curves", {
  cs <- dm_case()
  ce <- conditional_effects(cs$ds, effects = "x", resolution = 25)
  expect_s3_class(ce, "frmtmb_conditional_effects")
  df <- ce$x
  expect_equal(nrow(df), 25L)
  expect_true(all(is.finite(df$estimate__)))
  expect_true(all(df$lower__ <= df$estimate__ &
                    df$estimate__ <= df$upper__))
  expect_identical(attr(df, "band"), "posterior")
  # same density (flat priors), so the posterior-mean curve tracks the
  # maximum-likelihood curve; the yardstick is the wider of the two
  # bands' own standard errors, so a platform whose chain drifted still
  # judges wiring, not mixing
  cf <- conditional_effects(cs$fit, effects = "x", resolution = 25)
  expect_lt(max(abs(df$estimate__ - cf$x$estimate__) /
                  pmax(df$se__, cf$x$se__, 1e-8)), 5)
  # thinning changes the cost, not the shape
  ce5 <- conditional_effects(cs$ds, effects = "x", resolution = 25,
                             ndraws = 25)
  expect_equal(nrow(ce5$x), 25L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(ce, ask = FALSE, points = TRUE))
})

test_that("conditional_effects() runs on formula-route draws", {
  skip_sampler()
  set.seed(11)
  dd <- data.frame(x = stats::rnorm(50), g = factor(rep(1:5, 10)))
  dd$y <- stats::rnorm(50, 1 + 0.5 * dd$x, 1)
  dsf <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
               chains = 1, iter = 300, refresh = 0, seed = 3)))
  ce <- conditional_effects(dsf, effects = "x", resolution = 10)
  expect_s3_class(ce, "frmtmb_conditional_effects")
  expect_true(all(is.finite(ce$x$estimate__)))
  # the embedded fit alone cannot draw this curve: no ML estimates
  expect_error(conditional_effects(dsf$fit), "needs a fitted model")
})

test_that("conditional_effects() on draws refuses what it cannot mean", {
  cs <- dm_case()
  expect_error(conditional_effects(cs$ds, method = "predict"),
               "no method =")
  expect_error(conditional_effects(cs$ds, band = "boot"),
               "no band =")
  # laplace-shaped draws: random effects in the model, no b[ columns
  ld <- cs$ds
  ld$draws <- ld$draws[, !startsWith(colnames(ld$draws), "b["),
                       drop = FALSE]
  expect_error(conditional_effects(ld), "laplace = TRUE")
})
