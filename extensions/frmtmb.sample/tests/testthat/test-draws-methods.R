# The brmsfit-shaped method surface on frm_sample() draws: shape
# conversions, posterior summaries, structural delegations, sampler
# diagnostics, and the refusals that replace "could not find function"
# for a ported brms script.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

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
  lab <- c(frmtmb.sample:::all_par_labels(fit), "lp__")
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
  idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
  per <- vapply(seq_len(ndraws(cs$ds)), function(i) {
    coef(frmtmb.sample:::draws_fit_at(cs$ds, i, idx))$g[["(Intercept)"]]
  }, numeric(6L))
  expect_equal(unname(cf$g[, "Estimate", "(Intercept)"]),
               unname(rowMeans(per)), tolerance = 1e-12)

  # a slope with no group-level term is the same in every group
  expect_equal(length(unique(round(cf$g[, "Estimate", "x"], 12))), 1L)
})

## ---- sampler diagnostics and plots ----------------------------------

test_that("rhat() and neff_ratio() are brms's, on this package's names", {
  # brms's rhat.brmsfit is summarise_draws(rhat = posterior::rhat), the
  # rank-normalized split-R-hat, and its neff_ratio.brmsfit is
  # min(ess_bulk, ess_tail) / ndraws. Both used to be bayesplot on the
  # stanfit, which is rstan's classic split-R-hat and rstan's n_eff and
  # which reported STAN's parameter names. dev/brmsmatch-findings.md
  # has the size of both gaps, measured on 4 chains of 500.
  cs <- dm_case()
  skip_if_not_installed("posterior")
  a <- posterior::as_draws_array(cs$ds)
  want_r <- posterior::summarise_draws(a, rhat = posterior::rhat)
  expect_equal(rhat(cs$ds),
               stats::setNames(want_r$rhat, want_r$variable))
  want_e <- posterior::summarise_draws(a, ess_bulk = posterior::ess_bulk,
                                       ess_tail = posterior::ess_tail)
  expect_equal(neff_ratio(cs$ds),
               stats::setNames(pmin(want_e$ess_bulk, want_e$ess_tail) /
                                 posterior::ndraws(a), want_e$variable))

  # the frmtmb names fall out of computing through posterior on the
  # relabeled draws, so these two line up with every other accessor
  expect_equal(names(rhat(cs$ds)), variables(cs$ds))
  expect_equal(names(neff_ratio(cs$ds)), variables(cs$ds))
  expect_false(is.na(rhat(cs$ds)["x"]))
  expect_false(is.na(neff_ratio(cs$ds)["x"]))

  # both move away from the sampler's own numbers, and summary() moves
  # with them rather than being left behind
  skip_if_not_installed("bayesplot")
  expect_false(isTRUE(all.equal(unname(rhat(cs$ds)),
                                unname(bayesplot::rhat(cs$ds$stanfit)))))
  expect_false(isTRUE(all.equal(
    unname(neff_ratio(cs$ds)),
    unname(bayesplot::neff_ratio(cs$ds$stanfit)))))
})

test_that("summary() reports brms's three diagnostics, agreeing with rhat()", {
  # brms:::summary.brmsfit adds
  #   Rhat = posterior::rhat, Bulk_ESS = ess_bulk, Tail_ESS = ess_tail
  # so in brms summary(fit)[, "Rhat"] and rhat(fit) are one number.
  # This used to be rstan's classic split-R-hat under the same header
  # beside a rank-normalized rhat(), below 1 in one column and above it
  # in the other on the fit dev/reviews/20260915-brmsmatch.md measured.
  cs <- dm_case()
  skip_if_not_installed("posterior")
  s <- summary(cs$ds)
  expect_equal(colnames(s), c("mean", "sd", "2.5%", "97.5%",
                              "Rhat", "Bulk_ESS", "Tail_ESS"))
  expect_false("n_eff" %in% colnames(s))
  # the agreement itself, which is the point
  expect_equal(s[, "Rhat"], rhat(cs$ds)[rownames(s)])
  a <- posterior::subset_draws(posterior::as_draws_array(cs$ds),
                               variable = rownames(s))
  d <- posterior::summarise_draws(a, ess_bulk = posterior::ess_bulk,
                                  ess_tail = posterior::ess_tail)
  expect_equal(unname(s[, "Bulk_ESS"]), d$ess_bulk)
  expect_equal(unname(s[, "Tail_ESS"]), d$ess_tail)
  # and neff_ratio() is the smaller of the two over the draw count.
  #
  # THE DIRECTION OF THIS ONE IS LOAD-BEARING. Written as a division,
  # both sides are the same division of the same two numbers and the
  # relation is EXACT: identical() holds. Turned around into "pmin
  # equals neff_ratio times ndraws" it becomes a round trip through
  # x/N*N, which on this fit is off by 0.51 ulp on one row, 5.684e-14.
  # Do not "simplify" it into the multiplication: that form needs a
  # tolerance and this one does not.
  expect_equal(unname(neff_ratio(cs$ds)[rownames(s)]),
               unname(pmin(s[, "Bulk_ESS"], s[, "Tail_ESS"]) /
                        ndraws(cs$ds)))
})

test_that("rhat() and neff_ratio() take brms's OTHER `pars` rule", {
  # brms's rhat.brmsfit passes `variable = pars` to as_draws_array(),
  # so NULL is every variable and a string is an EXACT name. It does
  # NOT go through brms's extract_pars(), which is the rule mcmc_plot()
  # and posterior_interval() follow on this same argument name.
  cs <- dm_case()
  skip_if_not_installed("posterior")
  expect_equal(names(rhat(cs$ds, "x")), "x")
  expect_equal(names(rhat(cs$ds, c("x", "Intercept"))),
               c("x", "Intercept"))
  expect_equal(names(neff_ratio(cs$ds, "^b\\[", regex = TRUE)),
               grep("^b\\[", variables(cs$ds), value = TRUE))
  # NULL is brms's own default for this slot and means every variable
  expect_equal(rhat(cs$ds, NULL), rhat(cs$ds))
  expect_equal(neff_ratio(cs$ds, NULL), neff_ratio(cs$ds))
  # an exact name that is not there is an error, as in brms; a regular
  # expression without regex = TRUE is such a name
  expect_error(rhat(cs$ds, "^x$"))
  expect_error(neff_ratio(cs$ds, "nosuchvariable"))
  # the extract_pars rule is the OTHER methods', and they still have it
  expect_error(posterior_interval(cs$ds, 0.9),
               "must be NA or a character vector")
  expect_equal(rownames(posterior_interval(cs$ds, "^x$")), "x")
})

test_that("the bayesplot accessors read the stanfit", {
  cs <- dm_case()
  skip_if_not_installed("bayesplot")
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
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  expect_equal(unname(raw[3L, , ]),
               unname(mixture_probs(frmtmb.sample:::draws_fit_at(ds, 3L, idx))),
               tolerance = 1e-12)

  st <- pp_mixture(ds)
  expect_equal(dim(st), c(100L, 4L, 2L))
  expect_equal(unname(st[, "Estimate", 1L]),
               unname(apply(raw[, , 1L], 2L, mean)), tolerance = 1e-12)
  # the data are well separated, so the assignment is nearly certain
  expect_gt(mean(apply(st[, "Estimate", ], 1, max)), 0.95)

  # brms's own signature: `summary` is the eighth argument and the
  # second is `newdata`, with log, robust and probs alongside it
  expect_equal(pp_mixture(ds, log = TRUE, summary = FALSE), log(raw))
  rb <- pp_mixture(ds, robust = TRUE)
  expect_equal(unname(rb[, "Estimate", 1L]),
               unname(apply(raw[, , 1L], 2L, stats::median)),
               tolerance = 1e-12)
  pr <- pp_mixture(ds, probs = c(0.1, 0.9))
  expect_equal(dimnames(pr)[[2L]],
               c("Estimate", "Est.Error", "Q10", "Q90"))
  expect_error(pp_mixture(ds, dd), "does not take newdata")
  expect_equal(pp_mixture(ds, draw_ids = c(2L, 4L), summary = FALSE),
               raw[c(2L, 4L), , , drop = FALSE])
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
  if (sampler_gates_on()) {
    expect_lt(max(abs(df$estimate__ - cf$x$estimate__) /
                    pmax(df$se__, cf$x$se__, 1e-8)), 5)
  }
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

## ---- hypothesis() naming notes --------------------------------------

test_that("hypothesis() on draws gives the reserved-name note, once", {
  # A covariate named `sigma` shadows the residual SD, and hypothesis()
  # says which one it read. The note is ARMED per user-level call and
  # emitted deep inside hyp_env_vals(), so the method the user reached
  # has to arm it.
  #
  # Core used to arm it in its GENERIC. Core's exported `hypothesis`
  # now resolves to brms's generic whenever brms is loaded, so the
  # arming moved into core's own methods, and this method, which had
  # relied on the generic, lost the note in EVERY session, brms or not.
  # Measured before the fix: 1 note on the frmtmb_fit, 0 on its draws.
  #
  # fake_draws() rather than the sampler: the note is emitted while the
  # hypothesis is parsed against the fit, before any draw is read, so
  # zero draws exercise exactly the path in question and this test
  # needs no Stan build to run.
  set.seed(3)
  n <- 200
  dd <- data.frame(sigma = stats::rnorm(n),
                   g = factor(rep(1:10, length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.7 * dd$sigma +
                         stats::rnorm(10, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ sigma + (1 | g)), family = gaussian(), data = dd)
  pat <- "reads 'sigma' as the coefficient"

  # the control, on the fit, so a missing note on draws cannot be a
  # construction that never shadows anything
  on_fit <- capture_messages(hypothesis(fit, "sigma = 0"))
  expect_equal(sum(grepl(pat, on_fit, fixed = TRUE)), 1L)

  ds <- fake_draws(fit)
  on_draws <- capture_messages(
    suppressWarnings(hypothesis(ds, "sigma = 0")))
  expect_equal(sum(grepl(pat, on_draws, fixed = TRUE)), 1L)
})

test_that("every hypothesis() method arms the note itself", {
  # The structural half of the test above, and the reason it exists.
  #
  # Core guards that no borrowed generic carries work in its body
  # (tests/testthat/test-generic-collision.R), because the exported
  # generic may be brms's and then that work does not run. That guard
  # stops work coming BACK into a generic. It cannot see work that LEFT
  # the generic and did not reach every method, which is exactly what
  # happened here: the arming moved into core's two methods and this
  # package's draws method was not one of them. So this asserts where
  # the work WENT, over every hypothesis() method both packages define.
  meths <- c(
    ls(asNamespace("frmtmb"), all.names = TRUE, pattern = "^hypothesis[.]"),
    ls(asNamespace("frmtmb.sample"), all.names = TRUE,
       pattern = "^hypothesis[.]"))
  # the guard is only a guard if it found the methods it is about
  expect_true(all(c("hypothesis.frmtmb_fit", "hypothesis.frmtmb_multiple",
                    "hypothesis.frmtmb_draws") %in% meths))
  arms <- vapply(meths, function(m) {
    ns <- if (exists(m, envir = asNamespace("frmtmb.sample"),
                     inherits = FALSE)) "frmtmb.sample" else "frmtmb"
    f <- get(m, envir = asNamespace(ns), inherits = FALSE)
    b <- paste(deparse(body(f)), collapse = " ")
    # The disarm must sit in on.exit(): a disarm called straight after the
    # arm leaves the note switched off for the parse it was meant to cover,
    # and a bare call would satisfy a check that only asks it to appear.
    grepl("hyp_shadow_arm()", b, fixed = TRUE) &&
      grepl("on.exit(hyp_shadow_disarm(", b, fixed = TRUE)
  }, NA)
  expect_equal(names(arms)[!arms], character())
})
