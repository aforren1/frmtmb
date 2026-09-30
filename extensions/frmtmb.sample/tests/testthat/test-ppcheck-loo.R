# pp_check()'s loo_* types on draws. bayesplot's ppc_loo_* functions
# need the PSIS weights (`lw`) or the psis object, which
# brms:::pp_check.brmsfit() builds with loo(save_psis = TRUE) on the same
# draws as the predictions. The draws method built neither, so 4 of 6
# types failed on 3 of 3 models (dev/arcovsample-rev-07-ppcheck.R, seed
# 31). Against a compiled brms fit, with brms's own draws handed to both,
# the weights bayesplot receives are identical, max |diff| 0, for all
# draws and for a subset (dev/sampfix-05-ppcheck-brms.R).

ppc_data <- function() {
  set.seed(31L)
  dd <- data.frame(g = factor(rep(1:5, each = 6L)), t = rep(1:6, 5L))
  dd$x <- stats::rnorm(nrow(dd))
  dd$y <- 0.6 + 0.5 * dd$x + stats::rnorm(nrow(dd), 0, 0.8)
  dd
}

# draws around a plausible point, without a sampler: the weights are a
# function of log_lik() on whatever draws there are
ppc_draws <- function(fit, nd = 200L) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(0, nd, length(lab), dimnames = list(NULL, lab))
  set.seed(77L)
  for (j in seq_along(lab)) {
    nm <- lab[j]
    m[, j] <- if (nm == "lp__") 0
    else if (grepl("^sigma", nm)) exp(stats::rnorm(nd, log(0.8), 0.05))
    else if (grepl("^thetaac", nm)) stats::rnorm(nd, 0.4, 0.05)
    else stats::rnorm(nd, 0.5, 0.1)
  }
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

ppc_loo_types <- c("loo_pit_overlay", "loo_pit_qq", "loo_intervals",
                   "loo_ribbon")

# a small design with independent draws makes loo warn about k and
# posterior cap the relative efficiency
ppc_ok_warnings <- c("Pareto k", "pareto_k", "Relative effective",
                     "ESS has been capped")

test_that("every loo_* type draws, on three kinds of model", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("loo")
  dd <- ppc_data()
  forms <- list(bf(y ~ x), bf(y ~ x + (1 | g)), bf(y ~ x + ar(t, g)))
  for (f in forms) {
    fit <- frm(f, family = gaussian(), data = dd, dry_run = "objective")
    ds <- ppc_draws(fit)
    for (ty in ppc_loo_types) {
      p <- allow_warnings(suppressMessages(pp_check(ds, type = ty)),
                          ppc_ok_warnings)
      expect_s3_class(p, "ggplot")
    }
  }
})

test_that("bayesplot receives loo's PSIS weights of log_lik()", {
  skip_if_not_installed("loo")
  dd <- ppc_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
             dry_run = "objective")
  ds <- ppc_draws(fit)
  # the arguments bayesplot would be called with, in place of the plot
  local_mocked_bindings(draws_bayesplot_fun = function(nm, what) {
    real <- get(nm, asNamespace("bayesplot"))
    have <- intersect(c("y", "yrep", "lw", "psis_object"),
                      names(formals(real)))
    got <- lapply(have, function(a) {
      bquote(if (!missing(.(as.name(a)))) .(as.name(a)))
    })
    names(got) <- have
    f <- real
    body(f) <- as.call(c(as.name("list"), got))
    f
  })
  expect_weights <- function(ids) {
    ll <- log_lik(ds, draw_ids = ids)
    cid <- if (is.null(ids)) rep(1L, nrow(ll)) else rep(1L, length(ids))
    ps <- loo::psis(-ll, r_eff = loo::relative_eff(exp(ll), chain_id = cid))
    a <- allow_warnings(suppressMessages(
      pp_check(ds, type = "loo_pit_overlay", draw_ids = ids)),
      ppc_ok_warnings)
    expect_identical(a$lw, stats::weights(ps, log = TRUE))
    expect_identical(a$y, dd$y)
    b <- allow_warnings(suppressMessages(
      pp_check(ds, type = "loo_intervals", draw_ids = ids)),
      ppc_ok_warnings)
    expect_identical(stats::weights(b$psis_object, log = TRUE),
                     stats::weights(ps, log = TRUE))
    expect_identical(nrow(a$yrep), nrow(ll))
  }
  allow_warnings(expect_weights(NULL), ppc_ok_warnings)
  allow_warnings(expect_weights(seq(1L, 200L, by = 3L)), ppc_ok_warnings)
})

test_that("the weights use the chain structure only when every draw is used", {
  skip_if_not_installed("loo")
  dd <- ppc_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
             dry_run = "objective")
  ds <- ppc_draws(fit)
  # brms's r_eff_log_lik(): relative_eff() on the chains when all draws
  # are used, on one chain for a subset. Two chains of 100 here; the
  # draws are independent, so the two rules differ only through the
  # relative efficiency, which is enough to tell them apart
  local_mocked_bindings(draws_nchains = function(x) 2L)
  local_mocked_bindings(draws_bayesplot_fun = function(nm, what) {
    function(y, yrep, lw, ...) lw
  })
  ll <- log_lik(ds)
  rule <- function(ll, cid) {
    stats::weights(loo::psis(-ll, r_eff = loo::relative_eff(
      exp(ll), chain_id = cid)), log = TRUE)
  }
  two <- allow_warnings(rule(ll, rep(1:2, each = 100L)), ppc_ok_warnings)
  one <- allow_warnings(rule(ll, rep(1L, 200L)), ppc_ok_warnings)
  expect_gt(max(abs(two - one)), 0)
  lw <- allow_warnings(suppressMessages(
    pp_check(ds, type = "loo_pit_overlay")), ppc_ok_warnings)
  expect_identical(lw, two)
  ids <- seq(1L, 200L, by = 3L)
  lw_sub <- allow_warnings(suppressMessages(
    pp_check(ds, type = "loo_pit_overlay", draw_ids = ids)),
    ppc_ok_warnings)
  expect_identical(lw_sub, allow_warnings(rule(ll[ids, , drop = FALSE],
                                rep(1L, length(ids))), ppc_ok_warnings))
})

test_that("the loo_* types refuse where log_lik() refuses, naming both", {
  dd <- ppc_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
             dry_run = "objective")
  ds <- ppc_draws(fit)
  expect_error(suppressMessages(
    pp_check(ds, type = "loo_pit_overlay", newdata = dd)),
    "pp_check(type = 'loo_pit_overlay') weights the draws", fixed = TRUE)
  lap <- ds
  lap$draws <- ds$draws[, !startsWith(colnames(ds$draws), "r_"),
                        drop = FALSE]
  # laplace-shaped draws: the weights are refused before any prediction
  expect_error(suppressMessages(pp_check(lap, type = "loo_intervals")),
               "log_lik() refused: log_lik() needs draws of the random effects",
               fixed = TRUE)
  # and re_formula = NA, which the predictions could honor, is still a
  # log_lik() refusal: its density is conditional on the group effects
  expect_error(suppressMessages(
    pp_check(lap, type = "loo_intervals", re_formula = NA)),
    "log_lik() refused: log_lik() does not take re_formula", fixed = TRUE)
})

test_that("a type bayesplot no longer lists is refused, as brms refuses it", {
  skip_if_not_installed("bayesplot")
  dd <- ppc_data()
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd,
             dry_run = "objective")
  expect_error(pp_check(ppc_draws(fit), type = "loo_pit"),
               "Type 'loo_pit' is not a valid ppc type", fixed = TRUE)
})
