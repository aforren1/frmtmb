# A new group level in a conditional-effects grid is drawn by brms's
# sample_new_levels rule: "uncertainty" (brms's default) lends an
# observed level's effects per draw or bootstrap replicate, "old_levels"
# lends one level for the whole call, "gaussian" draws from the
# estimated covariance.

snl_fit <- function() {
  set.seed(11)
  dg <- data.frame(x = stats::rnorm(48), g = factor(rep(1:6, each = 8)))
  dg$y <- 1 + 0.5 * dg$x + stats::rnorm(6, sd = 0.8)[dg$g] +
    stats::rnorm(48, sd = 0.5)
  frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
}

snl_draws <- function(fit, plan, k = 20L) {
  vapply(seq_len(k), function(i) {
    cache <- new.env(parent = emptyenv())
    ce_plan_eval(fit, plan, 1L, cache, FALSE, "y", "mu", NULL, TRUE)
    unlist(as.list(cache), use.names = FALSE)
  }, 1)
}

test_that("a new level lends an observed level's effects by default", {
  fit <- snl_fit()
  lp <- find_linpred(fit, "y", "mu")
  gb <- ce_grids_build(fit, fit$spec$responses[[1L]], lp, "x", "y", "mu", 3,
                       list(), NULL, na_vars = ce_group_vars(fit))
  modes <- fit$estimates$b[fit$frame$re_blocks[[1L]]$b_idx]
  set.seed(3)
  plan <- ce_level_plan(fit, gb$grids, gb$base, NULL)
  expect_identical(plan$snl, "uncertainty")
  u <- snl_draws(fit, plan)
  # every draw is one of the six fitted modes, and more than one of them
  # is lent over 20 draws (the chance of one mode 20 times is 6^-19)
  expect_true(all(u %in% modes))
  expect_gt(length(unique(u)), 1L)
  plan_o <- ce_level_plan(fit, gb$grids, gb$base, NULL, "old_levels")
  o <- snl_draws(fit, plan_o)
  expect_true(all(o %in% modes))
  expect_length(unique(o), 1L)
  plan_g <- ce_level_plan(fit, gb$grids, gb$base, NULL, "gaussian")
  gs <- snl_draws(fit, plan_g)
  expect_false(any(gs %in% modes))
})

test_that("a fit's display draws a new level by the gaussian law", {
  # a refit holds each level's mode, not a draw of it, so the fit keeps
  # "gaussian", as predict() and fitted() do, and refuses the rules that
  # lend a level
  fit <- snl_fit()
  ic <- list(x = c(-1, 0, 1))
  ce_d <- suppressMessages(conditional_effects(
    fit, "x", re_formula = NULL, band = "boot", boot = 4, seed = 2,
    int_conditions = ic))
  ce_g <- suppressMessages(conditional_effects(
    fit, "x", re_formula = NULL, band = "boot", boot = 4, seed = 2,
    int_conditions = ic, sample_new_levels = "gaussian"))
  expect_identical(ce_d$x, ce_g$x)
  for (snl in c("uncertainty", "old_levels")) {
    expect_error(conditional_effects(fit, "x", re_formula = NULL,
                                     band = "boot", boot = 4,
                                     int_conditions = ic,
                                     sample_new_levels = snl),
                 "cannot carry|carry none of its conditional variance",
                 info = snl)
    # a band that draws nothing refuses a rule it would ignore
    expect_error(conditional_effects(fit, "x", re_formula = NULL,
                                     int_conditions = ic,
                                     sample_new_levels = snl),
                 "this band draws no new level", fixed = TRUE, info = snl)
  }
  # the guard's absent case: the Wald band takes its own law, named
  ce_w <- conditional_effects(fit, "x", re_formula = NULL,
                              int_conditions = ic)
  expect_identical(conditional_effects(fit, "x", re_formula = NULL,
                                       int_conditions = ic,
                                       sample_new_levels = "gaussian"),
                   ce_w)
  expect_error(conditional_effects(fit, "x", re_formula = NULL,
                                   sample_new_levels = "nearest"),
               "sample_new_levels")
})

test_that("a fit's boot band at a new level is as wide as the Wald band", {
  # strong shrinkage, where a lent mode is much narrower than a level's
  # effect: 30 groups of 2 rows, sd_g 1, sigma 2, data seed 100 (ML sd_g
  # 0.941, the 30 modes' sd 0.603). At 400 refits brms 2.23.0's default
  # band was 4.096 3.963 3.922 wide, frmtmb's Wald band 3.967 3.840
  # 3.940, its gaussian boot band 4.036 4.344 4.560, and a boot band that
  # lent the replicate's modes 2.903 2.815 2.787 (ratio to Wald 0.71 to
  # 0.73; dev/reviews/2026-11-13-parity2.md, F1). Over 30 boot seeds on
  # these data the minimum ratio over the three x values has mean 0.985
  # and sd 0.086, so 0.8 is about 2.1 sd below its mean (0 of 30 seeds
  # at or below it; dev/reviews/2026-11-13-parity2-short.md), and above
  # the lent band's 0.73. The fixed boot seed (2, ratio 0.961) carries
  # it, and 20 ulp perturbations of y keep 0.961 on both BLAS.
  skip_on_cran()
  set.seed(100)
  d <- data.frame(x = stats::rnorm(60), g = factor(rep(1:30, each = 2)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(30)[d$g] + stats::rnorm(60, sd = 2)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
  ic <- list(x = c(-1, 0, 1))
  wid <- function(ce) ce$x$upper__ - ce$x$lower__
  wb <- wid(suppressMessages(conditional_effects(
    fit, "x", re_formula = NULL, band = "boot", boot = 200, seed = 2,
    int_conditions = ic)))
  ww <- wid(conditional_effects(fit, "x", re_formula = NULL,
                                int_conditions = ic))
  expect_gt(min(wb / ww), 0.8)
})
