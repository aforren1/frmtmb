# frmtmb.sample answers brms's draws-side calls with brms's OUTPUT: the
# same object brms's own installed method returns on the same draws,
# compared with identical(). dev/brmsnames-match.R is the full harness
# and dev/brmsnames-findings.md its record; this file keeps the part a
# regression would break.
#
# brms is made to answer on these draws through a shim: brms's own
# `brm(empty = TRUE)` object for the same model, given a copy of the
# sampler's stanfit whose `sim` holds the post-warmup draws under the
# names frmtmb.sample puts on them. brms's methods read nothing else.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

bo_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    skip_if_not_installed("brms")
    skip_if_not_installed("posterior")
    if (is.null(cache)) {
      set.seed(11)
      dd <- data.frame(x = stats::rnorm(120), g = factor(rep(1:8, 15)))
      u <- cbind(stats::rnorm(8, 0, 0.7), stats::rnorm(8, 0, 0.4))
      dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + u[dd$g, 1] +
                             u[dd$g, 2] * dd$x, 1)
      fit <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 400, refresh = 0, seed = 3)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

# brms's side is the model frmtmb fitted, formula for formula: sigma
# has no formula in either, so brms's `sigma` and frmtmb.sample's are
# the same natural-scale parameter
bo_shim <- function(cs, M) {
  b <- suppressWarnings(suppressMessages(
    brms::brm(y ~ x + (1 + x | g), data = cs$dd, empty = TRUE,
              backend = "rstan")))
  sf <- cs$ds$stanfit
  nc <- sf@sim$chains
  n <- nrow(M) %/% nc
  sim <- sf@sim
  sim$samples <- lapply(seq_len(nc), function(ch) {
    rows <- (ch - 1L) * n + seq_len(n)
    stats::setNames(lapply(seq_len(ncol(M)), function(j) unname(M[rows, j])),
                    colnames(M))
  })
  sim$iter <- n
  sim$warmup <- 0
  sim$warmup2 <- rep(0, nc)
  sim$n_save <- rep(n, nc)
  sim$permutation <- lapply(seq_len(nc), function(i) seq_len(n))
  sim$pars_oi <- sim$fnames_oi <- colnames(M)
  sim$dims_oi <- stats::setNames(rep(list(integer(0)), ncol(M)),
                                 colnames(M))
  sim$n_flatnames <- ncol(M)
  sf@sim <- sim
  b$fit <- sf
  b
}

# brms's ranef() reshapes its r_ draws by position in brms's storage
# order, every level of one coefficient before the next; frmtmb stores
# level by level, so the shim for group-level calls reorders them
bo_brms_order <- function(M) {
  r <- grep("^r_", colnames(M))
  nm <- colnames(M)[r]
  lev <- as.integer(sub("^r_g[[]([0-9]+),.*$", "\\1", nm))
  cf <- sub("^r_g[[][0-9]+,(.*)[]]$", "\\1", nm)
  o <- order(match(cf, unique(cf)), lev)
  M[, r] <- M[, r[o]]
  colnames(M)[r] <- nm[o]
  M
}

test_that("the accessors and summaries are brms's, identical()", {
  cs <- bo_case()
  ds <- cs$ds
  M <- ds$draws
  sh <- bo_shim(cs, M)
  # the names themselves, which a shim built from these draws cannot
  # check: brms's b_ and r_ spellings (dev/brmsnames-naming.R derives
  # them from brms's own parse)
  expect_true(all(c("b_Intercept", "b_x", "sigma",
                    "r_g[1,Intercept]", "r_g[8,x]") %in% variables(ds)))
  expect_false("b_sigma_Intercept" %in% variables(ds))
  expect_identical(variables(ds), brms:::variables.brmsfit(sh))
  expect_identical(as.matrix(ds), brms:::as.matrix.brmsfit(sh))
  expect_identical(as.array(ds, variable = "b_x"),
                   brms:::as.array.brmsfit(sh, variable = "b_x"))
  expect_identical(as.data.frame(ds), brms:::as.data.frame.brmsfit(sh))
  expect_identical(as_draws_df(ds), brms:::as_draws_df.brmsfit(sh))
  expect_identical(posterior_summary(ds),
                   brms:::posterior_summary.brmsfit(sh))
  expect_identical(fixef(ds), brms:::fixef.brmsfit(sh))
  expect_identical(fixef(ds, FALSE), brms:::fixef.brmsfit(sh, FALSE))
  expect_identical(rhat(ds), brms:::rhat.brmsfit(sh))
  h <- c("x > 0", "x = 0.5")
  expect_identical(hypothesis(ds, h)$hypothesis[1:5],
                   brms:::hypothesis.brmsfit(sh, h)$hypothesis[1:5])
  expect_identical(hypothesis(ds, h)$samples,
                   brms:::hypothesis.brmsfit(sh, h)$samples)

  # the inverse case: one draw moved, and identical() says so, so the
  # assertions above cannot pass by comparing an object with itself
  Mp <- M
  Mp[1L, "b_x"] <- Mp[1L, "b_x"] + 1e-9
  expect_false(identical(fixef(ds, FALSE),
                         brms:::fixef.brmsfit(bo_shim(cs, Mp), FALSE)))
})

test_that("ranef(), coef() and VarCorr() are brms's, draws and summaries", {
  cs <- bo_case()
  ds <- cs$ds
  shb <- bo_shim(cs, bo_brms_order(ds$draws))
  expect_identical(ranef(ds), brms:::ranef.brmsfit(shb))
  expect_identical(ranef(ds, FALSE), brms:::ranef.brmsfit(shb, FALSE))
  expect_identical(coef(ds), brms:::coef.brmsfit(shb))
  expect_identical(coef(ds, FALSE, TRUE), brms:::coef.brmsfit(shb, FALSE,
                                                               TRUE))

  # brms reads sd_ and cor_ draws; these draws store theta, so the shim
  # carries the standard deviations and correlation frmtmb derives from
  # it, and what is compared is brms's layout and summary path on them
  vd <- VarCorr(ds, summary = FALSE)
  M <- cbind(ds$draws,
             sd_g__Intercept = unname(vd$g$sd[, "Intercept"]),
             sd_g__x = unname(vd$g$sd[, "x"]),
             cor_g__Intercept__x = vd$g$cor[, "Intercept", "x"])
  expect_equal(unname(vd$g$sd[, "Intercept"]),
               unname(exp(ds$draws[, "theta_1"])))
  sv <- bo_shim(cs, M)
  expect_identical(VarCorr(ds), brms:::VarCorr.brmsfit(sv))
  expect_identical(VarCorr(ds, NULL, FALSE),
                   brms:::VarCorr.brmsfit(sv, NULL, FALSE))
})
