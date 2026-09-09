# The static refusal of a tmbstan build that samples the wrong density,
# and the skip that keeps the gated tier from erroring on such a build.
#
# This file replaces the block that lived in test-sampling-ported.R.
# That block asserted the refusal PATTERN by grepping a marker out of a
# string it had just written, which is true of the literal and says
# nothing about the guard: the stop() branch was never executed and the
# comment said so. Here the broken build is constructed instead, with
# with_tmbstan_hpp() in helper-sampling.R, so both branches run on a
# machine whose own tmbstan is healthy.

skip_on_cran()

test_that("this installation's tmbstan is a healthy build", {
  skip_if_not_installed("tmbstan")
  # A canary, not a skip. If this machine's tmbstan samples a standard
  # normal, every draws assertion in this package would be measuring
  # garbage, and the tier below skips rather than reports it. That
  # state must be loud, so this block FAILS rather than skipping.
  expect_false(frmtmb.sample:::tmbstan_build_broken())
  expect_silent(frmtmb.sample:::check_tmbstan_build("frm_sample()"))
})

test_that("a generated model.hpp with the unpatched overload is refused", {
  with_tmbstan_hpp(broken = TRUE, {
    expect_true(frmtmb.sample:::tmbstan_build_broken())
    expect_error(
      frmtmb.sample:::check_tmbstan_build("frm_sample()"),
      "samples a standard normal", fixed = TRUE)
    expect_error(
      frmtmb.sample:::check_tmbstan_build("as_tmbstan()"),
      "as_tmbstan()", fixed = TRUE)
  })
  # the state the fixture borrowed is handed back
  expect_false(frmtmb.sample:::tmbstan_build_broken())
})

test_that("a patched model.hpp is not refused", {
  # the arm that matters more than the refusal: a check that fires on a
  # healthy build would refuse every correct installation.
  #
  # Two overloads, both patched. That is NOT the stanc 2.32 output,
  # which emits ONE overload; it is what a FIXED autogen would emit
  # under 2.39, which is the harder case for the detector and the one
  # worth pinning. The 2.32 shape is covered by the block below,
  # against this installation's own model.hpp.
  with_tmbstan_hpp(broken = FALSE, {
    expect_false(frmtmb.sample:::tmbstan_build_broken())
    expect_silent(frmtmb.sample:::check_tmbstan_build("frm_sample()"))
  })
})

test_that("skip_sampler() skips on a broken build and runs on a clean one", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  # a testthat skip is a condition, so it can be caught and asserted on
  # without this block itself skipping
  caught <- with_tmbstan_hpp(broken = TRUE, {
    tryCatch({
      skip_sampler()
      NULL
    }, skip = function(cnd) cnd)
  })
  expect_s3_class(caught, "skip")
  expect_match(conditionMessage(caught), "StanHeaders", fixed = TRUE)

  # and the arm every guard in the 0.55.1 round got wrong: a gate that
  # skipped unconditionally would look identical to this one on a
  # healthy machine, so the clean build must reach the code AFTER the
  # skip
  ran <- with_tmbstan_hpp(broken = FALSE, {
    tryCatch({
      skip_sampler()
      "ran"
    }, skip = function(cnd) "skipped")
  })
  expect_identical(ran, "ran")
})

test_that("the detector reads model.hpp rather than the runtime StanHeaders", {
  skip_if_not_installed("tmbstan")
  # The version of StanHeaders INSTALLED is irrelevant: model.hpp is
  # generated when tmbstan is compiled and shipped inside it. This
  # machine carries StanHeaders 2.39.1 beside a clean CRAN binary
  # tmbstan, which is the case a version comparison would have called
  # broken, and is also why "install a binary tmbstan build" is the
  # working remedy the refusal names.
  skip_if_not_installed("StanHeaders")
  sh <- utils::packageVersion("StanHeaders")
  broken <- frmtmb.sample:::tmbstan_build_broken()
  hpp <- system.file("model.hpp", package = "tmbstan")
  expect_true(nzchar(hpp))
  expect_identical(
    broken,
    any(grepl("std_normal_lpdf<propto__>(y)",
              readLines(hpp, warn = FALSE), fixed = TRUE)))
  if (sh >= "2.39.0" && !broken) {
    succeed(paste("StanHeaders", sh,
                  "installed beside a clean tmbstan build"))
  }
})

test_that("the reverse-mode gradient is the model's, not a standard normal", {
  skip_sampler()
  skip_if_not_installed("rstan")
  # The check that does not depend on the marker string.
  #
  # The two generated log_prob_impl overloads are selected by
  # stan::require_not_st_var (the one autogen patches) and
  # stan::require_st_var (the one it misses), so on an affected build
  # the entry points that matter read a standard normal instead of the
  # model. HMC takes those, which is why the defect is silent. This
  # asserts on them directly.
  #
  # A chain is started only because a stanfit is the object that
  # carries those entry points; NOTHING here reads its draws. The
  # identity holds at every point in the parameter space, so the
  # assertion does not depend on chain luck, on the seed reproducing,
  # or on platform determinism, and it says nothing about how the
  # placeholder is spelled. That last part is its whole value: an
  # upstream that renames the placeholder leaves
  # tmbstan_build_broken() answering FALSE on a sampler that is still
  # wrong, and this is what would still notice (measured:
  # dev/tmbstan-findings.md).
  #
  # It runs only on builds the static check passes, because
  # skip_sampler(), the first line of this block, skips first. That is
  # the point: this covers the gap the static check leaves, not the
  # same ground twice.
  #
  # The identity assumes the sampled parameters are UNBOUNDED. Stan
  # leaves lp__ untouched when both bounds are infinite, which is what
  # makes these two gradients the same transported adjoint rather than
  # two computations that happen to agree. A bounded parameter adds a
  # Jacobian term and the equality below would not hold.
  set.seed(4021)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- stats::rnorm(80, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  # fixed effects on purpose. With a random block fit$obj marginalizes
  # by Laplace while tmbstan samples the random effects too, so the two
  # parameter vectors are not the same object and the identity below
  # does not apply. The defect is in the generated model rather than in
  # any one tape, so one fixed-effects model tests it fully.
  sf <- suppressWarnings(suppressMessages(
    as_tmbstan(fit, chains = 1L, iter = 20L, warmup = 10L,
               refresh = 0, seed = 11)))
  np <- rstan::get_num_upars(sf)
  expect_equal(np, length(fit$obj$par))

  set.seed(4030)
  pts <- list(rep(0, np), fit$obj$par + 0.3, fit$obj$par - 0.7,
              stats::rnorm(np))
  for (p in pts) {
    u <- as.numeric(p)
    # obj$fn is the NEGATIVE log posterior and log_prob is the log
    # posterior, hence the sign
    g_stan <- as.numeric(rstan::grad_log_prob(sf, u))
    g_obj <- -as.numeric(fit$obj$gr(u))
    # what an affected build returns, in closed form rather than from
    # memory: its density is the standard normal kernel -sum(u^2)/2,
    # so its gradient is -u. dev/prior-dropping-investigation.md
    # recorded exactly that from the affected container, where stan_lp
    # came back as -mu^2/2 and stan_gr as -mu.
    g_bad <- -u
    scale <- max(abs(g_obj))
    expect_gt(scale, 0)
    # No absolute tolerance. Both sides are ratios this run measures,
    # and the separation is not marginal: over four model shapes and
    # twenty points the correct gap was 0 and bitwise, 0 ulp, while the
    # defect's was 0.963 to 0.994 of the gradient's own size. A factor
    # of a million is far inside that.
    expect_lt(max(abs(g_stan - g_obj)) / scale,
              max(abs(g_bad - g_obj)) / scale / 1e6)
  }

  # the value entry point discriminates too, per the same recorded
  # container table, so it is asserted rather than left to rot
  u <- as.numeric(fit$obj$par + 0.3)
  expect_equal(as.numeric(rstan::log_prob(sf, u)),
               -as.numeric(fit$obj$fn(u)))
})
