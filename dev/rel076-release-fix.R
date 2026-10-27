# Release fixes of the 0.76.0 consolidation, sourced by
# dev/rel076-resolve.R (local = TRUE, so rd(), wr() and edit_file() are
# the resolver's), after the textual merge and the hand resolutions. Each
# edit checks that its old text is there exactly once, so a merge rerun
# applies them all or stops at the first that no longer fits. They come
# from the lanes' reviews (dev/reviews/2026-10-27-<lane>.md), named by
# the coordinator as fixes for the release. The NEWS edit of fact the
# acprior review asked for is in dev/rel076-news-split.R, because the
# NEWS merge reads the lanes' own text after this file runs.
#
#   (sourced; not run on its own)

# ---- R1. ?set_prior: a prior on the sigma intercept under sigma ~ x -------
# acprior's review, fix round 2 check, consolidation note 2: under
# `sigma ~ x` with a first-order cov = TRUE ar(), ma() or arma(), brms's
# Intercept_sigma is the log INNOVATION sd (sigma = exp(Intercept_sigma
# + ...) enters normal_time_het_lpdf() with the AR correlation divided by
# 1 - ar^2), while frmtmb's sigma intercept is its log MARGINAL sd; a
# prior there is read on frmtmb's scale (the backlog has the fix). Its
# slopes are the same quantity in both, since log c is a shift that does
# not depend on x (dev/review-acprior-r2-brmssigx.R in the lane's tree).
edit_file("R/priors.R",
"#' innovation sd moves with `ar` and is no box on the marginal sd frmtmb
#' holds.
#'
#' @section The LKJ prior:",
"#' innovation sd moves with `ar` and is no box on the marginal sd frmtmb
#' holds. Under `sigma ~ x` with such a term, a prior on the intercept of
#' `sigma` (`class = \"Intercept\", dpar = \"sigma\"`) is still read at
#' frmtmb's log marginal sd, where brms reads it at the log innovation
#' sd; the slopes (`class = \"b\", dpar = \"sigma\"`) are the same in both.
#'
#' @section The LKJ prior:")

# ---- R2. R/influence.R: the comment above dfbetas.frmtmb_influence() ------
# mocode's review, fix round 2 check, nit 4: since that lane, get_coef()
# carries the mo() simplex for marginaleffects, and dfbeta() reads the
# vector without it (interop_coef_vector(simplex = FALSE)), so the
# comment named the wrong vector
edit_file("R/influence.R",
"  # dfbeta()'s columns are get_coef()'s vector, which is every
  # estimated coefficient plus an ordinal fit's thresholds, not brms's
  # population-level block. interop_vcov() is the covariance of that
  # exact vector; vcov_estimated() is one block of it, and dividing by
  # a shorter vector recycles rather than erroring",
"  # dfbeta()'s columns are interop_coef_vector(simplex = FALSE): every
  # estimated coefficient plus an ordinal fit's thresholds, without the
  # mo() simplex that get_coef() carries for marginaleffects, and not
  # brms's population-level block. interop_vcov(simplex = FALSE) is the
  # covariance of that exact vector; vcov_estimated() is one block of
  # it, and dividing by a shorter vector recycles rather than erroring")

# ---- R3. test-mo-coding.R: the marginaleffects bound as a ratio -----------
# mocode's review, fix round 2 check, nit 4: `max(abs(se / ex - 1)) <
# 1e-4` was a fixed bound. The quantity the test exists to tell apart is
# what the simplex adds to a standard error: the run measures it as the
# distance from the exact standard error to the one that leaves the
# simplex out (0.75.0's, from the same basis and covariance), 0.306 at
# k = 2 on this seed, against a marginaleffects error of about 1e-7
# (dev/rel076-me-probe.R, dev/rel076-log/me-probe.txt). The basis and
# covariance are checked to give core's own standard error first, so the
# leave-out arm is the same computation with two columns removed.
edit_file("tests/testthat/test-mo-coding.R",
"  se <- marginaleffects::predictions(fit, newdata = nd)$std.error
  # marginaleffects' Jacobian is a finite difference: 1e-6 relative here
  expect_lt(max(abs(se / ex - 1)), 1e-4)",
"  se <- marginaleffects::predictions(fit, newdata = nd)$std.error
  # the standard error that leaves the simplex out, 0.75.0's, from the
  # basis and covariance of core's delta method: the yardstick for how
  # close marginaleffects must be
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
  jc <- get_joint_cov(fit)
  A <- as.matrix(lb$A)
  Vb <- as.matrix(jc$V)[lb$coef_pos, lb$coef_pos, drop = FALSE]
  zc <- which(jc$names[lb$coef_pos] == \"zeta1\")
  expect_length(zc, 2L)
  own <- frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit
  expect_equal(unname(sqrt(rowSums((A %*% Vb) * A))), unname(own),
               tolerance = sqrt(.Machine$double.eps))
  A0 <- A[, -zc, drop = FALSE]
  se0 <- sqrt(rowSums((A0 %*% Vb[-zc, -zc, drop = FALSE]) * A0))
  # marginaleffects' Jacobian is a finite difference, about 1e-7
  # relative here; the simplex moves the standard error by 0.31
  expect_lt(max(abs(se / ex - 1)), 1e-4 * max(abs(se0 / ex - 1)))")

# ---- after the release review: R4 to R9 -------------------------------------
# dev/reviews/2026-10-27-release.md (FIX, no blocker or major); the
# fixes and their reasons are dev/rel076-release-fix2.R
source("dev/rel076-release-fix2.R", local = TRUE)
