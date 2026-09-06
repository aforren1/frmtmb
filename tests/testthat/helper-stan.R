# A Stan-file identity harness.
#
# helper-brms.R compares frmtmb with the Stan program brms GENERATES.
# This file compares frmtmb with a Stan program somebody WROTE: the
# reference implementation of a published model. A test says
#
#   here is a Stan program, here is its data, here is the map from
#   frmtmb's estimate to its parameters,
#
# and gets back the assertion that Stan's log_prob at that point equals
# frmtmb's log density there plus a stated constant.
#
# helper-rl.R does this for one model with the mechanism written out in
# place. Everything below is that mechanism with the model taken out.
# test-bcm-*.R is the caller.
#
# ---------------------------------------------------------------------
# The constant
# ---------------------------------------------------------------------
#
# The claim is
#
#   log_prob(theta_hat) - frmtmb_log_density(theta_hat) == const
#
# and `const` is a property of the pair of programs, never a number
# read off a failing run. There are three cases and no fourth:
#
#   1. The Stan program carries a prior frmtmb does not. `const` is the
#      total log prior density AT THE ESTIMATE. State it as a function
#      of the parameters so the test reads as the arithmetic it is:
#
#        const = function(pars, data) dbeta(pars$theta, 1, 1, log = TRUE)
#
#      A flat prior on a bounded parameter is a constant of -log(width),
#      which is zero only when the width is one.
#   2. frmtmb carries the same prior through set_prior(). Both sides
#      hold the same penalty and `const` is 0.
#   3. Neither carries a prior. `const` is 0.
#
# If the measured constant is not the stated one, the two programs are
# not the same model. That is the whole point of the harness, so the
# failure message reports both numbers.
#
# One adaptation is made to every reference program: a sampling
# statement `y ~ dist(...)` becomes `target += dist_lpdf(y | ...)`.
# Stan drops the normalizing constant of a `~` statement and frmtmb's
# families do not, so `~` would push a data-dependent term such as
# lchoose(n, k) into `const` and hide a real difference behind it.
#
# ---------------------------------------------------------------------
# The gradient
# ---------------------------------------------------------------------
#
# The value alone is a weak check: a wrong map that permutes two
# coefficients of the same size still lands close. The gradient at the
# same point is what catches it.
#
#   grad = "all"    every coordinate vanishes. Valid when the two
#                   programs hold the same function of the parameters,
#                   which is case 2 or 3 above, or case 1 with a prior
#                   whose density is flat over the estimate.
#   grad = "inner"  only the coordinates of `inner` vanish. This is the
#                   check for a fit with random effects: frmtmb's outer
#                   estimate maximizes the LAPLACE-APPROXIMATED
#                   marginal, so the joint's gradient there is not zero,
#                   while the inner block is at its conditional mode and
#                   its own gradient is.
#   grad = "none"   the Stan program's priors bend the surface, or a
#                   parameter sits on a boundary the optimizer stops at
#                   rather than zeroes, and the test claims the value
#                   only. Say why in the test.
#
# The default tolerance is 1e-3, which is helper-brms.R's: it is the
# scale at which frmtmb's own optimizer declares convergence, so a
# tighter one tests the optimizer rather than the map.
#
# The default is "inner" when `inner` names something and "all"
# otherwise.

# The gate. Compiling one program costs about a minute, so this tier is
# opt-in exactly like the brms fit tier: outside R CMD check BOTH
# variables are needed.
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
skip_unless_stan_identity <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rstan")
  if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
    testthat::skip("set FRMTMB_BRMS_FIT_TESTS=true to run the Stan tier")
  }
}

# helper-brms.R owns the content-addressed program cache, and with it
# the session cache that works around a DSO which cannot be re-read from
# its own RDS in the session that wrote it. testthat sources every
# helper into one environment, so calling it here is one cache for the
# whole suite rather than two fighting over a directory.
stan_program <- function(code) brms_stan_model(code)

# A stanfit with no draws. `chains = 0` builds every object log_prob()
# needs and runs no sampler, which is the cheapest way to reach Stan's
# density from R.
stan_density <- function(code, data) {
  suppressMessages(rstan::sampling(stan_program(code), data = data,
                                   chains = 0))
}

# The names Stan's `parameters` block declares, one per block rather
# than one per coordinate. constrained_param_names() would also return
# transformed parameters and generated quantities, which a caller must
# NOT supply; the unconstrained list is exactly the block.
stan_upar_names <- function(sf) {
  sf@.MISC$stan_fit_instance$unconstrained_param_names(FALSE, FALSE)
}

stan_par_blocks <- function(sf) {
  unique(sub("[.].*$", "", stan_upar_names(sf)))
}

# The coordinates of the unconstrained vector that belong to the named
# blocks. This is how `inner` becomes a gradient index, and it is exact
# for a simplex or a Cholesky factor, where the unconstrained length is
# not the constrained one.
stan_upar_index <- function(sf, blocks) {
  which(sub("[.].*$", "", stan_upar_names(sf)) %in% blocks)
}

# frmtmb's own log density at its estimate: the log likelihood plus any
# penalty set_prior() put on it, with the random effects at their
# conditional modes. This is the JOINT, which is what a Stan program
# that declares the random effects as parameters computes. For a fit
# with no random effects it is also the marginal, so logLik(fit) agrees
# and either spelling may be used.
frm_joint_lp <- function(fit) {
  -fit$obj$env$f(fit$obj$env$last.par.best)
}

# Population-level coefficients of one linear predictor, unnamed and in
# frmtmb's own column order, which is model.matrix()'s.
frm_b <- function(fit, dpar = "mu") {
  b <- fixef(fit)[[dpar]]
  if (is.null(b)) {
    stop("frm_b(): the fit has no '", dpar, "' linear predictor; it has ",
         paste(names(fixef(fit)), collapse = ", "), call. = FALSE)
  }
  unname(b)
}

# One group's conditional modes as a plain matrix, levels by
# coefficients. ranef() carries the level labels as row names and a Stan
# program indexes subjects 1..S, so a caller that builds its own index
# must build it from the same factor the frame saw.
frm_u <- function(fit, group = 1) {
  u <- ranef(fit)[[group]]
  matrix(as.numeric(as.matrix(u)), nrow = nrow(u))
}

# One random-effect block, keyed by the TERM LABEL VarCorr() prints.
#
# That label is the only unambiguous key. ranef() names every block
# after its GROUPING FACTOR, so a model with a random intercept on two
# distributional parameters has two blocks both called "id" whose
# columns are both called "(Intercept)"; `||` splits one term into
# several blocks with the same name again. VarCorr() distinguishes them
# ("1 | id" against "phi: 1 | id", "0 + half | id" against
# "0 + bias | id"), and its list runs in the same order as ranef()'s, so
# one index serves both.
frm_block_index <- function(fit, term) {
  nm <- names(VarCorr(fit))
  i <- match(term, nm)
  if (is.na(i)) {
    stop("frm_block_index(): no random-effect term '", term,
         "'; the fit has ", paste(nm, collapse = ", "), call. = FALSE)
  }
  i
}

# The conditional modes of one column of one block.
frm_u_term <- function(fit, term, col = 1L) {
  unname(ranef(fit)[[frm_block_index(fit, term)]][, col])
}

# That column's standard deviation.
frm_sd_term <- function(fit, term, col = 1L) {
  v <- unname(VarCorr(fit)[[frm_block_index(fit, term)]])
  sqrt(v[col, col])
}

# The lower Cholesky factor of one group's covariance. A Stan program
# that wants a correlated random effect and no constrained parameter
# takes this as DATA, which is helper-rl.R's arrangement and the one
# test-bcm-mpt.R uses: the gradient claim there is about the latent
# traits, and the outer gradient of the joint is not zero at a Laplace
# optimum anyway.
frm_chol <- function(fit, group = 1) {
  t(chol(unname(VarCorr(fit)[[group]])))
}

# The identity, at one point.
#
#   code    the Stan program, as a string
#   data    its data list
#   pars    the values of its `parameters` block, in frmtmb's estimate.
#           A function is called with the fit, which keeps the map next
#           to the fit that feeds it.
#   fit     the frmtmb fit, for the default `ours` and for the report
#   ours    frmtmb's log density; defaults to frm_joint_lp(fit)
#   const   the stated difference, a number or function(pars, data)
#   inner   the parameter blocks holding the inner (random-effect)
#           coordinates
#   grad    "all", "inner" or "none"; see the note at the top
stan_lp_check <- function(code, data, pars, fit = NULL, ours = NULL,
                          const = 0, inner = character(), grad = NULL,
                          tol = 1e-6, tol_grad = 1e-3, roundtrip = TRUE) {
  sf <- stan_density(code, data)
  if (is.function(pars)) {
    pars <- pars(fit)
  }
  # A missing block surfaces from rstan as a length complaint about the
  # whole vector, which names nothing. Refuse first, by name.
  missing <- setdiff(stan_par_blocks(sf), names(pars))
  if (length(missing)) {
    stop("stan_lp_check(): the map does not cover ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  extra <- setdiff(names(pars), stan_par_blocks(sf))
  if (length(extra)) {
    stop("stan_lp_check(): the map names ", paste(extra, collapse = ", "),
         ", which the parameters block does not declare", call. = FALSE)
  }

  upars <- rstan::unconstrain_pars(sf, pars)
  # A round trip through Stan's own constraint machinery is the unit
  # test on the map alone: a simplex that does not sum to one, or a
  # covariance that is not positive definite, comes back changed even
  # though log_prob() would still return a number.
  if (roundtrip) {
    back <- rstan::constrain_pars(sf, upars)
    for (nm in names(pars)) {
      testthat::expect_lt(
        max(abs(as.numeric(back[[nm]]) - as.numeric(pars[[nm]]))), 1e-10,
        label = paste0("round trip of ", nm))
    }
  }

  # `adjust_transform = FALSE` is the whole reason a constant can be
  # STATED rather than measured. With it, both sides are densities in
  # the CONSTRAINED parameterization, which is where frmtmb's log
  # density lives; `TRUE` would add the log Jacobian of Stan's own
  # constraining transform to one side only, and every constant in
  # every caller would carry that Jacobian instead of the prior.
  lp <- rstan::log_prob(sf, upars, adjust_transform = FALSE,
                        gradient = FALSE)
  if (is.null(ours)) {
    ours <- frm_joint_lp(fit)
  }
  if (is.function(const)) {
    const <- const(pars, data)
  }
  measured <- lp - ours
  testthat::expect_lt(
    abs(measured - const), tol * max(1, abs(ours)),
    label = sprintf("stated constant %.10g, measured %.10g", const,
                    measured))

  grad <- grad %||% if (length(inner)) "inner" else "all"
  gr <- rstan::grad_log_prob(sf, upars, adjust_transform = FALSE)
  if (identical(grad, "inner")) {
    idx <- stan_upar_index(sf, inner)
    testthat::expect_gt(length(idx), 0)
    gr <- gr[idx]
  }
  if (!identical(grad, "none")) {
    testthat::expect_lt(max(abs(gr)), tol_grad)
  }
  # a program with an empty parameters block (Planes) has no gradient,
  # and max() of nothing is -Inf with a warning
  mg <- if (length(gr)) max(abs(gr)) else NA_real_

  # A row that drifts should report a number, not only a fail, and the
  # vignette's table is built from these. Printing is opt-in:
  #   options(frmtmb.bcm_report = TRUE)
  if (isTRUE(getOption("frmtmb.bcm_report", FALSE))) {
    cat(sprintf("BCM const stated %.10g measured %.10g grad %.3g ours %.10g\n",
                const, measured, mg, ours))
  }
  invisible(list(lp = lp, ours = ours, const = const,
                 measured_const = measured, max_grad = mg,
                 sf = sf, pars = pars))
}
