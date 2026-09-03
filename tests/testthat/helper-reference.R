# Comparison helpers against reference packages. Tolerances are absolute.

#' @srrstats {G5.4} Correctness is tested by comparison against exact
#'   reference implementations, not against stored snapshots of our own
#'   output. Every model class the package supports is checked against a
#'   package that implements the same likelihood: glmmTMB for the
#'   covariance structures, the zero-inflated and hurdle families and the
#'   distributional models; lme4 for the classical LMM and GLMM surface;
#'   mgcv for smooths and functional terms; MASS for `polr` ordinal and
#'   `glm.nb`; survival for censored responses; quantreg for
#'   `asym_laplace`; nnet for `multinomial`; mclust for the multivariate
#'   gaussian mixtures; nlme for nonlinear standard errors; GLMMadaptive
#'   as a second adaptive-quadrature implementation; mice for
#'   multiple-imputation pooling; and `stats::lm`/`glm` for the
#'   no-random-effect cases. The model-building layer is separately
#'   checked against `brms::make_standata()`.
#' @srrstats {G5.4a} The parts with no reference implementation are tested
#'   against a hand-written likelihood instead: the test builds the
#'   log-density directly with `RTMB::MakeADFun()`, optimizes it, and
#'   requires the package to reach the same optimum. About twenty test
#'   files use this route, which separates correctness of the
#'   implementation from correctness of the method by deriving the answer
#'   twice, independently.
#' @srrstats {G5.4b} These are tests against previous implementations of
#'   existing methods. The reference packages are called live, at the
#'   version installed, rather than compared against transcribed numbers,
#'   so the comparison cannot go stale silently.
#' @srrstats {G5.5} Correctness tests run with a fixed random seed.
#'   `set.seed()` is called with a literal integer in 55 of the 61 test
#'   files, and the shared simulators defined here (`sim_ar1_data()`,
#'   `sim_pois_glmm()`) seed themselves from a default argument, so a
#'   test that uses one is deterministic even without its own call. The
#'   remaining files are structurally deterministic or use a standard
#'   data set.
#' @srrstats {G5.6,G5.6a} Parameter recovery is tested: data is simulated
#'   from known parameters and the fit is required to return them. This
#'   covers truncated and censored responses, `compois`, the prior and
#'   bounds machinery, `frm_simulate()`, nonlinear models, the spatial
#'   `car`/`bym2` fields, the latent-class mixtures, and the draws
#'   surface. Recovery is always judged within an explicit numeric
#'   tolerance, never by exact equality, and the tolerance is written at
#'   the call site.
#' @noRd
NULL

# Load glmmTMB's namespace once, quietly. On CI runners whose binary
# glmmTMB was built against an older TMB, its .onLoad warns about the
# version skew at the first skip_if_not_installed() of every file. The
# agreement tests themselves are the real check on glmmTMB's numbers:
# a skew that mattered would fail them at 1e-6.
suppressWarnings(requireNamespace("glmmTMB", quietly = TRUE))

expect_loglik_equal <- function(fit, ref, tol = 1e-6) {
  expect_lt(abs(as.numeric(stats::logLik(fit)) -
                  as.numeric(stats::logLik(ref))), tol)
}

expect_vector_equal <- function(x, y, tol) {
  expect_equal(length(x), length(y))
  expect_lt(max(abs(unname(x) - unname(y))), tol)
}

# Grouped AR(1) data shared across tests.
sim_ar1_data <- function(seed = 51, n_g = 60, n_t = 6, rho = 0.6,
                         sd_re = 1, sd_res = 0.5) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_t))
  tim <- factor(rep(seq_len(n_t), n_g))
  u <- replicate(n_g, {
    e <- rnorm(n_t, 0, sd_re)
    for (t in 2:n_t) e[t] <- rho * e[t - 1] + sqrt(1 - rho^2) * rnorm(1, 0, sd_re)
    e
  })
  data.frame(y = 1 + as.vector(u) + rnorm(n_g * n_t, 0, sd_res),
             g = g, tim = tim)
}

# Simulated poisson GLMM data shared across tests.
sim_pois_glmm <- function(seed = 101, n_g = 30, n_per = 20) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- stats::rnorm(n_g * n_per)
  b0 <- stats::rnorm(n_g, 0, 0.5)
  b1 <- stats::rnorm(n_g, 0, 0.3)
  eta <- 0.3 + 0.4 * x + b0[g] + b1[g] * x
  data.frame(y = stats::rpois(length(eta), exp(eta)), x = x, g = g)
}

# Chain-agreement gates: assertions that compare NUTS chain output to a
# reference quantity (the ML fit, a Wald SE, another chain). A seeded
# Stan chain is not platform-deterministic, and pkgcheck's container
# repeatedly draws chains that fail gates every vetted platform passes
# (local per-file, local single-session, the R-CMD-check runners), so
# that one workflow turns the gates off with FRMTMB_SAMPLER_GATES=false.
# Structural sampler assertions (dimensions, names, exact per-draw
# identities) never use this switch and run everywhere.
sampler_gates_on <- function() {
  !identical(Sys.getenv("FRMTMB_SAMPLER_GATES"), "false")
}
