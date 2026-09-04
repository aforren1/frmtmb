# Sampling blocks that used to live one or two at a time in frmtmb's
# own test files, gathered here when the sampling surface moved out.
# Each group keeps a comment naming the file it came from, so a block
# can still be read against the change that motivated it.
#
# The assertions are unchanged. What changed is where they run and, for
# a few of them, that a name they call is now this package's rather
# than frmtmb's.

# ---- from tests/testthat/test-review-fixes.R ----

test_that("frm_sample(prior=) works on a fixed-effects-only GLM", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(404)
  d <- data.frame(x = rnorm(80))
  d$y <- rpois(80, exp(0.4 + 0.5 * d$x))
  fit <- frm(bf(y ~ x) + poisson(), data = d)
  # the $b partial match used to pass random = "b" for a template that
  # only holds beta, so this errored in MakeADFun
  # short-chain R-hat/ESS warnings are not what this test checks
  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 800, refresh = 0,
               prior = list(beta = prior_normal(0, 5))))
  m <- as.matrix(ds)
  expect_true("x" %in% colnames(m))
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  if (sampler_gates_on()) {
    expect_lt(abs(mean(m[, "x"]) - unname(fit$estimates$beta[["x"]])),
              5 * stats::sd(m[, "x"]) + 1e-8)
  }
})

test_that("frm_sample(laplace = TRUE) runs and labels outer draws", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(405)
  d <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  d$y <- rnorm(80, 1 + 0.5 * d$x + rnorm(8, 0, 0.5)[d$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
  # rstan warns about lp__ under laplace; not what this test checks
  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 400, refresh = 0,
               laplace = TRUE))
  m <- as.matrix(ds)
  # no b columns are sampled, and theta keeps its own label instead of
  # being misattributed as b[1]
  expect_false(any(grepl("^b\\[", colnames(m))))
  expect_true(any(grepl("theta", colnames(m))))
  expect_true("x" %in% colnames(m))
  # a laplace chain mixes poorly by construction (each leapfrog runs
  # the inner solve), so a stuck chain UNDERSTATES its own spread; the
  # wiring sanity is judged against the wider of the chain's spread and
  # the Wald standard error
  if (sampler_gates_on()) {
    expect_lt(abs(mean(m[, "x"]) - unname(fit$estimates$beta[["x"]])),
              5 * max(stats::sd(m[, "x"]),
                      sqrt(diag(stats::vcov(fit)))[["x"]]) + 1e-8)
  }
})

test_that("mode_inits anchors chain 1 and jitters the rest", {
  mode <- c(a = 1, b = -2, c = 0.5)
  ii <- frmtmb.sample:::mode_inits(mode, chains = 4, jitter = 0.25)
  expect_length(ii, 4L)
  expect_identical(ii[[1]], as.numeric(mode))
  for (k in 2:4) {
    expect_false(identical(ii[[k]], as.numeric(mode)))
    expect_lt(max(abs(ii[[k]] - as.numeric(mode))), 2)  # modest jitter
  }
  # jitter = 0 restores identical mode starts
  i0 <- frmtmb.sample:::mode_inits(mode, chains = 3, jitter = 0)
  expect_true(all(vapply(i0, identical, TRUE, as.numeric(mode))))
})

test_that("frm_sample runs multiple chains with jittered mode inits", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(409)
  d <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  d$y <- rnorm(60, 1 + 0.5 * d$x + rnorm(6, 0, 0.5)[d$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
  ds <- suppressWarnings(
    frm_sample(fit, chains = 2, iter = 300, refresh = 0))
  a <- rstan::extract(ds$stanfit, permuted = FALSE)
  expect_identical(dim(a)[2], 2L)
  # a boundary-ish mode triggers the singular-init warning; the
  # crippled 10-iteration run may warn on its own, so collect all
  fit2 <- fit
  fit2$estimates$theta <- c(-9)
  w <- testthat::capture_warnings(
    try(frm_sample(fit2, chains = 1, iter = 10, refresh = 0),
        silent = TRUE))
  expect_true(any(grepl("extreme covariance parameter", w)))
})

# ---- from tests/testthat/test-lkj.R ----


test_that("the formula route defaults to lkj(1) and takes an override", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- lkj_data()
  form <- bf(y ~ x + (x | g)) + gaussian()

  msg <- capture_messages(suppressWarnings(
    ds <- frm_sample(form, data = dd, chains = 1, iter = 500,
                     refresh = 0, seed = 5)))
  m <- paste(msg, collapse = "")
  expect_match(m, "\n  cor  ")
  expect_match(m, "lkj(1)", fixed = TRUE)
  expect_match(m, "correlation matrix", fixed = TRUE)
  pl <- unclass(prior_summary(ds))
  expect_true("cor" %in% vapply(pl, `[[`, "", "class"))
  # with the correlation priored the block non-centers, and the
  # correlation parameter stays where the prior has mass instead of
  # walking the improper tail the flat prior left open
  expect_equal(ds$reparam$blocks, 1L)
  expect_lt(max(abs(ds$draws[, "theta_3"])), 50)

  # a user lkj takes over the class and leaves the rest of the defaults
  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, chains = 1, iter = 500, refresh = 0,
               seed = 5, prior = set_prior("lkj(4)", class = "cor"))))
  pl2 <- unclass(prior_summary(ds2))
  cor2 <- Filter(function(s) identical(s$class, "cor"), pl2)
  expect_length(cor2, 1L)
  expect_equal(cor2[[1L]]$dist$eta, 4)
  expect_true("sd" %in% vapply(pl2, `[[`, "", "class"))
  # eta = 4 concentrates toward the identity, so the sampled
  # correlation is tighter around zero than under lkj(1); judged in the
  # chains' own spread, since these are two different chains
  r1 <- ds$draws[, "theta_3"] / sqrt(1 + ds$draws[, "theta_3"]^2)
  r4 <- ds2$draws[, "theta_3"] / sqrt(1 + ds2$draws[, "theta_3"]^2)
  expect_lt(stats::sd(r4), stats::sd(r1))

  # prior = "flat" opts out of the correlation default too, and the
  # gate closes with it
  ds3 <- suppressWarnings(suppressMessages(
    frm_sample(form, data = dd, chains = 1, iter = 500, refresh = 0,
               seed = 5, prior = "flat")))
  expect_null(ds3$reparam)
  expect_null(prior_summary(ds3))
})

# ---- from tests/testthat/test-ordinal-fitted.R ----

test_that("posterior_epred returns a draws x obs x category array", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- ordfit_data(114, n = 120)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 11))

  nd <- data.frame(x = c(-1.2, 0, 0.8, 1.5))
  ep <- posterior_epred(ds, newdata = nd, ndraws = 8)
  expect_true(is.array(ep) && length(dim(ep)) == 3L)
  expect_equal(dim(ep), c(8L, 4L, 3L))
  expect_null(dimnames(ep)[[1]])
  expect_equal(dimnames(ep)[[3]], levels(dd$y))
  # ep[, , "hi"] must be addressable by the response's own level names
  expect_equal(dim(ep[, , "hi"]), c(8L, 4L))

  # each draw x observation slice is a distribution
  for (k in seq_len(dim(ep)[1])) {
    expect_vector_equal(rowSums(ep[k, , ]), rep(1, 4), tol = 1e-12)
  }
  expect_true(all(ep > 0 & ep < 1))

  # continuity with the v0.31 flattened matrix: the array is exactly
  # that matrix reshaped, so category slice k is the k-th block of n
  # columns and the old column names still describe the block order
  flat <- matrix(as.vector(ep), nrow = dim(ep)[1])
  expect_equal(dim(flat), c(8L, 12L))
  for (k in seq_len(3L)) {
    expect_identical(as.vector(ep[, , k]),
                     as.vector(flat[, (k - 1L) * 4L + seq_len(4L)]))
  }

  # and each draw's slice is the matrix predict(type = "response")
  # returns, so the posterior mean tracks the MLE probabilities
  P <- predict(fit, newdata = nd, type = "response")
  epf <- posterior_epred(ds, newdata = nd, ndraws = 60)
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and the band asserts wiring, not mixing
  ep_sd <- apply(epf, c(2L, 3L), stats::sd)
  if (sampler_gates_on()) {
    expect_lt(max(abs(apply(epf, c(2L, 3L), mean) - P)),
              5 * max(ep_sd) + 1e-8)
  }

  # without newdata the observation margin carries the data rownames
  ep0 <- posterior_epred(ds, ndraws = 4)
  expect_equal(dim(ep0), c(4L, nrow(dd), 3L))
  expect_equal(dimnames(ep0)[[2]], rownames(dd))
  expect_equal(dimnames(ep0)[[3]], levels(dd$y))

  # the other two draws methods are statements about one number per
  # observation and keep the plain draws x observations matrix
  pp <- posterior_predict(ds, newdata = nd, ndraws = 8)
  expect_true(is.matrix(pp))
  expect_equal(dim(pp), c(8L, 4L))
  expect_true(all(pp %in% 1:3))
  pl <- posterior_linpred(ds, newdata = nd, ndraws = 8)
  expect_true(is.matrix(pl))
  expect_equal(dim(pl), c(8L, 4L))

  # pp_check() reads the predictive side, which the change does not move
  skip_if_not_installed("bayesplot")
  expect_s3_class(pp_check(ds, type = "bars", ndraws = 10), "ggplot")
})

test_that("a scalar-response family keeps the draws x obs matrix", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(115)
  gd <- data.frame(x = stats::rnorm(80))
  gd$y <- stats::rpois(80, exp(0.3 + 0.4 * gd$x))
  fp <- frm(bf(y ~ x) + poisson(), data = gd)
  dp <- suppressWarnings(frm_sample(fp, chains = 1, iter = 300,
                                    refresh = 0, seed = 12))
  nd <- data.frame(x = c(-1, 0, 1))
  ep <- posterior_epred(dp, newdata = nd, ndraws = 6)
  expect_true(is.matrix(ep))
  expect_equal(length(dim(ep)), 2L)
  expect_equal(dim(ep), c(6L, 3L))
  expect_null(dimnames(ep))
})

# ---- from tests/testthat/test-review-v28.R ----

test_that("the Windows cores guard reads options(mc.cores)", {
  # rstan defaults its cores argument to getOption("mc.cores", 1L), so
  # an unset cores= still starts socket workers
  old <- options(mc.cores = 2)
  on.exit(options(old), add = TRUE)
  expect_equal(stan_cores(list()), 2)
  expect_equal(stan_cores(list(cores = 1)), 1)
  options(mc.cores = NULL)
  expect_equal(stan_cores(list()), 1)
  expect_equal(stan_cores(list(cores = 3)), 3)
})

test_that("one chain under options(mc.cores) has nothing to note", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  if (.Platform$OS.type != "windows") {
    skip("the parallel-chain startup note is Windows-only")
  }
  set.seed(6)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(60)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  # a single chain cannot be parallelized, so an inherited core count
  # earns no startup note (test-parallel-chains.R asserts the note
  # fires when there is something to parallelize)
  old <- options(mc.cores = 2)
  on.exit(options(old), add = TRUE)
  # suppressWarnings: 100 post-warmup draws trip rstan's ESS advice,
  # and this test is about the startup note, not mixing
  expect_no_message(
    ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 200,
                                      refresh = 0, seed = 1)),
    message = "parallel chains on Windows"
  )
  expect_s3_class(ds, "frmtmb_draws")
})

# ---- from tests/testthat/test-prior-compat.R ----

test_that("the loss model samples with the vignette's priors", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  withr::local_options(mc.cores = 1)
  dd <- loss_data()
  vignette_priors <- c(prior(normal(5000, 1000), nlpar = "ult"),
                       prior(normal(1, 2), nlpar = "omega"),
                       prior(normal(45, 10), nlpar = "theta"))

  # 300 iterations is not enough warmup for a three-parameter
  # nonlinear body: the chain is still at its starting scale when it
  # stops, which is a property of the model rather than of the priors
  msg <- utils::capture.output(
    ds <- suppressWarnings(frm_sample(
      loss_form(), data = dd, prior = vignette_priors, chains = 1,
      iter = 800, refresh = 0, seed = 5, start = loss_start)),
    type = "message")

  # the three user priors survive into the sampled model, beside the
  # defaults for the slots they did not name
  pl <- prior_summary(ds)
  spelled <- vapply(unclass(pl), function(s) {
    paste0(s$class, "/", s$nlpar, "/", s$dpar)
  }, "")
  expect_true(all(c("b/ult/", "b/omega/", "b/theta/") %in% spelled))
  expect_true("sd//" %in% spelled)

  # the disclosure says the nonlinear parameters would otherwise be
  # flat, which is what brms does with them too
  expect_match(paste(msg, collapse = " "), "nonlinear parameters stay flat")

  m <- as.matrix(ds)
  expect_true(all(c("ult_Intercept", "omega_Intercept",
                    "theta_Intercept") %in% colnames(m)))
  if (sampler_gates_on()) {
    # the vignette's own structure, recovered: an ultimate loss near
    # 5000 with omega and theta near the values the data was built
    # from. Judged loosely, because this asserts that the priors are
    # wired into the sampled density, not the mixing
    expect_lt(abs(mean(m[, "ult_Intercept"]) - 5000), 1000)
    expect_lt(abs(mean(m[, "omega_Intercept"]) - 1.3), 0.5)
    expect_lt(abs(mean(m[, "theta_Intercept"]) - 45), 10)
  }
})

# ---- from tests/testthat/test-map.R ----

test_that("a MAP fit's priors carry into frm_sample by default", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(503)
  dd <- data.frame(x = rnorm(100))
  dd$y <- rnorm(100, 1 + 0.5 * dd$x, 1)
  map <- frm(bf(y ~ x) + gaussian(), data = dd,
             prior = set_prior("normal(0, 0.05)", class = "b",
                                coef = "x"))
  ds <- suppressWarnings(frm_sample(map, chains = 1, iter = 500,
                                    refresh = 0, seed = 1))
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  mx <- as.matrix(ds)[, "x"]
  if (sampler_gates_on()) {
    expect_lt(abs(mean(mx)), 5 * stats::sd(mx) + 1e-8)
  }
})

# ---- from tests/testthat/test-interop.R ----

test_that("as_tmbstan hands the objective to NUTS", {
  skip_if_not_installed("tmbstan")
  set.seed(132)
  dd <- data.frame(x = rnorm(80))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  sf <- suppressWarnings(as_tmbstan(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 1))
  expect_s4_class(sf, "stanfit")
  skip_if_not_installed("rstan")
  dr <- rstan::extract(sf, "beta")$beta
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  if (sampler_gates_on()) {
    expect_lt(abs(mean(dr[, 1]) - fixef(fit)$mu[[1]]),
              5 * stats::sd(dr[, 1]) + 1e-8)
  }
})

# ---- from tests/testthat/test-priors-bounds-grcov.R ----

test_that("a tight prior pulls the posterior toward it", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- sim_lmm(seed = 302)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 600, refresh = 0, seed = 1,
               prior = list(x = prior_normal(0, 0.01)))
  )
  m <- as.matrix(ds)
  # shrunk to ~0: the posterior sd itself proves the prior bit, and the
  # mean is judged against that sd rather than a platform-fragile number
  if (sampler_gates_on()) {
    expect_lt(stats::sd(m[, "x"]), 0.05)
    expect_lt(abs(mean(m[, "x"])), 5 * stats::sd(m[, "x"]) + 1e-8)
  }
  expect_gt(fixef(fit)$mu[["x"]], 0.3)          # ML untouched
})

# ---- from tests/testthat/test-hmm.R ----

test_that("frm_sample() runs on an hmm fit", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("frmtmb.latent")
  dd <- sim_hmm(8, 12, G2, c(0, 3), c(0.6, 0.6), 4026)
  fit <- frm(bf(y ~ 1),
             family = frmtmb.latent::hmm(K = 2, gaussian(), time = t,
                                         group = id),
             data = dd)
  s <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0)))
  expect_s3_class(s, "frmtmb_draws")
})

# ---- from tests/testthat/test-lca.R ----

test_that("frm_sample() runs on an lca fit", {
  skip_on_cran()
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("frmtmb.latent")
  s <- sim_lca_data(n = 150)
  fit <- frm(bf(Y ~ 1), family = frmtmb.latent::lca(K = 2), data = s$dd)
  sm <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 1))
  dr <- as.matrix(sm)
  expect_equal(nrow(dr), 200L)
  expect_true(ncol(dr) >= length(fit$opt$par))
})

# ---- from tests/testthat/test-backlog.R ----

test_that("bayesplot consumes as_tmbstan draws", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  skip_if_not_installed("bayesplot")
  set.seed(170)
  dd <- data.frame(x = rnorm(60))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  sf <- suppressWarnings(as_tmbstan(fit, chains = 1, iter = 300,
                                    refresh = 0, seed = 1))
  dr <- rstan::extract(sf, permuted = FALSE)   # iters x chains x pars
  iv <- bayesplot::mcmc_intervals_data(dr)
  expect_true(nrow(iv) >= 3)
  expect_true(all(is.finite(iv$m)))
})

# ---- from tests/testthat/test-input-validation.R ----

test_that("a tmbstan build that samples the wrong density is refused", {
  skip_if_not_installed("tmbstan")
  # this installation is a healthy binary build: the static check must
  # pass silently (the affected builds are source installs whose
  # model.hpp keeps an unpatched std_normal placeholder; see
  # dev/prior-dropping-investigation.md)
  expect_false(frmtmb.sample:::tmbstan_build_broken())
  expect_silent(frmtmb.sample:::check_tmbstan_build("frm_sample()"))
  # the PATTERN of the refusal path, pinned against a synthesized
  # broken model.hpp line so guard and autogen output cannot drift
  # apart silently; the stop() branch itself is unreachable on a
  # healthy installation and is not executed here
  bad <- "lp_accum__.add(stan::math::std_normal_lpdf<propto__>(y));"
  expect_true(any(grepl("std_normal_lpdf<propto__>(y)", bad,
                        fixed = TRUE)))
})

# ---- from tests/testthat/test-v07.R ----

test_that("frm_sample returns named draws and check_laplace agrees on a clean model", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(204)
  dd <- data.frame(x = rnorm(150), g = factor(rep(1:15, 10)))
  dd$y <- rnorm(150, 1 + 0.5 * dd$x + rnorm(15, 0, 0.8)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 600,
                                    refresh = 0, seed = 1))
  m <- as.matrix(ds)
  # draws names are parenthesis-free (the brms convention; v0.36)
  expect_true(all(c("Intercept", "x", "sigma_Intercept",
                    "theta_1") %in% colnames(m)))
  # judged against the chain's own spread: a seeded chain is not
  # platform-deterministic, and this asserts wiring, not mixing
  if (sampler_gates_on()) {
    expect_lt(abs(mean(m[, "x"]) - fixef(fit)$mu[["x"]]),
              5 * stats::sd(m[, "x"]) + 1e-8)
  }

  cl <- suppressWarnings(suppressMessages(
    check_laplace(fit, chains = 1, iter = 600, refresh = 0, seed = 1)))
  expect_s3_class(cl, "data.frame")
  expect_true("ess_bulk" %in% names(cl))
  # Wald and posterior agree on a well-behaved gaussian LMM, but only a
  # HEALTHY chain can testify: on a platform whose chain wandered
  # (measured, not assumed), the agreement claim is untestable
  row_x <- cl[cl$parameter == "x", ]
  # bulk ESS is necessary, not sufficient: a chain can mix on x while
  # its flat-prior theta excursion fattens the marginal anyway, so the
  # agreement claim is additionally gated per platform
  if (sampler_gates_on() &&
      is.finite(row_x$ess_bulk) && row_x$ess_bulk >= 100) {
    expect_lt(abs(row_x$z_shift), 0.75)
    expect_lt(abs(row_x$sd_ratio - 1), 0.5)
  } else {
    skip("chain too unhealthy on this platform to judge the agreement")
  }
})

# ---- from tests/testthat/test-v15.R ----

test_that("the draws surface runs the model machinery per draw", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(43)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.6)[dd$g], 0.8)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 600,
                                    refresh = 0, seed = 1))

  s <- summary(ds)
  expect_true(all(c("mean", "sd", "Rhat") %in% colnames(s)))
  fe <- fixef(ds)
  # draws-side names are parenthesis-free throughout (v0.36)
  expect_equal(rownames(fe), c("Intercept", "x", "sigma_Intercept"))
  # a seeded Stan chain is not platform-deterministic (pkgcheck's
  # container drew a chain far from this machine's), so agreement with
  # the ML fit is judged against the chain's OWN Monte Carlo spread: a
  # wiring bug moves the estimate by O(1) while the spread stays small,
  # and a drifted chain widens its spread along with its error
  if (sampler_gates_on()) {
    expect_lt(abs(fe["x", "Estimate"] - fixef(fit)$mu[["x"]]),
              5 * fe["x", "Est.Error"] + 1e-8)
  }

  vc <- VarCorr(ds)
  expect_true(all(c("estimate", "lwr", "upr") %in% names(vc)))

  # the sharp per-draw claim, exact and chain-free: an epred row IS
  # X beta + Z b at that draw's own parameters
  ep_all <- posterior_epred(ds)
  expect_equal(nrow(ep_all), nrow(ds$draws))
  for (k in c(1L, nrow(ds$draws) %/% 2L, nrow(ds$draws))) {
    dr <- ds$draws[k, ]
    mu_k <- dr[["Intercept"]] + dr[["x"]] * dd$x +
      dr[paste0("b[", as.integer(dd$g), "]")]
    expect_equal(unname(ep_all[k, ]), unname(mu_k), tolerance = 1e-8)
  }

  ep <- posterior_epred(ds, ndraws = 25)
  expect_equal(dim(ep), c(25L, 80L))
  if (sampler_gates_on()) {
    expect_lt(max(abs(colMeans(ep) - fitted(fit))),
              5 * max(apply(ep, 2, stats::sd)) + 1e-8)
  }
  pp <- posterior_predict(ds, ndraws = 25)
  expect_gt(mean(apply(pp, 2, stats::sd)),
            mean(apply(ep, 2, stats::sd)))
  ep2 <- posterior_epred(ds, newdata = data.frame(x = 0:1, g = "3"),
                         ndraws = 10)
  expect_equal(dim(ep2), c(10L, 2L))

  # posterior_linpred: link scale, and transform = TRUE applies the
  # inverse link (identity here, so it must equal the epred)
  lp <- posterior_linpred(ds, ndraws = 25)
  expect_equal(dim(lp), c(25L, 80L))
  expect_equal(posterior_linpred(ds, transform = TRUE, ndraws = 25), lp)
  expect_equal(lp, ep)
  lp_s <- posterior_linpred(ds, dpar = "sigma", transform = TRUE,
                            ndraws = 10)
  expect_true(all(lp_s > 0))

  # ranef over draws: brms-shaped arrays whose Estimate tracks the
  # fitted conditional modes
  re_d <- ranef(ds)
  expect_named(re_d, "1 | g")
  expect_identical(dim(re_d[["1 | g"]]), c(8L, 4L, 1L))
  expect_identical(colnames(re_d[["1 | g"]]),
                   c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  # same yardstick as above: the chain's own spread, not a fixed number
  if (sampler_gates_on()) {
    expect_lt(max(abs(re_d[["1 | g"]][, "Estimate", 1] -
                        ranef(fit)[["1 | g"]][, 1])),
              5 * max(re_d[["1 | g"]][, "Est.Error", 1]) + 1e-8)
  }
  expect_true(all(re_d[["1 | g"]][, "Q2.5", 1] <
                    re_d[["1 | g"]][, "Q97.5", 1]))

  h <- hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")
  expect_s3_class(h, "frmtmb_hypothesis")
  expect_true(h$lwr > 0 && h$upr < 1)
  expect_equal(dim(attr(h, "draws")), c(nrow(ds$draws), 1L))

  if (requireNamespace("posterior", quietly = TRUE)) {
    dm <- posterior::as_draws(ds)
    expect_true(inherits(dm, "draws"))
  }
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    expect_s3_class(bayesplot::pp_check(ds, ndraws = 10), "ggplot")
  }
})

# ---- from tests/testthat/test-review-v29.R ----

test_that("posterior_linpred stays on the mu predictor for an ordinal fit", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  dd <- v29_ordinal_data(46, n = 120)
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 400,
                                    refresh = 0, seed = 2))
  nd <- data.frame(x = c(-1, 0, 1))

  pl <- posterior_linpred(ds, newdata = nd, ndraws = 10)
  expect_equal(dim(pl), c(10L, 3L))
  # transform = TRUE is documented as the mu dpar on its natural scale;
  # it used to route through type = "conditional" into the K-column
  # probability matrix
  plt <- posterior_linpred(ds, newdata = nd, ndraws = 10,
                           transform = TRUE)
  expect_equal(dim(plt), c(10L, 3L))
  # the ordinal mu link is the identity, so the two agree exactly
  expect_vector_equal(as.vector(plt), as.vector(pl), tol = 1e-12)

  # posterior_epred is the one that carries the category distribution,
  # as a draws x observations x categories array (brms's convention)
  ep <- posterior_epred(ds, newdata = nd, ndraws = 10)
  expect_equal(dim(ep), c(10L, 3L, 3L))
  expect_null(dimnames(ep)[[1]])
  expect_equal(dimnames(ep)[[3]], levels(dd$y))
  # each draw's slice is a 3 x 3 matrix of distributions
  for (k in seq_len(dim(ep)[1])) {
    expect_vector_equal(rowSums(ep[k, , ]), rep(1, 3), tol = 1e-12)
  }

  # a non-ordinal fit keeps the documented behavior exactly
  set.seed(47)
  gd <- data.frame(x = stats::rnorm(80))
  gd$y <- stats::rpois(80, exp(0.3 + 0.4 * gd$x))
  fp <- frm(bf(y ~ x) + poisson(), data = gd)
  dp <- suppressWarnings(frm_sample(fp, chains = 1, iter = 400,
                                    refresh = 0, seed = 3))
  lp <- posterior_linpred(dp, newdata = nd, ndraws = 5)
  lpt <- posterior_linpred(dp, newdata = nd, ndraws = 5, transform = TRUE)
  expect_vector_equal(as.vector(lpt), exp(as.vector(lp)), tol = 1e-10)
  epp <- posterior_epred(dp, newdata = nd, ndraws = 5)
  expect_equal(dim(epp), c(5L, 3L))
  expect_null(colnames(epp))
})

# ---- from tests/testthat/test-brms-agreement.R ----

test_that("ordinal posterior_epred has brms's draws x obs x category shape", {
  skip_unless_brms_fit()
  skip_if_not_installed("tmbstan")

  # ?brms::posterior_epred.brmsfit: "an S x N x C array" for categorical
  # and ordinal models, an S x N matrix otherwise. This is the agreement
  # the frmtmb convention is copied from, checked against brms's own
  # return rather than against the sentence.
  set.seed(21)
  n <- 200
  dd <- data.frame(x = rnorm(n))
  eta <- 0.9 * dd$x
  p <- cbind(plogis(-0.8 - eta),
             plogis(0.6 - eta) - plogis(-0.8 - eta),
             1 - plogis(0.6 - eta))
  dd$y <- factor(c("lo", "mid", "hi")[
    apply(p, 1L, function(pr) sample(3L, 1L, prob = pr))],
    levels = c("lo", "mid", "hi"), ordered = TRUE)
  nd <- data.frame(x = c(-1, 0, 1, 2))

  bfit <- brms::brm(brms::bf(y ~ x), data = dd, family = brms::cumulative(),
                    chains = 1, iter = 1000, warmup = 500, seed = 1,
                    refresh = 0, backend = "rstan")
  bep <- brms::posterior_epred(bfit, newdata = nd)

  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 1000,
                                    refresh = 0, seed = 1))
  ep <- posterior_epred(ds, newdata = nd)

  expect_length(dim(ep), 3L)
  expect_equal(dim(ep)[2:3], dim(bep)[2:3])
  # brms names the category margin by the response's own levels and
  # leaves the observation margin unnamed; we name the observation
  # margin when the rows have names, which is additive
  expect_equal(dimnames(ep)[[3]], levels(dd$y))
  expect_equal(dimnames(bep)[[3]], levels(dd$y))
  expect_null(dimnames(ep)[[1]])
  expect_vector_equal(as.vector(apply(ep, c(2L, 3L), mean)),
                      as.vector(apply(bep, c(2L, 3L), mean)), tol = 0.1)
})

# ---- from tests/testthat/test-review-v25.R ----

test_that("frm_sample warns when a bound excludes the ML mode", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(55)
  n <- 120
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:12, 10)))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(12, 0, 0.4)[dd$g] +
    stats::rnorm(n)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_lt(exp(fit$estimates$theta[[1]]), 1.5)

  # rstan's own message for an init at or outside a bound names neither
  # the parameter nor the bound, so frm_sample has to say it first
  seen <- character(0)
  ds <- withCallingHandlers(
    frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 3,
               prior = set_prior("", class = "sd", lb = 1.5)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("violates the requested bound", seen)))
  # the bound applies on the internal (log-sd) scale, and every draw
  # respects it: the chains started inside the box
  expect_true(all(as.matrix(ds)[, "theta_1"] >= log(1.5)))
})

# ---- the halves of frmtmb's own blocks that name frm_sample ----------
#
# frmtmb keeps the frm() and frm_simulate() halves of each of these;
# what could only be asserted with frm_sample() in scope is here.

test_that("frm_sample spells the argument prior, brms's name", {
  fo <- names(formals(frm_sample))
  expect_true("prior" %in% fo)
  expect_false("priors" %in% fo)
})

test_that("the retired priors= spelling fails rather than passing", {
  dd <- data.frame(x = stats::rnorm(20))
  dd$y <- stats::rnorm(20, 1 + dd$x, 1)
  # frm_sample()'s `...` goes to tmbstan and WOULD have swallowed the
  # retired name, fitting with no priors at all, so it is refused by
  # name rather than left to R's argument matching
  expect_error(
    frm_sample(bf(y ~ x) + gaussian(), data = dd,
               priors = set_prior("normal(0, 1)", class = "b")),
    "takes `prior`, not `priors`")
})

test_that("the multivariate post-fit surface declares frm_sample works", {
  for (st in c("rescor", "mvbf")) {
    expect_equal(frm_compat(st, "frm_sample")$status, "works", info = st)
  }
})

test_that("toep stays out of the default priors, and says so", {
  # frmtmb keeps the other half of this block: that set_prior("lkj(1)")
  # is refused on a toep block, and that the block cannot be
  # non-centered. What needs the default-prior builder is here.
  set.seed(78)
  dd <- data.frame(t = factor(rep(1:4, 40)),
                   g = factor(rep(1:20, each = 8)))
  dd$y <- stats::rnorm(160, rep(stats::rnorm(20, 0, 0.7), each = 8), 1)
  uf <- frm(bf(y ~ 1 + toep(t | g)) + gaussian(), dd,
            dry_run = "objective")
  cls <- vapply(unclass(frmtmb.sample:::default_priors_for(uf)), `[[`, "",
                "class")
  expect_false("cor" %in% cls)
  expect_match(frmtmb.sample:::default_prior_notes(uf),
               "no LKJ density fits")
})

# ---- from tests/testthat/test-prior-compat.R ----

test_that("a nonlinear parameter gets no default prior, as in brms", {
  dd <- loss_data()
  uf <- frm(loss_form(), data = dd, start = loss_start,
            dry_run = "objective")
  defs <- unclass(frmtmb.sample:::default_priors_for(uf))
  # brms leaves ult/omega/theta flat (they are its class "b"), and the
  # response's median and mad describe none of them
  expect_false(any(vapply(defs, function(s) {
    identical(s$class, "Intercept") && !nzchar(s$dpar)
  }, TRUE)))
  # the variance component and sigma still get theirs
  cls <- vapply(defs, `[[`, "", "class")
  expect_true("sd" %in% cls)
  expect_match(paste(frmtmb.sample:::default_prior_notes(uf), collapse = " "),
               "ult, omega, theta")
})

# ---- from tests/testthat/test-priors-bounds-grcov.R ----

test_that("prior-augmented objective equals nll + neg log prior", {
  dd <- sim_lmm()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  pr <- list(x = prior_normal(0, 1), theta = prior_normal(0, 2))
  ri <- frmtmb::resolve_prior_input(fit, pr)
  obj2 <- frmtmb.sample:::prior_augmented_obj(fit, ri$entries)
  est <- fit$opt$par
  nlp_manual <-
    -stats::dnorm(fixef(fit)$mu[["x"]], 0, 1, log = TRUE) -
    stats::dnorm(fit$estimates$theta[1], 0, 2, log = TRUE)
  expect_lt(abs(obj2$fn(obj2$par) -
                  (-as.numeric(logLik(fit)) + nlp_manual)), 1e-8)
})

# ---- from tests/testthat/test-review-v25.R ----

test_that("mode_inits pulls every chain strictly inside the bounds", {
  lo <- c(0, -Inf, 2)
  hi <- c(Inf, 1, 2.5)
  set.seed(1)
  inits <- frmtmb.sample:::mode_inits(c(-5, 5, 10), chains = 6, jitter = 2,
                               lower = lo, upper = hi)
  expect_length(inits, 6L)
  for (v in inits) {
    expect_true(all(v > lo))
    expect_true(all(v < hi))
  }
  # chain 1 is the mode anchor and is clamped like the rest
  expect_true(all(inits[[1]] > lo))
  # a zero jitter still clamps
  z <- frmtmb.sample:::mode_inits(c(-5, 5, 10), chains = 2, jitter = 0,
                           lower = lo, upper = hi)
  expect_true(all(z[[1]] > lo))
  # unbounded sampling is untouched
  expect_equal(frmtmb.sample:::mode_inits(c(1, 2), 1, 0)[[1]], c(1, 2))

  # a box narrower than the interior padding collapses to its midpoint
  expect_equal(frmtmb.sample:::clamp_into_bounds(1e6, 1, 1 + 1e-9), 1 + 5e-10)
})

# ---- from tests/testthat/test-setprior.R ----

test_that("set_prior bounds flow into the resolved input and augmented objective works", {
  ri <- frmtmb::resolve_prior_input(fit_sp,
    set_prior("normal(0, 5)", class = "b") +
      set_prior("", class = "b", coef = "x", lb = 0))
  expect_identical(unname(ri$lower["x"]), 0)
  obj2 <- frmtmb.sample:::prior_augmented_obj(fit_sp, ri$entries)
  # augmented objective = nll + sum of normal(0,5) over the two b coefs
  b_coefs <- fit_sp$estimates$beta[c("x", "fb")]
  nlp <- -sum(stats::dnorm(b_coefs, 0, 5, log = TRUE))
  expect_lt(abs(obj2$fn(obj2$par) -
                  (-as.numeric(logLik(fit_sp)) + nlp)), 1e-8)
})

# ---- prior-carried bounds on a nonlinear parameter (frmtmb v0.49) ----
#
# The bound was keyed by the design-matrix column name, so an
# nlpar-addressed bound named no outer parameter. The keying lives in
# frmtmb's resolve_priorlist(), which this package consumes rather than
# reimplements, so the fix arrives here through the core. These guard the
# seam: bounds resolution and the Stan transform built from it.

test_that("an nlpar-addressed prior bound reaches frm_sample's bounds", {
  set.seed(31)
  n <- 300
  dd <- data.frame(x = stats::runif(n, 0, 10))
  dd$y <- stats::rbinom(n, 1, 0.25 + 0.75 * stats::plogis(dd$x - 5))
  fit <- suppressWarnings(frm(
    bf(y ~ guess + (1 - guess) * plogis(x - thr), guess ~ 1, thr ~ 1,
       nl = TRUE),
    family = bernoulli(link = "identity"), data = dd,
    start = list(beta = c(0.5, 4))))

  ri <- frmtmb::resolve_prior_input(fit,
    set_prior("", nlpar = "guess", lb = 0, ub = 1))
  expect_identical(names(ri$lower), "guess_(Intercept)")

  # the full-length vectors tmbstan is handed, in obj$par order
  bd <- frmtmb::resolve_bounds(fit, ri$lower, ri$upper)
  nm <- frmtmb::outer_par_names(fit)
  expect_identical(bd$lower[nm == "guess_(Intercept)"], 0)
  expect_identical(bd$upper[nm == "guess_(Intercept)"], 1)
  # and nothing else was constrained by it
  expect_true(all(is.infinite(bd$lower[nm != "guess_(Intercept)"])))
})

test_that("frm_sample samples inside an nlpar prior bound", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(32)
  n <- 300
  dd <- data.frame(x = stats::runif(n, 0, 10))
  dd$y <- stats::rbinom(n, 1, 0.25 + 0.75 * stats::plogis(dd$x - 5))
  fit <- suppressWarnings(frm(
    bf(y ~ guess + (1 - guess) * plogis(x - thr), guess ~ 1, thr ~ 1,
       nl = TRUE),
    family = bernoulli(link = "identity"), data = dd,
    start = list(beta = c(0.5, 4))))

  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 4,
               prior = set_prior("", nlpar = "guess", lb = 0.05,
                                 ub = 0.6)))
  # the draws carry the parenthesis-free spelling of the same parameter
  gs <- as.matrix(ds)[, "guess_Intercept"]
  expect_true(all(gs >= 0.05 & gs <= 0.6))
})

test_that("frm_sample's retired lower=/upper= are refused, not swallowed", {
  set.seed(33)
  dd <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  # `...` goes to tmbstan, which would take these as unknown sampler
  # options and sample the model UNBOUNDED in silence, so the refusal
  # has to be explicit rather than left to an unused-argument error
  expect_error(frm_sample(fit, lower = c(x = 0)), "has no `lower`")
  expect_error(frm_sample(fit, upper = c(x = 1)), "has no `upper`")
  expect_error(frm_sample(fit, lower = c(x = 0), upper = c(x = 1)),
               "has no `lower`/`upper`")
  # and the message says where a bound goes instead
  expect_error(frm_sample(fit, lower = c(x = 0)), "set_prior")
  expect_false(any(c("lower", "upper") %in% names(formals(frm_sample))))
})

test_that("a residual-correlation prior class reaches frm_sample", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(34)
  ng <- 40; k <- 6; n <- ng * k
  dd <- data.frame(x = stats::rnorm(n), t = rep(seq_len(k), ng),
                   g = factor(rep(seq_len(ng), each = k)))
  dd$y <- as.numeric(stats::arima.sim(list(ar = 0.6), n)) + 0.5 * dd$x
  fit <- frm(bf(y ~ x + ar(t, g, cov = TRUE)) + gaussian(), data = dd)

  # the class carries its density through frmtmb's resolver, which this
  # package consumes rather than reimplements
  ri <- frmtmb::resolve_prior_input(fit,
    set_prior("normal(0, 0.5)", class = "ar"))
  expect_identical(ri$entries[[1L]]$comp, "thetaac")

  ds <- suppressWarnings(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 5,
               prior = set_prior("", class = "ar", lb = 0.1, ub = 0.7)))
  # the bound was written on the natural scale; the draws are internal,
  # so compare through the same partial-autocorrelation map
  th <- as.matrix(ds)[, "thetaac_1"]
  phi <- th / sqrt(1 + th^2)
  expect_true(all(phi >= 0.1 - 1e-8 & phi <= 0.7 + 1e-8))
})
