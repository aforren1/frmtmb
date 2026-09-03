# Non-centered sampling: frm_sample(reparameterize = TRUE).
#
# The claim under test is that the route changes the SAMPLER's
# coordinates and nothing else. So the tests are (a) the factor identity
# L L' = Sigma for every structure that declares one, (b) the exact
# per-draw back-transform b = L(theta) z against a hand recomputation,
# (c) a draws surface that cannot tell the two routes apart, and (d) the
# gate: a block is non-centered only when every parameter it has is a
# standard deviation AND carries a prior.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

skip_sampler <- function() {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
}

rp_data <- function(seed = 9, n = 60L, ng = 6L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.5)[dd$g], 1)
  dd
}

# a prior on every theta, which is the second half of the eligibility
# gate; the formula interface's defaults supply the same thing
theta_priors <- function(fit) {
  frmtmb:::resolve_prior_input(fit, list(theta = prior_normal(0, 1)))
}

#' The z the sampler actually explored: the stanfit holds the parameters
#' as Stan saw them, in the same order the draws matrix labels, so the
#' `b[i]` positions of the raw matrix are the z draws.
raw_stan_matrix <- function(ds) {
  a <- rstan::extract(ds$stanfit, permuted = FALSE)
  m <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
  colnames(m) <- colnames(ds$draws)
  m
}

## ---- the factor accessors --------------------------------------------

test_that("every declared Cholesky accessor factorizes its own vcov", {
  blk <- function(cs, d, extra = list()) {
    c(list(covstruct = cs, dim = d, n_levels = 4L,
           cnms = paste0("t", seq_len(d))), extra)
  }
  cases <- list(
    list(blk("us", 1L), 0.3),
    list(blk("us", 3L), c(-0.2, 0.4, 0.1, 0.7, -0.3, 0.5)),
    list(blk("diag", 3L), c(0.1, -0.4, 0.2)),
    list(blk("homdiag", 3L), 0.2),
    list(blk("smooth", 5L), -0.3),
    list(blk("cs", 4L), c(0.1, -0.2, 0.3, 0, 0.8)),
    list(blk("cs", 4L), c(0.1, -0.2, 0.3, 0, -1.2)),
    list(blk("homcs", 4L), c(0.1, 0.9)),
    list(blk("ar1", 5L), c(0.2, 1.3)),
    list(blk("ar1", 5L), c(0.2, -1.3)),
    list(blk("hetar1", 5L), c(0.1, 0.2, -0.1, 0, 0.3, -0.7)),
    list(blk("equalto", 3L,
             list(aux_A = diag(3) * 0.7 + 0.3)), numeric(0)),
    list(blk("hsgp", 6L,
             list(aux_omega = matrix(seq_len(6), 6, 1), gp_iso = TRUE)),
         c(0.1, log(0.3)))
  )
  for (cse in cases) {
    bk <- cse[[1L]]
    th <- cse[[2L]]
    L <- frmtmb:::ncp_block_chol(bk, th)
    V <- frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(th, bk)
    expect_equal(L %*% t(L), unname(as.matrix(V)),
                 tolerance = 1e-12,
                 info = paste(bk$covstruct, bk$dim))
    # lower triangular, so the map is orientation-fixed and invertible
    expect_true(all(abs(L[upper.tri(L)]) < 1e-14))
    # and the map round-trips
    z <- stats::rnorm(bk$dim * bk$n_levels)
    b <- frmtmb:::ncp_scale_b(bk, z, th)
    expect_equal(frmtmb:::ncp_unscale_b(bk, b, th), z, tolerance = 1e-10)
  }
})

test_that("gr(cov = ) factors the level-major Kronecker covariance", {
  set.seed(3)
  A <- crossprod(matrix(stats::rnorm(16), 4)) + diag(4)
  d <- 2L
  nl <- 4L
  bk <- list(covstruct = "gr_cov", dim = d, n_levels = nl, aux_A = A,
             cnms = c("a", "b"),
             aux_kron = frmtmb:::kron_cov_index(d, nl))
  th <- c(0.2, -0.3, 0.4)
  LS <- frmtmb:::ncp_block_chol(bk, th)
  LA <- frmtmb:::covstruct_registry$gr_cov$chol_A(bk)
  Lfull <- kronecker(LA, LS)
  S <- frmtmb:::covstruct_registry$gr_cov$vcov(th, bk)
  K <- matrix(as.vector(A)[bk$aux_kron$ia] * as.vector(S)[bk$aux_kron$is],
              d * nl, d * nl)
  expect_equal(Lfull %*% t(Lfull), K, tolerance = 1e-10)
  z <- stats::rnorm(d * nl)
  expect_equal(frmtmb:::ncp_scale_b(bk, z, th),
               as.vector(Lfull %*% z), tolerance = 1e-10)
})

## ---- the gate --------------------------------------------------------

test_that("a block is eligible when its factor consumes all its theta", {
  el <- function(cs, d = 2L) {
    frmtmb:::ncp_eligible(list(covstruct = cs, dim = d, n_levels = 3L))
  }
  # every parameter is a standard deviation: nothing to expose
  for (cs in c("diag", "homdiag", "smooth", "equalto", "hsgp")) {
    expect_true(el(cs), info = cs)
  }
  expect_true(el("us", 1L))
  expect_true(el("gr_cov", 1L))

  # a correlation parameter, with a factor and (since 0.39) an LKJ
  # density to bound it: eligible, and the prior half of the gate then
  # decides (see ncp_plan below)
  for (cs in c("us", "cs", "homcs", "ar1", "hetar1", "gr_cov")) {
    expect_true(el(cs, 2L), info = cs)
    expect_equal(frmtmb:::block_n_cor(list(covstruct = cs, dim = 2L)),
                 1L, info = cs)
  }
  # no factor at all
  for (cs in c("us_t", "diag_t", "car", "spde", "gr_prec", "gp", "ou",
               "exp", "gau", "mat", "toep", "homtoep", "rr")) {
    expect_false(el(cs), info = cs)
  }
  # a Student-t latent is refused by its distribution, not by its
  # covariance structure: the same `us` factor exists, but the density
  # it belongs to is a scale mixture
  expect_false(frmtmb:::ncp_eligible(
    list(covstruct = "us", dim = 1L, n_levels = 3L, dist_nu = 5)))

  expect_match(frmtmb:::ncp_reason(list(covstruct = "car", dim = 1L)),
               "sparse CAR")
  expect_match(frmtmb:::ncp_reason(list(covstruct = "toep", dim = 2L)),
               "positive definite")
  expect_match(frmtmb:::ncp_reason(
    list(covstruct = "us", dim = 1L, dist_nu = 5)), "Student-t")
})

test_that("a block whose sd has no prior stays centered", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  # the structure qualifies
  expect_true(frmtmb:::ncp_eligible(fit$frame$re_blocks[[1L]]))
  # but with no prior on its variance, non-centering would hand the
  # chain the flat tail at sd = 0
  bare <- frmtmb:::ncp_plan(fit, TRUE, FALSE, NULL)
  expect_length(bare$idx, 0L)
  expect_match(bare$centered, "flat prior")
  # with one, it qualifies
  withp <- frmtmb:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
  expect_equal(withp$idx, 1L)
  expect_length(withp$centered, 0L)
  # a prior on beta alone is not a prior on the variance
  bo <- frmtmb:::resolve_prior_input(fit, list(beta = prior_normal(0, 5)))
  expect_length(frmtmb:::ncp_plan(fit, TRUE, FALSE, bo$entries)$idx, 0L)
})

test_that("laplace = TRUE ignores the reparameterization", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ent <- theta_priors(fit)$entries
  # the random effects are integrated out there, so there is no funnel
  # left to remove and the transform would only cost inner Newton steps
  expect_length(frmtmb:::ncp_plan(fit, TRUE, TRUE, ent)$idx, 0L)
  expect_length(frmtmb:::ncp_plan(fit, TRUE, TRUE, ent)$centered, 0L)
  expect_length(frmtmb:::ncp_plan(fit, FALSE, FALSE, ent)$idx, 0L)
})

test_that("a model with no random effects is unaffected", {
  set.seed(15)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- stats::rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_length(frmtmb:::ncp_plan(fit, TRUE, FALSE, NULL)$idx, 0L)
  expect_length(frmtmb:::ncp_plan(fit, TRUE, FALSE, NULL)$centered, 0L)
})

## ---- the objective ---------------------------------------------------

test_that("the non-centered objective is the centered one plus the log-Jacobian", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + diag(x | g)), family = gaussian(), data = dd)
  bk <- fit$frame$re_blocks[[1L]]

  frc <- fit$frame
  frn <- frc
  frn$ncp_blocks <- 1L
  f_c <- frmtmb:::build_objective(frc)
  f_n <- frmtmb:::build_objective(frn)

  set.seed(11)
  p <- fit$estimates
  p$b <- stats::rnorm(length(p$b), 0, 0.4)
  # a well-conditioned theta, set rather than perturbed: the identity is
  # exact in real arithmetic, so a near-singular block would be testing
  # floating point rather than the transform
  p$theta <- c(0.3, -0.2)
  pz <- p
  pz$b <- frmtmb:::ncp_unscale_b(bk, p$b, p$theta[bk$theta_idx])

  L <- frmtmb:::ncp_block_chol(bk, p$theta[bk$theta_idx])
  jac <- bk$n_levels * sum(log(diag(L)))
  # -log density of b equals -log density of z minus the log-Jacobian of
  # b = L z, per level; anything else means the two tapes describe
  # different models
  expect_equal(f_c(p) - f_n(pz), jac, tolerance = 1e-10)

  # and the fitted frame is never touched: a fit cannot acquire a
  # non-centered tape by having been sampled
  expect_null(fit$frame$ncp_blocks)
})

test_that("the identity holds for a CORRELATED block too", {
  dd <- rp_data(n = 90L, ng = 9L)
  fit <- suppressWarnings(
    frm(bf(y ~ x + (x | g)), family = gaussian(), data = dd))
  bk <- fit$frame$re_blocks[[1L]]
  expect_equal(bk$covstruct, "us")
  expect_equal(length(bk$theta_idx), 3L)   # two sds and a correlation

  frc <- fit$frame
  frn <- frc
  frn$ncp_blocks <- 1L
  f_c <- frmtmb:::build_objective(frc)
  f_n <- frmtmb:::build_objective(frn)

  set.seed(21)
  p <- fit$estimates
  p$b <- stats::rnorm(length(p$b), 0, 0.4)
  # a well-conditioned block, set rather than perturbed: the correlation
  # theta is the third entry
  p$theta <- c(0.3, -0.2, 0.7)
  pz <- p
  pz$b <- frmtmb:::ncp_unscale_b(bk, p$b, p$theta[bk$theta_idx])

  L <- frmtmb:::ncp_block_chol(bk, p$theta[bk$theta_idx])
  jac <- bk$n_levels * sum(log(diag(L)))
  expect_equal(f_c(p) - f_n(pz), jac, tolerance = 1e-10)
  # the correlation really is in the factor: a different correlation
  # gives a different Jacobian
  p2 <- p
  p2$theta[3L] <- -1.5
  pz2 <- p2
  pz2$b <- frmtmb:::ncp_unscale_b(bk, p2$b, p2$theta[bk$theta_idx])
  L2 <- frmtmb:::ncp_block_chol(bk, p2$theta[bk$theta_idx])
  expect_gt(abs(sum(log(diag(L2))) - sum(log(diag(L)))), 1e-3)
  expect_equal(f_c(p2) - f_n(pz2), bk$n_levels * sum(log(diag(L2))),
               tolerance = 1e-10)
})

test_that("reparameterize = must be a single TRUE or FALSE", {
  skip_sampler()
  dd <- rp_data()
  expect_error(frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                          data = dd, reparameterize = "yes"),
               "must be TRUE or FALSE")
  expect_error(frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                          data = dd, reparameterize = c(TRUE, TRUE)),
               "must be TRUE or FALSE")
})

test_that("the mode-anchored init starts the chains at the ML mode", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + diag(x | g)), family = gaussian(), data = dd)
  bk <- fit$frame$re_blocks[[1L]]
  p0 <- frmtmb:::ncp_start_pars(fit, 1L)
  # z0 = L^-1 b_hat, so the transform sends the start back onto the mode
  expect_equal(frmtmb:::ncp_scale_b(bk, p0[["b"]][bk$b_idx],
                                    p0$theta[bk$theta_idx]),
               unname(fit$estimates[["b"]][bk$b_idx]), tolerance = 1e-10)
  expect_equal(p0$beta, fit$estimates$beta)
  expect_equal(p0$theta, fit$estimates$theta)
})

## ---- the back-transform, exactly -------------------------------------

test_that("each b draw is L(theta of that draw) times its own z", {
  skip_sampler()
  dd <- rp_data()
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + diag(x | g)), family = gaussian(), data = dd,
               chains = 1, iter = 400, refresh = 0, seed = 4)))
  expect_equal(ds$reparam$blocks, 1L)

  idx <- frmtmb:::draws_par_index(ds$fit)
  raw <- raw_stan_matrix(ds)
  bk <- ds$fit$frame$re_blocks[[1L]]
  bcol <- idx[["b"]][bk$b_idx]
  tcol <- idx[["theta"]][bk$theta_idx]

  # the outer parameters are the sampler's own, untransformed
  expect_equal(unname(ds$draws[, idx[["beta"]]]),
               unname(raw[, idx[["beta"]]]))
  expect_equal(unname(ds$draws[, tcol]), unname(raw[, tcol]))

  for (i in c(1L, nrow(ds$draws) %/% 2L, nrow(ds$draws))) {
    th <- unname(raw[i, tcol])
    z <- unname(raw[i, bcol])
    # rebuilt by hand from theta, not by calling the transform: the
    # factor of an uncorrelated block is diag(exp(theta)), level-major
    b <- z * rep(exp(th), times = bk$n_levels)
    expect_equal(unname(ds$draws[i, bcol]), b, tolerance = 1e-10)
    # and the z really was standard normal in the sampler's coordinates,
    # i.e. not equal to b
    expect_gt(max(abs(z - b)), 1e-6)
  }
})

test_that("a correlated block non-centers once its correlation is priored", {
  skip_sampler()
  dd <- rp_data(n = 90L, ng = 9L)
  # the slope variance is at the boundary on this data, which is not
  # what is under test here: only the block's SHAPE and its priors
  # decide the gate
  fit <- suppressWarnings(
    frm(bf(y ~ x + (x | g)), family = gaussian(), data = dd))
  # with no prior at all the block stays centered, and the reason names
  # the correlation rather than the standard deviation
  bare <- frmtmb:::ncp_plan(fit, TRUE, FALSE, NULL)
  expect_length(bare$idx, 0L)
  expect_match(bare$centered, "improper")
  # a prior on the standard deviations alone is not enough: the
  # correlation is still flat
  sd_only <- frmtmb:::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd"))
  expect_length(frmtmb:::ncp_plan(fit, TRUE, FALSE,
                                  sd_only$entries)$idx, 0L)
  # sd + cor covers every theta the block has, which is the gate
  both <- frmtmb:::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd") +
      set_prior("lkj(1)", class = "cor"))
  plan <- frmtmb:::ncp_plan(fit, TRUE, FALSE, both$entries)
  expect_equal(plan$idx, 1L)
  expect_length(plan$centered, 0L)

  # and the formula route, whose defaults now cover both, non-centers it
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
               chains = 1, iter = 400, refresh = 0, seed = 19)))
  expect_equal(ds$reparam$blocks, 1L)

  # the per-draw back-transform, by hand: b = L(sd, cor at THAT draw) z,
  # with L rebuilt from theta rather than by calling the transform
  idx <- frmtmb:::draws_par_index(ds$fit)
  raw <- raw_stan_matrix(ds)
  bk <- ds$fit$frame$re_blocks[[1L]]
  bcol <- idx[["b"]][bk$b_idx]
  tcol <- idx[["theta"]][bk$theta_idx]
  for (i in c(1L, nrow(ds$draws) %/% 2L, nrow(ds$draws))) {
    th <- unname(raw[i, tcol])
    Lr <- diag(2)
    Lr[2L, 1L] <- th[3L]
    Lr <- Lr / sqrt(rowSums(Lr * Lr))     # row-normalized correlation
    L <- Lr * exp(th[1:2])                # times the sds, by row
    z <- unname(raw[i, bcol])
    b <- as.vector(L %*% matrix(z, bk$dim, bk$n_levels))
    expect_equal(unname(ds$draws[i, bcol]), b, tolerance = 1e-10)
    expect_gt(max(abs(z - b)), 1e-6)
  }

  # the FIT route is a likelihood diagnostic, so its priors are flat and
  # the same block stays centered there
  expect_message(
    ds2 <- suppressWarnings(
      frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 19)),
    "sampling stays centered")
  expect_null(ds2$reparam)
  expect_equal(unname(ds2$draws), unname(raw_stan_matrix(ds2)))
})

test_that("a fit with flat priors stays centered, and a prior turns it on", {
  skip_sampler()
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  expect_message(
    ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 300,
                                      refresh = 0, seed = 23)),
    "flat prior")
  expect_null(ds$reparam)
  expect_equal(unname(ds$draws), unname(raw_stan_matrix(ds)))

  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 23,
               priors = list(theta_1 = prior_normal(0, 1)))))
  expect_equal(ds2$reparam$blocks, 1L)
  expect_gt(max(abs(unname(ds2$draws) - unname(raw_stan_matrix(ds2)))),
            1e-6)
})

## ---- the draws surface cannot tell the routes apart -------------------

rp_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      dd <- rp_data()
      mk <- function(ncp) {
        suppressWarnings(suppressMessages(
          frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
                     chains = 1, iter = 800, refresh = 0, seed = 7,
                     reparameterize = ncp)))
      }
      cache <<- list(dd = dd, cen = mk(FALSE), ncp = mk(TRUE))
    }
    cache
  }
})

test_that("the draws matrix has the same columns in the same order", {
  cs <- rp_case()
  expect_identical(colnames(cs$ncp$draws), colnames(cs$cen$draws))
  expect_identical(dim(cs$ncp$draws), dim(cs$cen$draws))
  expect_true(any(startsWith(colnames(cs$ncp$draws), "b[")))
  # no z is exposed anywhere on the object's own surface
  expect_false(any(grepl("^z", colnames(cs$ncp$draws))))
  expect_null(cs$cen$reparam)
  expect_equal(cs$ncp$reparam$blocks, 1L)
})

test_that("log_lik() row sums are the taped likelihood on non-centered draws", {
  cs <- rp_case()
  ll <- log_lik(cs$ncp)
  idx <- frmtmb:::draws_par_index(cs$ncp$fit)
  i <- 42L
  sh <- frmtmb:::draws_fit_at(cs$ncp, i, idx)
  pars <- sh$estimates
  # the CENTERED objective, because the draw is on the centered scale:
  # that is the whole point of back-transforming at extraction
  nll <- frmtmb:::build_objective(sh$frame)(pars)
  bk <- sh$frame$re_blocks[[1L]]
  re_prior <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll(
    pars[["b"]][bk$b_idx], pars$theta[bk$theta_idx], bk)
  expect_equal(sum(ll[i, ]), as.numeric(-nll) - as.numeric(re_prior),
               tolerance = 1e-10)
})

test_that("posterior_epred() on non-centered draws is X beta + Z b per draw", {
  cs <- rp_case()
  ep <- posterior_epred(cs$ncp)
  expect_equal(nrow(ep), nrow(cs$ncp$draws))
  for (k in c(1L, nrow(cs$ncp$draws) %/% 2L, nrow(cs$ncp$draws))) {
    dr <- cs$ncp$draws[k, ]
    mu_k <- dr[["Intercept"]] + dr[["x"]] * cs$dd$x +
      dr[paste0("b[", as.integer(cs$dd$g), "]")]
    expect_equal(unname(ep[k, ]), unname(mu_k), tolerance = 1e-8)
  }
})

test_that("the whole draws method surface runs on non-centered draws", {
  cs <- rp_case()
  ds <- cs$ncp
  expect_true(is.matrix(summary(ds)))
  expect_true(is.matrix(fixef(ds)))
  expect_s3_class(VarCorr(ds), "data.frame")
  expect_true(is.list(ranef(ds)))
  expect_s3_class(hypothesis(ds, "sd_g__Intercept > 0"), "data.frame")
  ce <- conditional_effects(ds, effects = "x", resolution = 10)
  expect_equal(nrow(ce[[1L]]), 10L)
  pp <- posterior_predict(ds, ndraws = 10)
  expect_equal(dim(pp), c(10L, nrow(cs$dd)))
  expect_true(is.matrix(as.matrix(ds)))
  expect_s3_class(prior_summary(ds), "frmtmb_priorlist")
})

test_that("the two routes agree within their own Monte Carlo spread", {
  cs <- rp_case()
  # a seeded Stan chain is not platform-deterministic and the two routes
  # are two different chains by construction, so agreement is judged
  # against the chains' OWN spread: a wiring bug moves an estimate by
  # O(1) while the spread stays small
  keep <- setdiff(colnames(cs$cen$draws),
                  c("lp__", grep("^b\\[", colnames(cs$cen$draws),
                                 value = TRUE)))
  for (nm in keep) {
    a <- cs$cen$draws[, nm]
    b <- cs$ncp$draws[, nm]
    n_eff <- 50   # deliberately pessimistic, so the gate is about bias
    mcse <- sqrt(stats::var(a) / n_eff + stats::var(b) / n_eff)
    expect_lt(abs(mean(a) - mean(b)), 6 * mcse + 1e-8)
    expect_lt(abs(stats::sd(a) - stats::sd(b)),
              0.5 * max(stats::sd(a), stats::sd(b)) + 1e-8)
  }
  # the group-level values agree too, which is what the back-transform
  # is for
  for (nm in grep("^b\\[", colnames(cs$cen$draws), value = TRUE)) {
    a <- cs$cen$draws[, nm]
    b <- cs$ncp$draws[, nm]
    expect_lt(abs(mean(a) - mean(b)),
              0.5 * max(stats::sd(a), stats::sd(b)) + 1e-8)
  }
})

## ---- the same, on a correlated block ---------------------------------

rp_cor_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(31)
      ng <- 12L
      dd <- data.frame(x = stats::rnorm(120),
                       g = factor(rep(seq_len(ng), length.out = 120)))
      u <- matrix(stats::rnorm(2 * ng), 2)
      dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + 0.7 * u[1, dd$g] +
                             0.4 * u[2, dd$g] * dd$x, 1)
      mk <- function(ncp) {
        suppressWarnings(suppressMessages(
          frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
                     chains = 1, iter = 1000, refresh = 0, seed = 33,
                     reparameterize = ncp)))
      }
      cache <<- list(dd = dd, cen = mk(FALSE), ncp = mk(TRUE))
    }
    cache
  }
})

test_that("a correlated non-centered run keeps the draws surface exactly", {
  cs <- rp_cor_case()
  expect_equal(cs$ncp$reparam$blocks, 1L)
  expect_null(cs$cen$reparam)
  expect_identical(colnames(cs$ncp$draws), colnames(cs$cen$draws))
  expect_identical(dim(cs$ncp$draws), dim(cs$cen$draws))
  expect_false(any(grepl("^z", colnames(cs$ncp$draws))))
  # the centered route's draws ARE the sampled parameters; the
  # non-centered one's are not
  expect_equal(unname(cs$cen$draws), unname(raw_stan_matrix(cs$cen)))
  expect_gt(max(abs(unname(cs$ncp$draws) -
                      unname(raw_stan_matrix(cs$ncp)))), 1e-6)
  expect_true(is.matrix(summary(cs$ncp)))
  expect_s3_class(VarCorr(cs$ncp), "data.frame")
})

test_that("log_lik() row sums hold on a correlated non-centered run", {
  cs <- rp_cor_case()
  ll <- log_lik(cs$ncp)
  idx <- frmtmb:::draws_par_index(cs$ncp$fit)
  i <- 42L
  sh <- frmtmb:::draws_fit_at(cs$ncp, i, idx)
  pars <- sh$estimates
  nll <- frmtmb:::build_objective(sh$frame)(pars)
  bk <- sh$frame$re_blocks[[1L]]
  re_prior <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll(
    pars[["b"]][bk$b_idx], pars$theta[bk$theta_idx], bk)
  expect_equal(sum(ll[i, ]), as.numeric(-nll) - as.numeric(re_prior),
               tolerance = 1e-10)
})

test_that("the two routes agree on the correlated model as well", {
  cs <- rp_cor_case()
  # two different chains by construction, so agreement is judged in the
  # chains' OWN Monte Carlo spread (the house pattern): a wiring bug
  # moves an estimate by O(1) while the spread stays small
  keep <- setdiff(colnames(cs$cen$draws),
                  c("lp__", grep("^b\\[", colnames(cs$cen$draws),
                                 value = TRUE)))
  for (nm in keep) {
    a <- cs$cen$draws[, nm]
    b <- cs$ncp$draws[, nm]
    n_eff <- 50   # deliberately pessimistic, so the gate is about bias
    mcse <- sqrt(stats::var(a) / n_eff + stats::var(b) / n_eff)
    expect_lt(abs(mean(a) - mean(b)), 6 * mcse + 1e-8, label = nm)
    expect_lt(abs(stats::sd(a) - stats::sd(b)),
              0.5 * max(stats::sd(a), stats::sd(b)) + 1e-8, label = nm)
  }
  for (nm in grep("^b\\[", colnames(cs$cen$draws), value = TRUE)) {
    a <- cs$cen$draws[, nm]
    b <- cs$ncp$draws[, nm]
    expect_lt(abs(mean(a) - mean(b)),
              0.5 * max(stats::sd(a), stats::sd(b)) + 1e-8)
  }
})

## ---- priors ride on the outer parameters -----------------------------

test_that("a prior on an sd applies identically under both routes", {
  skip_sampler()
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  pr <- list(theta_1 = prior_normal(0, 0.3))
  ri <- frmtmb:::resolve_prior_input(fit, pr)

  # the prior lane never sees b (or z): it addresses theta, which is the
  # same parameter on either route
  expect_true(all(vapply(ri$entries, function(e) e$comp, "") == "theta"))

  frc <- fit$frame
  frn <- frc
  frn$ncp_blocks <- 1L
  nlp <- frmtmb:::neg_log_prior_fn(ri$entries)
  bk <- frc$re_blocks[[1L]]
  p <- fit$estimates
  set.seed(2)
  p$b <- stats::rnorm(length(p$b), 0, 0.4)
  p$theta <- 0.35
  pz <- p
  pz$b <- frmtmb:::ncp_unscale_b(bk, p$b, p$theta[bk$theta_idx])
  # the prior term is the same NUMBER at the same model point
  expect_equal(nlp(p), nlp(pz))

  # and the augmented objectives still differ by the log-Jacobian alone
  fc <- frmtmb:::build_objective(frc)
  fn <- frmtmb:::build_objective(frn)
  L <- frmtmb:::ncp_block_chol(bk, p$theta[bk$theta_idx])
  jac <- bk$n_levels * sum(log(diag(L)))
  expect_equal((fc(p) + nlp(p)) - (fn(pz) + nlp(pz)), jac,
               tolerance = 1e-10)

  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 5,
               priors = pr)))
  expect_equal(ds$reparam$blocks, 1L)
  # the prior is doing its job on the sampled theta, non-centered or not
  tcol <- frmtmb:::draws_par_index(ds$fit)[["theta"]]
  expect_lt(stats::sd(ds$draws[, tcol]), 1)
})

test_that("the formula route keeps its default priors under the new default", {
  skip_sampler()
  dd <- rp_data()
  expect_message(
    ds <- suppressWarnings(
      frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
                 chains = 1, iter = 300, refresh = 0, seed = 6)),
    "default priors")
  cls <- vapply(unclass(prior_summary(ds)), function(s) s$class, "")
  expect_true(all(c("Intercept", "sd") %in% cls))
  expect_equal(ds$reparam$blocks, 1L)

  # priors = "flat" opts out of the defaults, and with them out the
  # gate closes: nothing left to make the sd tail integrable
  expect_message(
    ds2 <- suppressWarnings(
      frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
                 chains = 1, iter = 300, refresh = 0, seed = 6,
                 priors = "flat")),
    "flat prior")
  expect_null(ds2$reparam)
})

## ---- blocks with no non-centered form --------------------------------

test_that("a Student-t block samples centered and says so", {
  skip_sampler()
  set.seed(12)
  dd <- data.frame(x = stats::rnorm(60),
                   g = factor(rep(seq_len(6), 10)))
  dd$y <- stats::rnorm(60, 1 + 0.5 * dd$x + stats::rt(6, 5)[dd$g], 1)
  ff <- bf(y ~ x + (1 | gr(g, dist = "student")))
  fit <- frm(ff, family = gaussian(), data = dd)
  plan <- frmtmb:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
  expect_length(plan$idx, 0L)
  expect_match(plan$centered, "Student-t latent")

  expect_message(
    ds <- suppressWarnings(
      frm_sample(ff, family = gaussian(), data = dd, chains = 1,
                 iter = 300, refresh = 0, seed = 8)),
    "sampling stays centered")
  expect_null(ds$reparam)
  # centered means the draws ARE the sampled parameters
  expect_equal(unname(ds$draws), unname(raw_stan_matrix(ds)))
})

test_that("a mixed model non-centers only the block that qualifies", {
  skip_sampler()
  set.seed(13)
  dd <- data.frame(x = stats::rnorm(90),
                   g = factor(rep(seq_len(6), 15)),
                   h = factor(rep(seq_len(5), each = 18)))
  dd$y <- stats::rnorm(90, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.6)[dd$g] +
                         stats::rt(5, 5)[dd$h], 1)
  ff <- bf(y ~ x + (1 | g) + (1 | gr(h, dist = "student")))
  fit <- frm(ff, family = gaussian(), data = dd)
  cs <- vapply(fit$frame$re_blocks, function(bk) bk$covstruct, "")
  plan <- frmtmb:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
  expect_equal(plan$idx, which(cs == "us"))
  expect_length(plan$centered, 1L)
  expect_match(plan$centered, "Student-t latent")

  expect_message(
    ds <- suppressWarnings(
      frm_sample(ff, family = gaussian(), data = dd, chains = 1,
                 iter = 300, refresh = 0, seed = 14)),
    "these blocks")
  expect_equal(ds$reparam$blocks, which(cs == "us"))

  # the transformed block moved, the centered one did not
  raw <- raw_stan_matrix(ds)
  idx <- frmtmb:::draws_par_index(ds$fit)
  bk_n <- fit$frame$re_blocks[[which(cs == "us")]]
  bk_c <- fit$frame$re_blocks[[which(cs != "us")]]
  expect_gt(max(abs(ds$draws[, idx[["b"]][bk_n$b_idx]] -
                      raw[, idx[["b"]][bk_n$b_idx]])), 1e-6)
  expect_equal(unname(ds$draws[, idx[["b"]][bk_c$b_idx]]),
               unname(raw[, idx[["b"]][bk_c$b_idx]]))
})

test_that("reparameterize = FALSE reproduces the centered route exactly", {
  cs <- rp_case()
  expect_equal(unname(cs$cen$draws), unname(raw_stan_matrix(cs$cen)))
})

## ---- structures reached through the formula interface ----------------

test_that("smooth and diag blocks non-center through the sampler", {
  skip_sampler()
  set.seed(16)
  dd <- data.frame(x = stats::runif(80, -3, 3),
                   g = factor(rep(seq_len(8), 10)))
  dd$y <- stats::rnorm(80, sin(dd$x) + stats::rnorm(8, 0, 0.5)[dd$g], 0.5)
  ff <- bf(y ~ s(x, k = 6) + diag(1 + x | g))
  fit <- frm(ff, family = gaussian(), data = dd)
  cs <- vapply(fit$frame$re_blocks, function(bk) bk$covstruct, "")
  expect_true(all(cs %in% c("smooth", "diag")))
  plan <- frmtmb:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
  expect_equal(plan$idx, seq_along(cs))
  expect_length(plan$centered, 0L)

  ds <- suppressWarnings(suppressMessages(
    frm_sample(ff, family = gaussian(), data = dd, chains = 1,
               iter = 300, refresh = 0, seed = 17)))
  expect_equal(ds$reparam$blocks, seq_along(cs))
  idx <- frmtmb:::draws_par_index(ds$fit)
  raw <- raw_stan_matrix(ds)
  for (i in seq_along(cs)) {
    bk <- ds$fit$frame$re_blocks[[i]]
    bcol <- idx[["b"]][bk$b_idx]
    tcol <- idx[["theta"]][bk$theta_idx]
    k <- 3L
    L <- frmtmb:::ncp_block_chol(bk, unname(raw[k, tcol]))
    z <- unname(raw[k, bcol])
    b <- as.vector(L %*% matrix(z, bk$dim, bk$n_levels))
    expect_equal(unname(ds$draws[k, bcol]), b, tolerance = 1e-10,
                 info = cs[i])
  }
  # and the surface still works with two blocks of different structures
  expect_true(is.matrix(posterior_epred(ds, ndraws = 5)))
})
