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
  frmtmb::resolve_prior_input(fit, set_prior("lognormal(0, 1)",
                                             class = "sd"))
}

#' The z the sampler actually explored: the stanfit holds the parameters
#' as Stan saw them, in the same order the draws matrix labels, so the
#' `b[i]` positions of the raw matrix are the z draws.
raw_stan_matrix <- function(ds) {
  a <- rstan::extract(ds$stanfit, permuted = FALSE)
  m <- do.call(rbind, lapply(seq_len(dim(a)[2]), function(ch) a[, ch, ]))
  # the draws matrix also holds brms's sd_ and cor_ columns, derived
  # after sampling, which the stanfit does not
  dv <- frmtmb.sample:::draws_natural_cols(ds$fit)$derived$names
  colnames(m) <- setdiff(colnames(ds$draws), dv)
  # the stored draws hold an unmodeled sigma on its natural scale, as
  # brms does; the sampler's own values are its log
  frmtmb.sample:::draws_to_natural(m, ds$fit)
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
    V <- frmtmb:::cs_entry(bk$covstruct)$vcov(th, bk)
    expect_equal(L %*% t(L), unname(as.matrix(V)),
                 tolerance = 1e-12,
                 info = paste(bk$covstruct, bk$dim))
    # lower triangular, so the map is orientation-fixed and invertible
    expect_true(all(abs(L[upper.tri(L)]) < 1e-14))
    # and the map round-trips
    z <- stats::rnorm(bk$dim * bk$n_levels)
    b <- frmtmb::ncp_scale_b(bk, z, th)
    expect_equal(frmtmb::ncp_unscale_b(bk, b, th), z, tolerance = 1e-10)
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
  LA <- frmtmb:::cs_entry("gr_cov")$chol_A(bk)
  Lfull <- kronecker(LA, LS)
  S <- frmtmb:::cs_entry("gr_cov")$vcov(th, bk)
  K <- matrix(as.vector(A)[bk$aux_kron$ia] * as.vector(S)[bk$aux_kron$is],
              d * nl, d * nl)
  expect_equal(Lfull %*% t(Lfull), K, tolerance = 1e-10)
  z <- stats::rnorm(d * nl)
  expect_equal(frmtmb::ncp_scale_b(bk, z, th),
               as.vector(Lfull %*% z), tolerance = 1e-10)
})

## ---- the gate --------------------------------------------------------

test_that("a block is eligible when its factor consumes all its theta", {
  el <- function(cs, d = 2L) {
    frmtmb::ncp_eligible(list(covstruct = cs, dim = d, n_levels = 3L))
  }
  # every parameter is a standard deviation: nothing to expose
  for (cs in c("diag", "homdiag", "smooth", "equalto", "hsgp", "gp")) {
    expect_true(el(cs), info = cs)
  }
  expect_true(el("us", 1L))
  expect_true(el("gr_cov", 1L))

  # a correlation parameter, with a factor and (since 0.39) an LKJ
  # density to bound it: eligible, and the prior half of the gate then
  # decides (see ncp_plan below)
  for (cs in c("us", "cs", "homcs", "ar1", "hetar1", "gr_cov")) {
    expect_true(el(cs, 2L), info = cs)
    expect_equal(frmtmb::block_n_cor(list(covstruct = cs, dim = 2L)),
                 1L, info = cs)
  }
  # no factor at all. car() and spde() are frmtmb.spatial's, which
  # asks the same of them in its own suite
  for (cs in c("us_t", "diag_t", "gr_prec", "ou",
               "exp", "gau", "mat", "toep", "homtoep", "rr")) {
    expect_false(el(cs), info = cs)
  }
  # a Student-t latent is refused by its distribution, not by its
  # covariance structure: the same `us` factor exists, but the density
  # it belongs to is a scale mixture
  expect_false(frmtmb::ncp_eligible(
    list(covstruct = "us", dim = 1L, n_levels = 3L, dist_nu = 5)))

  expect_match(frmtmb.sample:::ncp_reason(list(covstruct = "car", dim = 1L)),
               "sparse CAR")
  expect_match(frmtmb.sample:::ncp_reason(list(covstruct = "toep", dim = 2L)),
               "positive definite")
  expect_match(frmtmb.sample:::ncp_reason(
    list(covstruct = "us", dim = 1L, dist_nu = 5)), "Student-t")
})

test_that("a block whose sd has no prior stays centered", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  # the structure qualifies
  expect_true(frmtmb::ncp_eligible(fit$frame$re_blocks[[1L]]))
  # but with no prior on its variance, non-centering would hand the
  # chain the flat tail at sd = 0
  bare <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, NULL)
  expect_length(bare$idx, 0L)
  expect_match(bare$centered, "carries no density")
  # with one, it qualifies
  withp <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
  expect_equal(withp$idx, 1L)
  expect_length(withp$centered, 0L)
  # a prior on beta alone is not a prior on the variance
  bo <- frmtmb::resolve_prior_input(fit, set_prior("normal(0, 5)",
                                                   class = "b"))
  expect_length(frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, bo$entries)$idx, 0L)
})

test_that("laplace = TRUE ignores the reparameterization", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ent <- theta_priors(fit)$entries
  # the random effects are integrated out there, so there is no funnel
  # left to remove and the transform would only cost inner Newton steps
  expect_length(frmtmb.sample:::ncp_plan(fit, TRUE, TRUE, ent)$idx, 0L)
  expect_length(frmtmb.sample:::ncp_plan(fit, TRUE, TRUE, ent)$centered, 0L)
  expect_length(frmtmb.sample:::ncp_plan(fit, FALSE, FALSE, ent)$idx, 0L)
})

test_that("a model with no random effects is unaffected", {
  set.seed(15)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- stats::rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_length(frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, NULL)$idx, 0L)
  expect_length(frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, NULL)$centered, 0L)
})

## ---- the objective ---------------------------------------------------

test_that("the non-centered objective is the centered one plus the log-Jacobian", {
  dd <- rp_data()
  fit <- frm(bf(y ~ x + diag(x | g)), family = gaussian(), data = dd)
  bk <- fit$frame$re_blocks[[1L]]

  frc <- fit$frame
  frn <- frc
  frn$ncp_blocks <- 1L
  f_c <- frmtmb::build_objective(frc)
  f_n <- frmtmb::build_objective(frn)

  set.seed(11)
  p <- fit$estimates
  p$b <- stats::rnorm(length(p$b), 0, 0.4)
  # a well-conditioned theta, set rather than perturbed: the identity is
  # exact in real arithmetic, so a near-singular block would be testing
  # floating point rather than the transform
  p$theta <- c(0.3, -0.2)
  pz <- p
  pz$b <- frmtmb::ncp_unscale_b(bk, p$b, p$theta[bk$theta_idx])

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
  f_c <- frmtmb::build_objective(frc)
  f_n <- frmtmb::build_objective(frn)

  set.seed(21)
  p <- fit$estimates
  p$b <- stats::rnorm(length(p$b), 0, 0.4)
  # a well-conditioned block, set rather than perturbed: the correlation
  # theta is the third entry
  p$theta <- c(0.3, -0.2, 0.7)
  pz <- p
  pz$b <- frmtmb::ncp_unscale_b(bk, p$b, p$theta[bk$theta_idx])

  L <- frmtmb:::ncp_block_chol(bk, p$theta[bk$theta_idx])
  jac <- bk$n_levels * sum(log(diag(L)))
  expect_equal(f_c(p) - f_n(pz), jac, tolerance = 1e-10)
  # the correlation really is in the factor: a different correlation
  # gives a different Jacobian
  p2 <- p
  p2$theta[3L] <- -1.5
  pz2 <- p2
  pz2$b <- frmtmb::ncp_unscale_b(bk, p2$b, p2$theta[bk$theta_idx])
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
  p0 <- frmtmb.sample:::ncp_start_pars(fit, 1L)
  # z0 = L^-1 b_hat, so the transform sends the start back onto the mode
  expect_equal(frmtmb::ncp_scale_b(bk, p0[["b"]][bk$b_idx],
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

  idx <- frmtmb.sample:::draws_par_index(ds$fit)
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
  bare <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, NULL)
  expect_length(bare$idx, 0L)
  expect_match(bare$centered, "improper")
  # a prior on the standard deviations alone is not enough: the
  # correlation is still flat
  sd_only <- frmtmb::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd"))
  expect_length(frmtmb.sample:::ncp_plan(fit, TRUE, FALSE,
                                  sd_only$entries)$idx, 0L)
  # sd + cor covers every theta the block has, which is the gate
  both <- frmtmb::resolve_prior_input(
    fit, set_prior("exponential(1)", class = "sd") +
      set_prior("lkj(1)", class = "cor"))
  plan <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, both$entries)
  expect_equal(plan$idx, 1L)
  expect_length(plan$centered, 0L)

  # and the formula route, whose defaults now cover both, non-centers it
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
               chains = 1, iter = 400, refresh = 0, seed = 19)))
  expect_equal(ds$reparam$blocks, 1L)

  # the per-draw back-transform, by hand: b = L(sd, cor at THAT draw) z,
  # with L rebuilt from theta rather than by calling the transform
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
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

  # the FIT route takes the same defaults, so it non-centers the same
  # block, reports the same plan, and back-transforms the same way
  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 19)))
  expect_equal(ds2$reparam$blocks, ds$reparam$blocks)
  expect_equal(ds2$reparam$labels, ds$reparam$labels)
  expect_equal(ds2$reparam$centered, ds$reparam$centered)
  raw2 <- raw_stan_matrix(ds2)
  for (i in c(1L, nrow(ds2$draws))) {
    th <- unname(raw2[i, tcol])
    Lr <- diag(2)
    Lr[2L, 1L] <- th[3L]
    Lr <- Lr / sqrt(rowSums(Lr * Lr))
    L <- Lr * exp(th[1:2])
    z <- unname(raw2[i, bcol])
    expect_equal(unname(ds2$draws[i, bcol]),
                 as.vector(L %*% matrix(z, bk$dim, bk$n_levels)),
                 tolerance = 1e-10)
  }

  # the opt-out is brms's (flat): flat on the sds and on the Cholesky
  # factor, which closes both tails, so the block non-centers there too
  ds3 <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 19,
               prior = "flat")))
  expect_equal(ds3$reparam$blocks, 1L)
})

test_that("a fit non-centers on its defaults and under prior = flat", {
  skip_sampler()
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  # the class "sd" default covers this block's only theta, so the gate
  # opens on the fit route exactly as it does from the formula
  ds0 <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 23)))
  expect_equal(ds0$reparam$blocks, 1L)
  expect_gt(max(abs(unname(ds0$draws) - unname(raw_stan_matrix(ds0)))),
            1e-6)

  # the opt-out is brms's (flat) on the sd, which vanishes at sd = 0 in
  # theta's coordinates, so the gate stays open
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 23,
               prior = "flat")))
  expect_equal(ds$reparam$blocks, 1L)
  expect_gt(max(abs(unname(ds$draws) - unname(raw_stan_matrix(ds)))),
            1e-6)

  # reparameterize = FALSE still runs, and still hands back centered
  # draws, defaults or not
  ds1 <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 23,
               reparameterize = FALSE)))
  expect_null(ds1$reparam)
  expect_equal(unname(ds1$draws), unname(raw_stan_matrix(ds1)))
  expect_equal(colnames(ds1$draws), colnames(ds0$draws))

  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 23,
               prior = set_prior("lognormal(0, 1)", class = "sd"))))
  expect_equal(ds2$reparam$blocks, 1L)
  expect_gt(max(abs(unname(ds2$draws) - unname(raw_stan_matrix(ds2)))),
            1e-6)
})

test_that("check_laplace() samples flat and centered whatever the default", {
  skip_sampler()
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  # observable without re-deriving draws: the announcement it makes is
  # the flat-prior centered one, not the default-prior banner
  msg <- capture_messages(suppressWarnings(
    cl <- check_laplace(fit, chains = 1, iter = 300, refresh = 0,
                        seed = 23)))
  expect_length(grep("default priors", msg), 0L)
  expect_length(grep("carries no density", msg), 1L)
  # and the argument that carries it is explicit rather than inherited
  expect_true(".diagnostic" %in% names(formals(frm_sample)))
  expect_match(paste(deparse(body(check_laplace)), collapse = " "),
               ".diagnostic = TRUE", fixed = TRUE)
  # a MAP fit is checked against ITS OWN objective. Its penalty is
  # taped into fit$obj, so the diagnostic route leaves that tape alone
  # rather than resolving the prior a second time; only the user-facing
  # prior = "flat" rebuilds without it, and says so
  mf <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            prior = set_prior("normal(0, 0.3)", class = "b"))
  msg2 <- capture_messages(suppressWarnings(
    check_laplace(mf, chains = 1, iter = 300, refresh = 0, seed = 23)))
  expect_length(grep("default priors", msg2), 0L)
  expect_length(grep("drops the prior", msg2), 0L)
  expect_s3_class(cl, "data.frame")

  dg <- suppressWarnings(suppressMessages(
    frm_sample(mf, chains = 1, iter = 300, refresh = 0, seed = 23,
               .diagnostic = TRUE)))
  msg3 <- capture_messages(suppressWarnings(
    fl <- frm_sample(mf, chains = 1, iter = 300, refresh = 0, seed = 23,
                     prior = "flat")))
  expect_length(grep("drops the prior", msg3), 1L)
  # same seed, same everything but the penalty: the draws must differ
  expect_gt(max(abs(unname(dg$draws) - unname(fl$draws))), 1e-8)
  expect_null(dg$reparam)
  # prior = "flat" is brms's (flat) on the sd, which non-centers
  expect_equal(fl$reparam$blocks, 1L)
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
  expect_true(any(startsWith(colnames(cs$ncp$draws), "r_g[")))
  # no z is exposed anywhere on the object's own surface
  expect_false(any(grepl("^z", colnames(cs$ncp$draws))))
  expect_null(cs$cen$reparam)
  expect_equal(cs$ncp$reparam$blocks, 1L)
})

test_that("log_lik() row sums are the taped likelihood on non-centered draws", {
  cs <- rp_case()
  ll <- log_lik(cs$ncp)
  idx <- frmtmb.sample:::draws_par_index(cs$ncp$fit)
  i <- 42L
  sh <- frmtmb.sample:::draws_fit_at(cs$ncp, i, idx)
  pars <- sh$estimates
  # the CENTERED objective, because the draw is on the centered scale:
  # that is the whole point of back-transforming at extraction
  nll <- frmtmb::build_objective(sh$frame)(pars)
  bk <- sh$frame$re_blocks[[1L]]
  re_prior <- frmtmb:::cs_entry(bk$covstruct)$nll(
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
    mu_k <- dr[["b_Intercept"]] + dr[["b_x"]] * cs$dd$x +
      dr[paste0("r_g[", levels(cs$dd$g)[as.integer(cs$dd$g)], ",Intercept]")]
    expect_equal(unname(ep[k, ]), unname(mu_k), tolerance = 1e-8)
  }
})

test_that("the whole draws method surface runs on non-centered draws", {
  cs <- rp_case()
  ds <- cs$ncp
  expect_true(is.matrix(summary(ds)))
  expect_true(is.matrix(fixef(ds)))
  expect_named(VarCorr(ds)$g, "sd")
  expect_true(is.list(ranef(ds)))
  hy <- hypothesis(ds, "sd_g__Intercept > 0", class = NULL)
  # brms's SHAPE under frmtmb's own class: a frmtmb object must not
  # answer is(x, "<brms class>") TRUE (rule 2)
  expect_s3_class(hy, "frmtmb_hypothesis")
  expect_false(inherits(hy, "brmshypothesis"))
  expect_named(hy, c("hypothesis", "samples", "prior_samples",
                     "class", "alpha"))
  ce <- conditional_effects(ds, effects = "x", resolution = 10)
  expect_equal(nrow(ce[[1L]]), 10L)
  pp <- posterior_predict(ds, ndraws = 10)
  expect_equal(dim(pp), c(10L, nrow(cs$dd)))
  expect_true(is.matrix(as.matrix(ds)))
  expect_s3_class(prior_summary(ds), "frmtmb_priorlist")
})

test_that("the two routes agree within their own Monte Carlo spread", {
  skip_if_not(sampler_gates_on(), "chain-agreement gates are off")
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
      # 30 groups and sds 0.8 and 0.6, so that the CENTERED chain mixes
      # and its Monte Carlo error can be estimated from the chain: on
      # the 12 groups and sds 0.7 and 0.4 used before, a centered chain
      # of 1000 stuck (ESS 1 to 2 on some seeds) and no gate on its
      # draws was sound. The effects correlate (0.6), so that a wrong
      # factor shows: with uncorrelated effects the factor is nearly
      # diagonal, and a transposed one (b = L' z) passed the gate on 11
      # of 12 seeds; at 0.6 it fails on 12 of 12 (review of lane
      # moparity, r1). Over 24 seeds on 0.69.0 and on this build the
      # smallest centered ESS was 45 and the largest |z| 3.6
      # (dev/moparity-rp-gate.R, dev/moparity-findings.md, punch
      # round 2)
      set.seed(31)
      ng <- 30L
      dd <- data.frame(x = stats::rnorm(240),
                       g = factor(rep(seq_len(ng), length.out = 240)))
      u <- matrix(stats::rnorm(2 * ng), 2)
      u[2, ] <- 0.6 * u[1, ] + 0.8 * u[2, ]
      dd$y <- stats::rnorm(240, 1 + 0.5 * dd$x + 0.8 * u[1, dd$g] +
                             0.6 * u[2, dd$g] * dd$x, 1)
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
  expect_true(is.matrix(VarCorr(cs$ncp)[[1L]]$sd))
})

test_that("log_lik() row sums hold on a correlated non-centered run", {
  cs <- rp_cor_case()
  ll <- log_lik(cs$ncp)
  idx <- frmtmb.sample:::draws_par_index(cs$ncp$fit)
  i <- 42L
  sh <- frmtmb.sample:::draws_fit_at(cs$ncp, i, idx)
  pars <- sh$estimates
  nll <- frmtmb::build_objective(sh$frame)(pars)
  bk <- sh$frame$re_blocks[[1L]]
  re_prior <- frmtmb:::cs_entry(bk$covstruct)$nll(
    pars[["b"]][bk$b_idx], pars$theta[bk$theta_idx], bk)
  expect_equal(sum(ll[i, ]), as.numeric(-nll) - as.numeric(re_prior),
               tolerance = 1e-10)
})

test_that("the two routes agree on the correlated model as well", {
  skip_if_not(sampler_gates_on(), "chain-agreement gates are off")
  skip_if_not_installed("posterior")
  cs <- rp_cor_case()
  # two different chains by construction, so each mean and sd is
  # compared in units of the two chains' own Monte Carlo standard
  # errors, from each chain's ESS. A wiring bug moves an estimate by a
  # sizable part of a posterior sd, which is about 9.5 such units here
  # at the median ESS. Over 24 sampler seeds on 0.69.0 and on this
  # build (132 statistics each), the largest |z| was 3.6
  # (dev/moparity-rp-gate.R); a normal z passes 5 with probability
  # 5.7e-7
  z_lim <- 5
  # posterior caps an ESS above N log10(N) (an antithetic chain), which
  # can only make the MCSE larger and the gate looser
  mcse <- function(f, v) allow_warnings(f(v), "ESS has been capped")
  for (nm in setdiff(colnames(cs$cen$draws), "lp__")) {
    a <- cs$cen$draws[, nm]
    b <- cs$ncp$draws[, nm]
    z_mean <- (mean(a) - mean(b)) /
      sqrt(mcse(posterior::mcse_mean, a)^2 +
             mcse(posterior::mcse_mean, b)^2)
    z_sd <- (stats::sd(a) - stats::sd(b)) /
      sqrt(mcse(posterior::mcse_sd, a)^2 + mcse(posterior::mcse_sd, b)^2)
    expect_lt(abs(z_mean), z_lim, label = paste(nm, "mean"))
    expect_lt(abs(z_sd), z_lim, label = paste(nm, "sd"))
  }
})

## ---- priors ride on the outer parameters -----------------------------

test_that("a prior on an sd applies identically under both routes", {
  skip_sampler()
  dd <- rp_data()
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  # lognormal(0, 0.3) on the sd is normal(0, 0.3) on theta_1, its log
  pr <- set_prior("lognormal(0, 0.3)", class = "sd")
  ri <- frmtmb::resolve_prior_input(fit, pr)

  # the prior lane never sees b (or z): it addresses theta, which is the
  # same parameter on either route
  expect_true(all(vapply(ri$entries, function(e) e$comp, "") == "theta"))

  frc <- fit$frame
  frn <- frc
  frn$ncp_blocks <- 1L
  nlp <- frmtmb::neg_log_prior_fn(ri$entries)
  bk <- frc$re_blocks[[1L]]
  p <- fit$estimates
  set.seed(2)
  p$b <- stats::rnorm(length(p$b), 0, 0.4)
  p$theta <- 0.35
  pz <- p
  pz$b <- frmtmb::ncp_unscale_b(bk, p$b, p$theta[bk$theta_idx])
  # the prior term is the same NUMBER at the same model point
  expect_equal(nlp(p), nlp(pz))

  # and the augmented objectives still differ by the log-Jacobian alone
  fc <- frmtmb::build_objective(frc)
  fn <- frmtmb::build_objective(frn)
  L <- frmtmb:::ncp_block_chol(bk, p$theta[bk$theta_idx])
  jac <- bk$n_levels * sum(log(diag(L)))
  expect_equal((fc(p) + nlp(p)) - (fn(pz) + nlp(pz)), jac,
               tolerance = 1e-10)

  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 400, refresh = 0, seed = 5,
               prior = pr)))
  expect_equal(ds$reparam$blocks, 1L)
  # the prior is doing its job on the sampled theta, non-centered or not
  tcol <- frmtmb.sample:::draws_par_index(ds$fit)[["theta"]]
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

  # prior = "flat" replaces the defaults with brms's (flat), and flat
  # on the sd itself still makes the sd tail integrable: the gate stays
  # open
  ds2 <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
               chains = 1, iter = 300, refresh = 0, seed = 6,
               prior = "flat")))
  expect_equal(ds2$reparam$blocks, 1L)
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
  plan <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
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
  plan <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
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
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
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
  plan <- frmtmb.sample:::ncp_plan(fit, TRUE, FALSE, theta_priors(fit)$entries)
  expect_equal(plan$idx, seq_along(cs))
  expect_length(plan$centered, 0L)

  ds <- suppressWarnings(suppressMessages(
    frm_sample(ff, family = gaussian(), data = dd, chains = 1,
               iter = 300, refresh = 0, seed = 17)))
  expect_equal(ds$reparam$blocks, seq_along(cs))
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
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

## ---- an exact gp() field (dev/moparity-findings.md, item 7) ----------

test_that("an exact gp() field non-centers as brms's sdgp * L zgp", {
  # On 0.69.0 (rellib-r7) a gp block had no factor and was sampled
  # centered: on lane surface's y ~ gp(x) the formula route ran 1286 of
  # 2000 transitions into the maximum tree depth, with R-hat 1.52 on
  # b_Intercept (dev/moparity-gp-ncp.R)
  set.seed(5)
  n <- 60
  d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
  d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
  bk <- fit$frame$re_blocks[[1L]]
  expect_identical(bk$covstruct, "gp")
  expect_true(frmtmb::ncp_eligible(bk))
  th <- fit$estimates$theta[bk$theta_idx]
  L <- frmtmb:::ncp_block_chol(bk, th)
  V <- frmtmb:::cs_entry("gp")$vcov(th, bk)
  expect_equal(L %*% t(L), unname(as.matrix(V)),
               tolerance = 1e-10 * max(abs(V)))
  z <- stats::rnorm(bk$dim)
  b <- frmtmb::ncp_scale_b(bk, z, th)
  expect_equal(frmtmb::ncp_unscale_b(bk, b, th), z,
               tolerance = 1e-8 * max(abs(z)))
  # the non-centered objective is the centered one plus log |L|
  frn <- fit$frame
  frn$ncp_blocks <- 1L
  f_c <- frmtmb::build_objective(fit$frame)
  f_n <- frmtmb::build_objective(frn)
  p <- fit$estimates
  pz <- p
  pz$b <- frmtmb::ncp_unscale_b(bk, p$b, th)
  jac <- sum(log(diag(L)))
  expect_equal(f_c(p) - f_n(pz), jac, tolerance = 1e-8 * abs(jac))
  # and on the tape, where the factor is taken by RTMB's chol()
  obj_n <- RTMB::MakeADFun(f_n, pz, silent = TRUE)
  expect_equal(obj_n$fn(obj_n$par), f_n(pz), tolerance = 1e-10 *
                 abs(f_n(pz)))
  g <- obj_n$gr(obj_n$par)
  k <- which(names(obj_n$par) == "theta")[2L]
  fd_at <- function(h) {
    e <- replace(numeric(length(obj_n$par)), k, h)
    (obj_n$fn(obj_n$par + e) - obj_n$fn(obj_n$par - e)) / (2 * h)
  }
  fd <- fd_at(1e-5)
  # the difference's own error, measured: a step half as long changes
  # it by about its truncation error, and the rounding of fn bounds the
  # rest
  fd_err <- abs(fd - fd_at(5e-6)) +
    .Machine$double.eps * abs(obj_n$fn(obj_n$par)) / 5e-6
  expect_lt(abs(g[k] - fd), 10 * fd_err)
})

## ---- an exact gp() field at the sampler's nugget ----------------------

gp_low_noise_data <- function(noise) {
  set.seed(20261102)
  n <- 30
  x <- round(stats::runif(n, 0, 6), 3)
  data.frame(x = x, y = sin(2 * x) + noise * stats::rnorm(n))
}

test_that("the non-centered exact gp() carries a nugget of 1e-12, relative", {
  skip_sampler()
  # frm() keeps 1e-6 because its centered field does not converge reliably
  # below it, and the non-centered field has no such limit. At 1e-6 the sampled
  # field took its variance from sigma, as frm() does: the median sigma read
  # 0.18 of brms's at true sigma 1e-3, in 2 chains of 1500 on 40 rows, 1996 s
  # (dev/reviews/2026-11-02-gpgap.md, M3).
  d <- gp_low_noise_data(1e-3)
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
  # frm() itself is left where it was
  bk0 <- fit$frame$re_blocks[[1L]]
  expect_null(bk0$gp_nugget)
  ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 40, refresh = 0,
                                    seed = 19))
  bk <- ds$fit$frame$re_blocks[[1L]]
  expect_identical(bk$gp_nugget, 1e-12)
  # the draws krige through it: at 1e-12 it still predicts at new positions
  ep <- posterior_epred(ds, newdata = data.frame(x = c(1.23456, 5.4321)))
  expect_false(anyNA(ep))
  th <- fit$estimates$theta[bk$theta_idx]
  v <- frmtmb:::cs_entry("gp")$vcov(th, bk)
  v0 <- frmtmb:::cs_entry("gp")$vcov(th, bk0)
  sd2 <- exp(2 * th[[1L]])
  expect_equal(diag(v), rep(sd2 * (1 + 1e-12), bk$dim),
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(diag(v0), rep(sd2 * (1 + 1e-6), bk$dim),
               tolerance = 1e-14, ignore_attr = TRUE)

  # What the run's tape maximizes is the model at that nugget, not frm()'s:
  # the closed form, whose covariance sigma conditions, reads each one. A floor
  # at 1e-6 would lose a quarter of the one-parameter 95 percent boundary at
  # least, and the fit's tape keeps it.
  fr <- ds$fit$frame
  fr$ncp_blocks <- 1L
  tp <- fr$par_template
  tp$b[] <- 0
  tape <- RTMB::MakeADFun(frmtmb::build_objective(fr), tp, random = "b",
                          map = fr$map, silent = TRUE)
  op <- suppressWarnings(stats::nlminb(
    tape$par, tape$fn, tape$gr,
    control = list(iter.max = 800, eval.max = 1200)))
  u <- sort(unique(d$x))
  j <- match(d$x, u)
  nll <- function(p, nu) {
    s2 <- exp(2 * p[2])
    k <- s2 * (exp(-outer(u, u, "-")^2 / (2 * exp(2 * p[3]))) +
                 diag(nu, length(u)))
    w <- k[j, j]
    diag(w) <- diag(w) + exp(2 * p[4])
    r <- chol(w)
    z <- backsolve(r, d$y - p[1], transpose = TRUE)
    0.5 * sum(z^2) + sum(log(diag(r))) + 0.5 * length(z) * log(2 * pi)
  }
  p0 <- c(mean(d$y), th[[1L]], th[[2L]], log(1e-3))
  at_new <- stats::nlminb(p0, nll, nu = 1e-12)$objective
  at_old <- stats::nlminb(p0, nll, nu = 1e-6)$objective
  q <- stats::qchisq(0.95, 1) / 8
  expect_gt(at_old - at_new, q)
  expect_lt(abs(op$objective - at_new), q)
  expect_gt(abs(op$objective - at_old), q)
})

test_that("the sampler's starting z comes from the fit's own nugget", {
  skip_sampler()
  # through the smaller nugget the mode's field has a z the size of the
  # near-null directions, which is no place to start a chain
  d <- gp_low_noise_data(0.3)
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
  fs <- gp_sample_nugget(fit, 1L)
  bk <- fs$frame$re_blocks[[1L]]
  z <- ncp_start_pars(fs, 1L)$b[bk$b_idx]
  expect_lt(sum(z^2), 3 * bk$dim)
  # and it rebuilds the fit's field through the fit's nugget, in the factor's
  # order
  th <- fit$estimates$theta[bk$theta_idx]
  bk1 <- bk
  bk1$gp_nugget <- NULL
  b <- as.numeric(frmtmb::ncp_scale_b(bk1, z, th))
  b0 <- as.numeric(fit$estimates$b[bk$b_idx])
  expect_lt(max(abs(b - b0)), 1e-6 * max(abs(b0)))
})

test_that("a kernel past the reach of LAPACK at 1e-12 still factors", {
  skip_sampler()
  # a draw the factor on the tape took must not lose a finished chain
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = gp_low_noise_data(0.3))
  fs <- gp_sample_nugget(fit, 1L)
  bk <- fs$frame$re_blocks[[1L]]
  # a length scale far above the range of x leaves a kernel of one rank
  th <- c(0, log(1e4))
  z <- rep(1, bk$dim)
  b <- ncp_scale_b(bk, z, th)
  expect_length(b, bk$dim)
  expect_false(anyNA(b))
  # the field is the constant the kernel allows, near sd * the first z's scale
  expect_lt(diff(range(b)), 1e-2 * abs(mean(b)))
})


test_that("laplace = TRUE reads the exact gp() through its factor too", {
  skip_sampler()
  # the route that mixes at low noise: the integral is exact for a gaussian
  # response. It read frm()'s nugget as well, and sdgp 1.49 where brms reads
  # 2.17 at true sigma 1e-3 on 15 rows (dev/gpgap-findings.md)
  d <- gp_low_noise_data(1e-3)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ gp(x)), family = gaussian(), data = d, laplace = TRUE,
               chains = 1, iter = 60, refresh = 0, seed = 19)))
  bk <- ds$fit$frame$re_blocks[[1L]]
  expect_identical(bk$gp_nugget, 1e-12)
  expect_identical(ds$reparam$blocks, 1L)
  # no field in the draws to map back, and every outer parameter finite
  expect_true(draws_is_laplace(ds))
  m <- as.matrix(ds$draws)
  expect_false(any(startsWith(colnames(m), "b[")))
  expect_false(anyNA(m))
  # another block stays centered on this route
  d$g <- factor(rep(1:5, length.out = nrow(d)))
  d$y <- 0.3 + stats::rnorm(5)[as.integer(d$g)] + stats::rnorm(nrow(d))
  dr <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ 1 + (1 | g)), family = gaussian(), data = d,
               laplace = TRUE, chains = 1, iter = 60, refresh = 0, seed = 19)))
  expect_null(dr$reparam)
})


test_that("the factor is taken in a farthest-first order of the positions", {
  skip_sampler()
  # the sort was what kept the drawn field from mixing: rounding in a factor
  # over sorted positions moves the field by 1.7e-5 at n 15, against 6.1e-11 in
  # this
  # order, and at noise 1e-3 that is a gradient error of 34 on a gradient of 100
  blk <- list(aux_D2 = list(outer(0:4, 0:4, "-")^2))
  p <- gp_farthest_first(blk)
  # the medoid, then an end, then the other, then what is left
  expect_identical(p[1:3], c(3L, 1L, 5L))
  expect_identical(sort(p), 1:5)
  # rows at one position take their turn at the end
  x <- c(0, 0, 2, 2, 1)
  tied <- gp_farthest_first(list(aux_D2 = list(outer(x, x, "-")^2)))
  expect_identical(sort(tied), 1:5)
  # every dimension counts
  a <- c(0, 0, 1)
  b <- c(0, 5, 0)
  two <- gp_farthest_first(list(aux_D2 = list(outer(a, a, "-")^2,
                                              outer(b, b, "-")^2)))
  expect_identical(two[1:2], c(1L, 2L))

  d <- gp_low_noise_data(0.3)
  fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
  fs <- gp_sample_nugget(fit, 1L)
  bk <- fs$frame$re_blocks[[1L]]
  ord <- bk$gp_order
  expect_identical(sort(ord), seq_len(bk$dim))
  # frm()'s own block keeps none, and its factor stays a triangle
  expect_null(fit$frame$re_blocks[[1L]]$gp_order)
  th <- fit$estimates$theta[bk$theta_idx]
  L <- frmtmb:::ncp_block_chol(bk, th)
  V <- frmtmb:::cs_entry("gp")$vcov(th, bk)
  sd2 <- exp(2 * th[[1L]])
  # L L' is the kernel at its own nugget, in the block's order, up to the
  # rounding the factor is taken with
  expect_lt(max(abs(tcrossprod(L) - unname(as.matrix(V)))),
            1e-13 * sd2 * bk$dim)
  # a row-permuted triangle, with its triangle in the factor's order
  lp <- L[ord, ]
  expect_lt(max(abs(lp[upper.tri(lp)])), 1e-12 * sqrt(sd2))
  # the map and its inverse
  z <- stats::rnorm(bk$dim)
  bv <- frmtmb::ncp_scale_b(bk, z, th)
  expect_lt(max(abs(frmtmb::ncp_unscale_b(bk, bv, th) - z)), 1e-3 * max(abs(z)))
  # the identity with the centered objective: log |L| is a determinant now, not
  # the sum of a diagonal, at a nugget the centered density can still read
  fr1 <- fs$frame
  fr1$re_blocks[[1L]]$gp_nugget <- 1e-6
  fr2 <- fr1
  fr2$ncp_blocks <- 1L
  f_c <- frmtmb::build_objective(fr1)
  f_n <- frmtmb::build_objective(fr2)
  L1 <- frmtmb:::ncp_block_chol(fr1$re_blocks[[1L]], th)
  p <- fit$estimates
  pz <- p
  pz$b <- as.numeric(solve(L1, p$b))
  jac <- as.numeric(determinant(L1, logarithm = TRUE)$modulus)
  expect_equal(f_c(p) - f_n(pz), jac, tolerance = 1e-8 * abs(jac))
  # the tape reads what the numeric route reads
  tp <- RTMB::MakeADFun(frmtmb::build_objective(fr2), pz, silent = TRUE)
  expect_equal(tp$fn(tp$par), f_n(pz), tolerance = 1e-10 * abs(f_n(pz)))
})

test_that("each kind of exact block carries its own order", {
  skip_sampler()
  set.seed(20261103)
  n <- 40
  d <- data.frame(x = round(stats::runif(n, 0, 6), 3),
                  z = round(stats::runif(n, 0, 6), 3),
                  f = factor(sample(c("a", "b"), n, TRUE)),
                  t = sample(seq(0, 6, by = 0.5), n, TRUE))
  d$y <- sin(2 * d$x) + 0.5 * cos(d$z) + 0.3 * stats::rnorm(n)
  d$yt <- sin(d$t) + 0.3 * stats::rnorm(n)
  specs <- list(bf(y ~ gp(x) + gp(z)), bf(y ~ gp(x, by = f)),
                bf(y ~ gp(x, z, iso = FALSE)),
                # gr = FALSE keeps one value per row, so rows tie in the kernel
                bf(yt ~ gp(t, gr = FALSE)))
  for (sp in specs) {
    fit <- suppressWarnings(frm(sp, family = gaussian(), data = d))
    ib <- which(vapply(fit$frame$re_blocks, function(b) {
      identical(b$covstruct, "gp")
    }, NA))
    fs <- gp_sample_nugget(fit, ib)
    for (i in ib) {
      bk <- fs$frame$re_blocks[[i]]
      expect_identical(sort(bk$gp_order), seq_len(bk$dim))
      th <- fit$estimates$theta[bk$theta_idx]
      expect_false(anyNA(frmtmb::ncp_scale_b(bk, rep(1, bk$dim), th)))
    }
  }
  # the draws' names and the field's order do not move: a field in sorted order,
  # which reads at the observed rows as frm()'s does
  dd <- d[, c("x", "y")]
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ gp(x)), family = gaussian(), data = dd, chains = 1,
               iter = 200, refresh = 0, seed = 19)))
  m <- as.matrix(ds$draws)
  bk <- ds$fit$frame$re_blocks[[1L]]
  expect_identical(grep("^b[[]", colnames(m), value = TRUE),
                   paste0("b[", seq_len(bk$dim), "]"))
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = dd))
  expect_gt(stats::cor(colMeans(m[, startsWith(colnames(m), "b[")]),
                       as.numeric(fit$estimates$b)), 0.9)
  at <- posterior_epred(ds, newdata = dd[1:3, "x", drop = FALSE])
  expect_lt(max(abs(at - posterior_epred(ds)[, 1:3])), 1e-8)
})

test_that("a fit whose nugget the run moves starts at random, and says so", {
  skip_sampler()
  # from frm()'s collapsed mode the chains stayed there: sigma 6.9e-8 against
  # 1.03e-3, step size 1e-15, 1580 divergent (review of 2026-11-02, F2)
  d <- gp_low_noise_data(1e-3)
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
  s0 <- as.numeric(stats::sigma(fit))
  expect_lt(s0, 1e-5)
  said <- character()
  ds <- withCallingHandlers(
    suppressWarnings(frm_sample(fit, chains = 1, iter = 300, refresh = 0,
                                seed = 19)),
    message = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("chains start at random", said, fixed = TRUE)))
  expect_gt(stats::median(ds$draws[, "sigma"]), 1e2 * s0)
  # asked for, the mode is what it gets
  said <- character()
  dm <- withCallingHandlers(
    suppressWarnings(frm_sample(fit, init = "last.par.best", chains = 1,
                                iter = 20, refresh = 0, seed = 19)),
    message = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_false(any(grepl("chains start at random", said, fixed = TRUE)))
  # no exact gp(), nothing moves: the mode stays the default
  said <- character()
  dd <- gp_low_noise_data(0.3)
  dd$g <- factor(rep(1:5, length.out = nrow(dd)))
  f2 <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ms <- withCallingHandlers(
    suppressWarnings(frm_sample(f2, chains = 1, iter = 20, refresh = 0,
                                seed = 19)),
    message = function(m) {
      said <<- c(said, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_false(any(grepl("chains start at random", said, fixed = TRUE)))
})

