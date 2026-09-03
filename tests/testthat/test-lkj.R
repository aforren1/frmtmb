# The LKJ prior on frmtmb's own correlation parameters.
#
# The claim under test is a DENSITY claim, so the tests are density
# tests, in three layers: (a) the normalizing constant against the
# published LKJ constant, (b) the change of variables against a numeric
# Jacobian of the map t -> C, at several dimensions and shapes, and (c)
# distribution level - sample the prior alone and check every pairwise
# correlation against its closed-form LKJ marginal, which is
# Beta(eta - 1 + d/2, eta - 1 + d/2) rescaled to (-1, 1).

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

# t -> the strictly-lower entries of the correlation matrix, frmtmb's
# own map (us_chol_L is what every covstruct uses)
cor_of_t <- function(t, d) {
  L <- frmtmb:::us_chol_L(t, d)
  C <- L %*% t(L)
  C[lower.tri(C)]
}

# log |det| of a numeric central-difference Jacobian, so the test does
# not depend on a differentiation package
log_abs_det_jac <- function(f, x, ...) {
  h <- 1e-5 * pmax(1, abs(x))
  J <- vapply(seq_along(x), function(i) {
    xp <- x; xp[i] <- xp[i] + h[i]
    xm <- x; xm[i] <- xm[i] - h[i]
    (f(xp, ...) - f(xm, ...)) / (2 * h[i])
  }, numeric(length(x)))
  log(abs(det(as.matrix(J))))
}

# LKJ (2009) closed form for the normalizing constant, independent of
# the per-row assembly frmtmb uses
lkj_logc_published <- function(eta, d) {
  k <- seq_len(d - 1L)
  log(2) * sum((2 * (eta - 1) + d - k) * (d - k)) +
    sum((d - k) * lbeta(eta + (d - 1 - k) / 2, eta + (d - 1 - k) / 2))
}

ld_at <- function(t, eta, d) {
  dist <- frmtmb:::lkj_dist(eta, list(kind = "chol", d = d,
                                      idx = seq_along(t)))
  as.numeric(frmtmb:::lkj_logdens(t, dist))
}

## ---- (a) the normalizing constant ------------------------------------

test_that("the assembled constant is the published LKJ constant", {
  for (d in 2:6) {
    for (eta in c(0.5, 0.9, 1, 1.5, 2, 5)) {
      expect_equal(frmtmb:::lkj_lognorm(eta, d),
                   lkj_logc_published(eta, d), tolerance = 1e-12,
                   info = paste("d =", d, "eta =", eta))
    }
  }
  # a block with no correlation has nothing to normalize
  expect_equal(frmtmb:::lkj_lognorm(1, 1L), 0)
})

## ---- (b) the change of variables -------------------------------------

test_that("the density is LKJ times the exact Jacobian of t -> C", {
  set.seed(6)
  for (d in 2:5) {
    for (eta in c(0.8, 1, 2, 3.5)) {
      for (rep in 1:3) {
        t <- stats::rnorm(d * (d - 1) / 2, 0, 1.2)
        L <- frmtmb:::us_chol_L(t, d)
        C <- L %*% t(L)
        ref <- -lkj_logc_published(eta, d) +
          (eta - 1) * as.numeric(determinant(C)$modulus) +
          log_abs_det_jac(cor_of_t, t, d = d)
        expect_equal(ld_at(t, eta, d), ref, tolerance = 1e-6,
                     info = paste("d =", d, "eta =", eta))
      }
    }
  }
})

test_that("d = 2 is the closed form exactly", {
  tg <- seq(-40, 40, length.out = 401)
  for (eta in c(1, 2, 5)) {
    lhs <- vapply(tg, ld_at, 0, eta = eta, d = 2L)
    # p(t) = (1 + t^2)^-(eta + 1/2) / (2^(2 eta - 1) B(eta, eta))
    rhs <- -(eta + 0.5) * log1p(tg^2) -
      (log(2) * (2 * eta - 1) + lbeta(eta, eta))
    expect_lt(max(abs(exp(lhs - rhs) - 1)), 1e-10)
    # and on the correlation scale it is (1 - rho^2)^(eta - 1)
    rho <- tg / sqrt(1 + tg^2)
    p_rho <- exp(lhs) / (1 + tg^2)^(-1.5)
    ref <- (1 - rho^2)^(eta - 1) /
      exp(log(2) * (2 * eta - 1) + lbeta(eta, eta))
    expect_lt(max(abs(p_rho / ref - 1)), 1e-10)
  }
  # eta = 1 is uniform on rho, so the implied density on t is exactly
  # (1 + t^2)^(-3/2) / 2
  expect_lt(max(abs(exp(vapply(tg, ld_at, 0, eta = 1, d = 2L)) /
                      ((1 + tg^2)^(-1.5) / 2) - 1)), 1e-12)
})

test_that("the one-parameter maps are proper densities", {
  # ar1 / hetar1: rho = t / sqrt(1 + t^2), the d = 2 form exactly
  for (eta in c(0.7, 1, 2, 4)) {
    da <- frmtmb:::lkj_dist(eta, list(kind = "ar1", d = 5L, idx = 1L))
    f <- function(t) exp(vapply(t, function(u) {
      as.numeric(frmtmb:::lkj_logdens(u, da))
    }, 0))
    expect_equal(stats::integrate(f, -Inf, Inf)$value, 1,
                 tolerance = 1e-6, info = paste("ar1 eta =", eta))
    expect_equal(f(0.7), exp(ld_at(0.7, eta, 2L)), tolerance = 1e-12)
  }
  # cs / homcs: the scaled logistic onto (-1/(d - 1), 1), renormalized
  # over that window
  for (d in c(2L, 3L, 6L)) {
    for (eta in c(1, 2)) {
      dc <- frmtmb:::lkj_dist(eta, list(kind = "cs", d = d, idx = 1L))
      f <- function(t) exp(vapply(t, function(u) {
        as.numeric(frmtmb:::lkj_logdens(u, dc))
      }, 0))
      expect_equal(stats::integrate(f, -Inf, Inf)$value, 1,
                   tolerance = 1e-6,
                   info = paste("cs d =", d, "eta =", eta))
    }
  }
  # at d = 2 the cs window is the whole range, so cs and ar1 are the
  # same density on rho, only differently parameterized
  dc <- frmtmb:::lkj_dist(2, list(kind = "cs", d = 2L, idx = 1L))
  t <- c(-2, -0.3, 0.5, 3)
  rho <- -1 + 2 / (1 + exp(-t))
  p_cs <- exp(vapply(t, function(u) {
    as.numeric(frmtmb:::lkj_logdens(u, dc))
  }, 0)) / (2 * exp(-t) / (1 + exp(-t))^2)
  expect_equal(p_cs, (1 - rho^2) / exp(log(2) * 3 + lbeta(2, 2)),
               tolerance = 1e-10)
  # the saturating logistic gives -Inf, never NaN, at eta = 1
  d1 <- frmtmb:::lkj_dist(1, list(kind = "cs", d = 4L, idx = 1L))
  expect_false(is.nan(as.numeric(frmtmb:::lkj_logdens(-800, d1))))
})

## ---- (c) distribution level ------------------------------------------

#' Sample the prior ALONE: no data, no likelihood, just the taped LKJ
#' density over one block's correlation parameters.
prior_rho <- function(eta, d, seed) {
  k <- d * (d - 1L) / 2L
  dist <- frmtmb:::lkj_dist(eta, list(kind = "chol", d = d,
                                      idx = seq_len(k)))
  obj <- RTMB::MakeADFun(
    function(p) -frmtmb:::lkj_logdens(p$t, dist),
    list(t = numeric(k)), silent = TRUE)
  # "random", not the mode: the mode init passes a length-1 vector as a
  # scalar and Stan refuses the k = 1 case
  fit <- suppressWarnings(tmbstan::tmbstan(
    obj, chains = 1, iter = 8000, refresh = 0, seed = seed,
    init = "random"))
  m <- as.matrix(fit)[, seq_len(k), drop = FALSE]
  m <- m[seq(1L, nrow(m), by = 4L), , drop = FALSE]   # thin
  # filled row by row, not apply()ed: with one correlation the
  # transposed apply() would be a 1 x n matrix of draws
  rho <- matrix(0, nrow(m), k)
  for (i in seq_len(nrow(m))) rho[i, ] <- cor_of_t(m[i, ], d)
  rho
}

test_that("every pairwise correlation has the LKJ marginal", {
  skip_if_not(sampler_gates_on(), "chain-agreement gates are off")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  seed <- 0
  for (d in c(2L, 3L, 4L)) {
    for (eta in c(1, 2)) {
      seed <- seed + 1L
      rho <- prior_rho(eta, d, seed)
      a <- eta - 1 + d / 2      # the closed-form marginal shape
      for (j in seq_len(ncol(rho))) {
        ks <- suppressWarnings(
          stats::ks.test(rho[, j], function(q) stats::pbeta((q + 1) / 2,
                                                            a, a)))
        # a chain of ~2000 thinned draws: the iid 5% critical value is
        # 0.030, and the gate is set at twice it so that autocorrelation
        # cannot fail a correct density, while a wrong exponent (the
        # failure mode this guards) moves D by 0.2 or more
        expect_lt(as.numeric(ks$statistic), 0.06,
                  label = paste0("KS D (d = ", d, ", eta = ", eta,
                                 ", pair ", j, ")"))
        # the marginal's first two moments, against the chain's own
        # Monte Carlo spread rather than a fixed number
        mcse <- stats::sd(rho[, j]) / sqrt(length(rho[, j]) / 10)
        expect_lt(abs(mean(rho[, j])), 5 * mcse)
        expect_lt(abs(stats::sd(rho[, j]) - sqrt(1 / (2 * a + 1))),
                  0.1 * sqrt(1 / (2 * a + 1)))
      }
    }
  }
})

## ---- the grammar ------------------------------------------------------

test_that("lkj() parses and belongs to class cor, both ways", {
  expect_equal(frmtmb:::parse_prior_dist("lkj(2)")$eta, 2)
  expect_equal(prior_lkj(2)$kind, "lkj")
  expect_error(prior_lkj(0), "finite positive")
  expect_error(prior_lkj(c(1, 2)), "finite positive")
  expect_error(set_prior("lkj(2)", class = "sd"), "belongs to class")
  expect_error(set_prior("normal(0, 1)", class = "cor"),
               "takes an lkj")
  expect_error(set_prior("", class = "cor", lb = 0), "takes an lkj")
  # a bound belongs to one parameter, so it is refused rather than
  # silently dropped on a class that names a whole matrix
  expect_error(set_prior("lkj(2)", class = "cor", ub = 0.5),
               "takes no lb/ub")
  pl <- set_prior("lkj(2)", class = "cor")
  expect_match(utils::capture.output(print(pl)), "lkj(2) class=cor",
               fixed = TRUE)
  # the objects are interchangeable with the string spelling
  expect_equal(unclass(set_prior(prior_lkj(2), class = "cor")),
               unclass(pl))
})

lkj_data <- function(seed = 77, n = 160L, ng = 16L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  u <- matrix(stats::rnorm(2 * ng), 2)
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + 0.8 * u[1, dd$g] +
                         0.4 * u[2, dd$g] * dd$x, 1)
  dd
}

test_that("class cor addresses a block, and names what it cannot fit", {
  dd <- lkj_data()
  fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
  bk <- fit$frame$re_blocks[[1L]]
  ri <- frmtmb:::resolve_prior_input(fit, set_prior("lkj(2)",
                                                   class = "cor"))
  expect_length(ri$entries, 1L)
  e <- ri$entries[[1L]]
  # the whole correlation segment of THAT block, and nothing else
  expect_equal(e$idx, bk$theta_idx[frmtmb:::block_cor_spec(bk)$idx])
  expect_equal(e$dist$kind, "lkj")
  expect_equal(e$dist$map, "chol")

  # group = narrows exactly as it does for class "sd"
  expect_length(frmtmb:::resolve_prior_input(
    fit, set_prior("lkj(2)", class = "cor", group = "g"))$entries, 1L)
  expect_error(frmtmb:::resolve_prior_input(
    fit, set_prior("lkj(2)", class = "cor", group = "zzz")),
    "No random-effect correlations")

  # an uncorrelated model has nothing to address, and the refusal says
  # what the model does have
  f0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(f0, set_prior("lkj(1)",
                                                          class = "cor")),
               "no correlation parameter")
  # ... including a model with no random effects at all
  f1 <- frm(bf(y ~ x) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(f1, set_prior("lkj(1)",
                                                          class = "cor")),
               "no random-effect blocks")
})

test_that("toep is refused by name, and stays out of the defaults", {
  set.seed(78)
  dd <- data.frame(t = factor(rep(1:4, 40)),
                   g = factor(rep(1:20, each = 8)))
  dd$y <- stats::rnorm(160, rep(stats::rnorm(20, 0, 0.7), each = 8), 1)
  fit <- suppressWarnings(
    frm(bf(y ~ 1 + toep(t | g)) + gaussian(), data = dd))
  expect_identical(fit$frame$re_blocks[[1L]]$covstruct, "toep")
  expect_error(frmtmb:::resolve_prior_input(fit, set_prior("lkj(1)",
                                                           class = "cor")),
               "positive definite")
  # the default builder skips it rather than failing, and says so
  uf <- frm(bf(y ~ 1 + toep(t | g)) + gaussian(), dd,
            dry_run = "objective")
  cls <- vapply(unclass(frmtmb:::default_priors_for(uf)), `[[`, "",
                "class")
  expect_false("cor" %in% cls)
  expect_match(frmtmb:::default_prior_notes(uf), "no LKJ density fits")
  # and it cannot be non-centered either, for the same reason
  expect_false(frmtmb:::ncp_eligible(fit$frame$re_blocks[[1L]]))
})

test_that("cs, ar1 and gr(cov =) take the prior their map calls for", {
  set.seed(79)
  dd <- data.frame(t = factor(rep(1:4, 40)),
                   g = factor(rep(1:20, each = 8)))
  dd$y <- stats::rnorm(160, rep(stats::rnorm(20, 0, 0.7), each = 8), 1)
  maps <- c(cs = "cs", homcs = "cs", ar1 = "ar1", hetar1 = "ar1")
  for (nm in names(maps)) {
    # ar1() wants the factor without an intercept; cs() takes it too
    ff <- stats::as.formula(paste0("y ~ 1 + ", nm, "(t + 0 | g)"))
    fit <- suppressWarnings(frm(bf(ff) + gaussian(), data = dd))
    expect_identical(fit$frame$re_blocks[[1L]]$covstruct, nm)
    e <- frmtmb:::resolve_prior_input(
      fit, set_prior("lkj(2)", class = "cor"))$entries[[1L]]
    expect_length(e$idx, 1L)
    expect_equal(e$dist$map, unname(maps[nm]), info = nm)
    # and the position really is the structure's correlation: the block
    # density changes when it moves, the standard deviations do not
    expect_false(e$idx %in%
                   fit$frame$re_blocks[[1L]]$theta_idx[
                     frmtmb:::covstruct_registry[[nm]]$sd_idx(4L)])
  }
})

test_that("later wins between class cor and the class theta hatch", {
  dd <- lkj_data()
  fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
  # the LKJ term is retired when a raw prior lands on a position it
  # covers, so the two densities are never added together
  ri <- frmtmb:::resolve_prior_input(
    fit, set_prior("lkj(2)", class = "cor") +
      set_prior("normal(0, 1)", class = "theta", coef = "theta_3"))
  kinds <- vapply(ri$entries, function(e) e$dist$kind, "")
  expect_false("lkj" %in% kinds)
  # and the other way round
  ri2 <- frmtmb:::resolve_prior_input(
    fit, set_prior("normal(0, 1)", class = "theta", coef = "theta_3") +
      set_prior("lkj(2)", class = "cor"))
  expect_equal(vapply(ri2$entries, function(e) e$dist$kind, ""), "lkj")
  # a class "sd" prior sits beside it: different positions, both kept
  ri3 <- frmtmb:::resolve_prior_input(
    fit, set_prior("lkj(2)", class = "cor") +
      set_prior("exponential(1)", class = "sd"))
  expect_setequal(vapply(ri3$entries, function(e) e$dist$kind, ""),
                  c("lkj", "exponential"))
})

test_that("the named-list spelling refuses an LKJ prior by name", {
  dd <- lkj_data()
  fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)
  expect_error(frmtmb:::resolve_prior_input(
    fit, list(theta_3 = prior_lkj(2))), "class = \"cor\"", fixed = TRUE)
})

## ---- the formula route ------------------------------------------------

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

## ---- what the prior does to a fit ------------------------------------

test_that("an LKJ MAP penalty pulls the correlation toward zero", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  ff <- bf(Reaction ~ Days + (Days | Subject)) + gaussian()
  f0 <- frm(ff, data = sleepstudy)
  cor_of <- function(f) {
    V <- VarCorr(f)[[1L]]
    V[1L, 2L] / sqrt(V[1L, 1L] * V[2L, 2L])
  }
  # eta > 1 concentrates toward the identity, so the MAP correlation
  # moves toward 0 and keeps moving as eta grows
  f2 <- frm(ff, data = sleepstudy,
            prior = set_prior("lkj(4)", class = "cor"))
  f3 <- frm(ff, data = sleepstudy,
            prior = set_prior("lkj(50)", class = "cor"))
  expect_lt(abs(cor_of(f2)), abs(cor_of(f0)) + 1e-8)
  expect_lt(abs(cor_of(f3)), abs(cor_of(f2)))
  expect_lt(abs(cor_of(f3)), 0.02)
  # a MAP fit records what it was penalized with, and the penalty is
  # what moved the estimate: the pure likelihood at the MAP estimate is
  # WORSE than at the ML one, which is the definition of a penalty
  expect_s3_class(prior_summary(f3), "frmtmb_priorlist")
  # f0's objective is the LIKELIHOOD alone (no prior term), and it is
  # worse at the MAP estimate than at the ML one
  expect_gt(as.numeric(f0$obj$fn(f3$opt$par)),
            as.numeric(f0$obj$fn(f0$opt$par)))
  # eta = 1 is uniform over correlation matrices, but NOT flat on the
  # parameter, so it still moves the mode; it is a proper prior
  f1 <- frm(ff, data = sleepstudy,
            prior = set_prior("lkj(1)", class = "cor"))
  expect_lt(abs(cor_of(f1)), abs(cor_of(f0)) + 1e-8)
})

test_that("frm() is untouched when no prior is given", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
             data = sleepstudy)
  # the recorded pre-LKJ value of this fit (v0.38.0, same machine and
  # data): the correlation default is a SAMPLING default and never
  # enters frm(), so the ML fit reproduces to the last digit
  expect_equal(as.numeric(logLik(fit)), -875.9697, tolerance = 1e-4)
  expect_null(fit$prior)
  # and the objective the fit holds carries no prior term
  expect_equal(as.numeric(fit$obj$fn(fit$opt$par)),
               -as.numeric(logLik(fit)), tolerance = 1e-8)
})
