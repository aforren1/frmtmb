# Cross-validation against brms itself. frmtmb reimplements the brms
# formula grammar, so brms is the reference for what a formula MEANS.
#
# Tier 1 (this file's bulk) compares design objects: brms::make_standata()
# builds every design matrix, index vector, and addition-term vector
# without touching Stan, so the comparison is fast and safe for
# NOT_CRAN runs. Tier 2 fits Stan models and is opt-in (see
# skip_unless_brms_fit()).
#
# Naming conventions differ and are never asserted: brms calls the
# intercept column "Intercept" (Stan identifiers cannot hold
# parentheses), suffixes dpar designs as X_<dpar>, and numbers grouping
# factors instead of naming them. Values are what must agree.

# ---------------------------------------------------------------------
# Tier 1: design agreement through brms's data-generating functions
# ---------------------------------------------------------------------

brms_agree_data <- function(seed = 20240501, n = 200) {
  set.seed(seed)
  d <- data.frame(
    x = rnorm(n),
    z = rnorm(n),
    inc = sample(0:3, n, TRUE),
    fac = factor(sample(c("a", "b", "c"), n, TRUE)),
    g = factor(rep(seq_len(n / 10), each = 10)),
    nt = rep(12, n),
    w = runif(n, 0.5, 2),
    sdy = runif(n, 0.1, 0.5)
  )
  d$yg <- 1 + 0.7 * d$x - 0.3 * d$z + rnorm(n)
  d$yc <- rpois(n, 3)
  d$yb <- rbinom(n, 12, 0.4)
  d$ord <- sample(1:4, n, TRUE)
  d$cl <- rep(c(TRUE, FALSE), length.out = n)
  d
}

test_that("population-level designs match brms standata X", {
  skip_unless_brms()
  dd <- brms_agree_data()

  cases <- list(
    list(frm = bf(yg ~ x + z) + gaussian(),
         brm = brms::bf(yg ~ x + z), fam = gaussian(), key = "yg.mu"),
    # factor contrasts: brms and frmtmb both go through model.matrix
    list(frm = bf(yg ~ fac * x) + gaussian(),
         brm = brms::bf(yg ~ fac * x), fam = gaussian(), key = "yg.mu"),
    # offsets must land in the same place and on the same scale
    list(frm = bf(yc ~ x + offset(log(w))) + poisson(),
         brm = brms::bf(yc ~ x + offset(log(w))), fam = poisson(),
         key = "yc.mu"),
    list(frm = bf(yg ~ 0 + fac) + gaussian(),
         brm = brms::bf(yg ~ 0 + fac), fam = gaussian(), key = "yg.mu")
  )
  for (cs in cases) {
    sd <- brms_standata(cs$brm, data = dd, family = cs$fam)
    fr <- frm(cs$frm, data = dd, dry_run = "frame")
    expect_design_equal(fr$linpreds[[cs$key]]$X, sd$X)
  }

  # offsets: brms exposes them as `offsets`, we store them on the linpred
  sd <- brms_standata(brms::bf(yc ~ x + offset(log(w))), data = dd,
                      family = poisson())
  fr <- frm(bf(yc ~ x + offset(log(w))) + poisson(), data = dd,
            dry_run = "frame")
  expect_vector_equal(fr$linpreds[["yc.mu"]]$offset, sd$offsets, tol = 1e-12)
})

test_that("random-effect designs match brms Z_*/J_* and level order", {
  skip_unless_brms()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  sd <- brms_standata(brms::bf(Reaction ~ Days + (Days | Subject)),
                      data = sleepstudy, family = gaussian())
  fr <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
            data = sleepstudy, dry_run = "frame")
  bk <- fr$re_blocks[[1]]
  expect_identical(bk$n_levels, as.integer(sd$N_1))
  expect_identical(bk$dim, as.integer(sd$M_1))
  expect_identical(length(bk$theta_idx), 2L + as.integer(sd$NC_1))
  expect_identical(bk$levels, levels(sleepstudy$Subject))

  # brms stores one column per RE coefficient plus a level index; our Z
  # is the sparse expansion, level-major within a block
  Z <- as.matrix(fr$linpreds[["Reaction.mu"]]$Z)
  Zb <- matrix(0, nrow(Z), ncol(Z))
  idx <- cbind(seq_len(nrow(Z)), (sd$J_1 - 1L) * bk$dim + 1L)
  Zb[idx] <- sd$Z_1_1
  Zb[cbind(idx[, 1], idx[, 2] + 1L)] <- sd$Z_1_2
  expect_design_equal(Z, Zb)
})

test_that("|ID| merging matches brms's cross-formula grouping block", {
  skip_unless_brms()
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")

  bform <- brms::bf(Reaction ~ Days + (1 | p | Subject)) +
    brms::bf(Days ~ 1 + (1 | p | Subject)) + brms::set_rescor(FALSE)
  sd <- brms_standata(bform, data = sleepstudy, family = gaussian())
  fr <- frm(bf(Reaction ~ Days + (1 | p | Subject)) +
              bf(Days ~ 1 + (1 | p | Subject)) + gaussian(),
            data = sleepstudy, dry_run = "frame")

  # both packages merge the two intercepts into ONE correlated block
  expect_length(fr$re_blocks, 1)
  expect_identical(fr$re_blocks[[1]]$dim, as.integer(sd$M_1))
  expect_identical(fr$re_blocks[[1]]$n_levels, as.integer(sd$N_1))
  expect_identical(as.integer(sd$NC_1), 1L)
})

test_that("distributional dpar designs match brms X_<dpar>", {
  skip_unless_brms()
  dd <- brms_agree_data()

  sd <- brms_standata(brms::bf(yg ~ x, sigma ~ z + fac), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg ~ x, sigma ~ z + fac) + gaussian(), data = dd,
            dry_run = "frame")
  expect_design_equal(fr$linpreds[["yg.mu"]]$X, sd$X)
  expect_design_equal(fr$linpreds[["yg.sigma"]]$X, sd$X_sigma)

  sd <- brms_standata(brms::bf(yc ~ x, zi ~ z), data = dd,
                      family = brms::zero_inflated_poisson())
  fr <- frm(bf(yc ~ x, zi ~ z) + zero_inflated_poisson(), data = dd,
            dry_run = "frame")
  expect_design_equal(fr$linpreds[["yc.zi"]]$X, sd$X_zi)

  sd <- brms_standata(brms::bf(yc ~ x, hu ~ z), data = dd,
                      family = brms::hurdle_poisson())
  fr <- frm(bf(yc ~ x, hu ~ z) + hurdle_poisson(), data = dd,
            dry_run = "frame")
  expect_design_equal(fr$linpreds[["yc.hu"]]$X, sd$X_hu)

  # the dpar sets themselves agree with what brms considers valid
  bt <- brms::brmsterms(brms::bf(yc ~ x, zi ~ z,
                                 family = brms::zero_inflated_poisson()))
  expect_identical(names(bt$dpars), c("mu", "zi"))
})

test_that("ordinal families agree on thresholds and the dropped intercept", {
  skip_unless_brms()
  dd <- brms_agree_data()

  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    sd <- brms_standata(brms::bf(ord ~ x + z), data = dd,
                        family = get(fam, asNamespace("brms"))())
    fr <- frm(bf(ord ~ x + z) + get(fam)(), data = dd, dry_run = "frame")
    # thresholds replace the intercept in both packages
    expect_length(fr$par_template$tau_raw, as.integer(sd$nthres))
    expect_design_equal(fr$linpreds[["ord.mu"]]$X, sd$X)
    expect_false("(Intercept)" %in% colnames(fr$linpreds[["ord.mu"]]$X))
  }
})

test_that("cs() category-specific effects match brms Xcs", {
  skip_unless_brms()
  skip_if_not_installed("brms")
  data(inhaler, package = "brms")

  sd <- brms_standata(brms::bf(rating ~ period + carry + cs(treat)),
                      data = inhaler, family = brms::sratio())
  fr <- frm(bf(rating ~ period + carry + cs(treat)) + sratio(),
            data = inhaler, dry_run = "frame")
  cs <- fr$linpreds[["rating.mu"]]$cs
  expect_length(cs, as.integer(sd$Kcs))
  expect_vector_equal(cs[[1]]$vals, sd$Xcs[, 1], tol = 1e-12)
  # brms keeps cs() out of X, as we do; both carry nthres coefficients
  expect_design_equal(fr$linpreds[["rating.mu"]]$X, sd$X)
  expect_length(fr$par_template[[cs[[1]]$par]], as.integer(sd$nthres))
})

test_that("mo() codes match brms Xmo and the simplex dimension matches Jmo", {
  skip_unless_brms()
  dd <- brms_agree_data()

  sd <- brms_standata(brms::bf(yg ~ mo(inc) + z), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg ~ mo(inc) + z) + gaussian(), data = dd,
            dry_run = "frame")
  lp <- fr$linpreds[["yg.mu"]]
  mo <- lp$mo[[1]]
  expect_vector_equal(mo$codes, sd$Xmo_1, tol = 1e-12)
  expect_identical(mo$D, as.integer(sd$Jmo[1]))

  # DIVERGENCE (convention): brms samples the simplex itself, so it has
  # Jmo free-ish entries under a Dirichlet; we hold the D - 1 free
  # softmax coordinates and rebuild the simplex, which is the same
  # manifold with one fewer stored number.
  expect_length(fr$par_template[[mo$zeta]], as.integer(sd$Jmo[1]) - 1L)

  # DIVERGENCE (convention): brms keeps mo() out of X entirely (Ksp
  # coefficients live in bsp); we carry a zero placeholder column so the
  # coefficient bookkeeping stays uniform. Drop it and X agrees.
  expect_identical(ncol(lp$X), ncol(sd$X) + as.integer(sd$Ksp))
  expect_design_equal(lp$X[, -mo$col, drop = FALSE], sd$X)

  # ordered factors: brms codes them 0..K-1 the same way
  dd$oinc <- factor(dd$inc, ordered = TRUE)
  sd <- brms_standata(brms::bf(yg ~ mo(oinc)), data = dd, family = gaussian())
  fr <- frm(bf(yg ~ mo(oinc)) + gaussian(), data = dd, dry_run = "frame")
  expect_vector_equal(fr$linpreds[["yg.mu"]]$mo[[1]]$codes, sd$Xmo_1, tol = 1e-12)
})

test_that("exact gp() matches brms's unique-position grouping", {
  skip_unless_brms()
  dd <- brms_agree_data(n = 60)
  dd$xd <- round(dd$x, 1)   # deliberate ties, to exercise the grouping

  sd <- brms_standata(brms::bf(yg ~ gp(xd)), data = dd, family = gaussian())
  fr <- frm(bf(yg ~ gp(xd)) + gaussian(), data = dd, dry_run = "frame")
  bk <- fr$re_blocks[[1]]
  expect_identical(bk$covstruct, "gp")
  expect_identical(bk$gp_D, as.integer(sd$Dgp_1))
  # one latent value per distinct coordinate, in both packages
  expect_identical(bk$dim, as.integer(sd$Nsubgp_1))

  # our Z is the 0/1 indicator of the position index; brms ships the
  # index itself. Level numbering is an internal detail, so compare the
  # partitions the two indices induce, row pair by row pair.
  Z <- as.matrix(fr$linpreds[["yg.mu"]]$Z)
  expect_true(all(rowSums(Z) == 1))
  ours <- apply(Z, 1, which.max)
  expect_true(all(outer(ours, ours, "==") == outer(sd$Jgp_1, sd$Jgp_1, "==")))
})

test_that("gp(k =) HSGP basis matches brms exactly", {
  skip_unless_brms()
  set.seed(2)
  ds <- data.frame(x = runif(200, 0, 10))
  ds$y <- sin(ds$x) + rnorm(200, 0, 0.4)

  sd <- brms_standata(brms::bf(y ~ gp(x, k = 25)), data = ds,
                      family = gaussian())
  fr <- frm(bf(y ~ gp(x, k = 25)) + gaussian(), data = ds,
            dry_run = "frame")
  # basis size is k^D in both packages
  expect_identical(fr$re_blocks[[1]]$dim, as.integer(sd$NBgp_1))
  expect_identical(fr$re_blocks[[1]]$covstruct, "hsgp")

  # same input convention: brms rescales the covariates by the largest
  # pairwise distance and takes L = c * max(1, range(scaled)), which
  # after that scaling is always exactly c. We do the same, so the
  # boundary, the basis and the spectral frequencies agree directly.
  gi <- fr$linpreds[["y.mu"]]$gps[[1]]
  expect_vector_equal(gi$L, as.numeric(sd$Lgp_1), tol = 1e-12)
  expect_equal(as.numeric(sd$Lgp_1), 1.25)
  # brms builds its basis over distinct positions and indexes it with
  # Jgp; our Z carries one basis row per observation
  expect_design_equal(fr$linpreds[["y.mu"]]$Z,
                      sd$Xgp_1[sd$Jgp_1, , drop = FALSE], tol = 1e-12)
  expect_vector_equal(as.vector(gi$omega), as.vector(sd$slambda_1),
                      tol = 1e-12)

  # multi-dimensional gp: one shared scale over all coordinates, one
  # boundary per dimension, and the same k^D basis
  dd <- brms_agree_data(n = 120)
  sd2 <- brms_standata(brms::bf(yg ~ gp(x, z, k = 5)), data = dd,
                       family = gaussian())
  fr2 <- frm(bf(yg ~ gp(x, z, k = 5)) + gaussian(), data = dd,
             dry_run = "frame")
  expect_identical(fr2$re_blocks[[1]]$dim, as.integer(sd2$NBgp_1))
  expect_identical(fr2$re_blocks[[1]]$gp_D, as.integer(sd2$Dgp_1))
  gi2 <- fr2$linpreds[["yg.mu"]]$gps[[1]]
  expect_vector_equal(gi2$L, as.numeric(sd2$Lgp_1), tol = 1e-12)
  expect_design_equal(fr2$linpreds[["yg.mu"]]$Z,
                      sd2$Xgp_1[sd2$Jgp_1, , drop = FALSE], tol = 1e-12)
  expect_vector_equal(as.vector(gi2$omega), as.vector(sd2$slambda_1),
                      tol = 1e-12)

  # c = is per covariate in brms, and so here
  sd3 <- brms_standata(brms::bf(yg ~ gp(x, z, k = 5, c = c(1.5, 2))),
                       data = dd, family = gaussian())
  fr3 <- frm(bf(yg ~ gp(x, z, k = 5, c = c(1.5, 2))) + gaussian(),
             data = dd, dry_run = "frame")
  expect_vector_equal(fr3$linpreds[["yg.mu"]]$gps[[1]]$L,
                      as.numeric(sd3$Lgp_1), tol = 1e-12)
  expect_design_equal(fr3$linpreds[["yg.mu"]]$Z,
                      sd3$Xgp_1[sd3$Jgp_1, , drop = FALSE], tol = 1e-12)

  # tied coordinates: brms's gr = TRUE collapses duplicate positions
  # before it computes the scale and the center, and so do we
  dt <- data.frame(x = round(runif(150, 0, 10), 1))
  dt$y <- sin(dt$x) + rnorm(150, 0, 0.3)
  sd4 <- brms_standata(brms::bf(y ~ gp(x, k = 12)), data = dt,
                       family = gaussian())
  fr4 <- frm(bf(y ~ gp(x, k = 12)) + gaussian(), data = dt,
             dry_run = "frame")
  expect_design_equal(fr4$linpreds[["y.mu"]]$Z,
                      sd4$Xgp_1[sd4$Jgp_1, , drop = FALSE], tol = 1e-12)
})

test_that("mgcv smooths agree with brms's Zs/Xs bases", {
  skip_unless_brms()
  dd <- brms_agree_data()

  # t2(): identical bases, element for element
  sd <- brms_standata(brms::bf(yg ~ t2(x, z)), data = dd, family = gaussian())
  fr <- frm(bf(yg ~ t2(x, z)) + gaussian(), data = dd, dry_run = "frame")
  sm <- fr$linpreds[["yg.mu"]]$smooths[[1]]
  expect_identical(as.integer(sm$nr), as.integer(sd$knots_1))
  expect_identical(length(sm$nr), as.integer(sd$nb_1))
  expect_identical(sm$nf, as.integer(sd$Ks))
  expect_design_equal(fr$linpreds[["yg.mu"]]$Z,
                      cbind(sd$Zs_1_1, sd$Zs_1_2, sd$Zs_1_3), tol = 1e-10)

  # s(): same dimensions, same column space, different basis rotation.
  # DIVERGENCE (convention): brms calls mgcv::smoothCon() with
  # diagonal.penalty = TRUE, we do not. The models are identical - the
  # spans below coincide - but the individual basis columns are a
  # reparameterization of each other, so coefficient values are not
  # comparable one to one.
  sd <- brms_standata(brms::bf(yg ~ s(x)), data = dd, family = gaussian())
  fr <- frm(bf(yg ~ s(x)) + gaussian(), data = dd, dry_run = "frame")
  sm <- fr$linpreds[["yg.mu"]]$smooths[[1]]
  expect_identical(as.integer(sm$nr), as.integer(sd$knots_1))
  expect_identical(sm$nf, as.integer(sd$Ks))
  expect_span_equal(fr$linpreds[["yg.mu"]]$Z, sd$Zs_1_1)
  expect_span_equal(cbind(fr$linpreds[["yg.mu"]]$X), cbind(sd$X, sd$Xs))
  # and the reparameterization is exactly diagonal.penalty
  scl <- mgcv::smoothCon(mgcv::s(x), data = dd, absorb.cons = TRUE,
                         diagonal.penalty = TRUE)
  re2 <- mgcv::smooth2random(scl[[1]], names(dd), type = 2)
  expect_design_equal(re2$rand[[1]], sd$Zs_1_1)

  # smooths in a dpar keep the same suffixing scheme
  sd <- brms_standata(brms::bf(yg ~ x, sigma ~ s(z)), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg ~ x, sigma ~ s(z)) + gaussian(), data = dd,
            dry_run = "frame")
  expect_identical(
    as.integer(fr$linpreds[["yg.sigma"]]$smooths[[1]]$nr),
    as.integer(sd$knots_sigma_1)
  )
})

test_that("addition terms match brms's standata vectors", {
  skip_unless_brms()
  dd <- brms_agree_data()

  sd <- brms_standata(brms::bf(yb | trials(nt) ~ x), data = dd,
                      family = binomial())
  fr <- frm(bf(yb | trials(nt) ~ x) + binomial(), data = dd,
            dry_run = "frame")
  expect_vector_equal(fr$aterm_values$yb$trials, sd$trials, tol = 1e-12)
  expect_vector_equal(fr$y$yb, sd$Y, tol = 1e-12)

  sd <- brms_standata(brms::bf(yg | weights(w) ~ x), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg | weights(w) ~ x) + gaussian(), data = dd,
            dry_run = "frame")
  expect_vector_equal(fr$aterm_values$yg$weights, sd$weights, tol = 1e-12)

  # censoring: brms codes -1/0/1/2, and maps logicals the same way
  sd <- brms_standata(brms::bf(yg | cens(cl) ~ x), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg | cens(cl) ~ x) + gaussian(), data = dd,
            dry_run = "frame")
  expect_vector_equal(fr$aterm_values$yg$cens, sd$cens, tol = 1e-12)

  dd$lb <- -4; dd$ub <- 6
  sd <- brms_standata(brms::bf(yg | trunc(lb = lb, ub = ub) ~ x), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg | trunc(lb = lb, ub = ub) ~ x) + gaussian(), data = dd,
            dry_run = "frame")
  expect_vector_equal(fr$aterm_values$yg$trunc_lb, sd$lb, tol = 1e-12)
  expect_vector_equal(fr$aterm_values$yg$trunc_ub, sd$ub, tol = 1e-12)

  sd <- brms_standata(brms::bf(yg | se(sdy) ~ x), data = dd,
                      family = gaussian())
  fr <- frm(bf(yg | se(sdy) ~ x) + gaussian(), data = dd,
            dry_run = "frame")
  expect_vector_equal(fr$aterm_values$yg$se, sd$se, tol = 1e-12)
  # se() alone: the residual SD is the known se, and brms says so by
  # shipping the constant sigma = 0
  expect_null(fr$aterm_values$yg$se_sigma)
  expect_identical(as.numeric(sd$sigma), 0)
  # se(x, sigma = TRUE): sigma becomes a real parameter, so brms drops
  # the constant from the data entirely
  fr <- frm(bf(yg | se(sdy, sigma = TRUE) ~ x) + gaussian(), data = dd,
            dry_run = "frame")
  sd <- brms_standata(brms::bf(yg | se(sdy, sigma = TRUE) ~ x), data = dd,
                      family = gaussian())
  expect_true(fr$aterm_values$yg$se_sigma)
  expect_null(sd$sigma)
})

test_that("censored survival data agrees on the whole design (kidney)", {
  skip_unless_brms()
  data(kidney, package = "brms")

  sd <- brms_standata(brms::bf(time | cens(censored) ~ age + sex + disease),
                      data = kidney, family = brms::lognormal())
  fr <- frm(bf(time | cens(censored) ~ age + sex + disease) + lognormal(),
            data = kidney, dry_run = "frame")
  expect_design_equal(fr$linpreds[["time.mu"]]$X, sd$X)
  expect_vector_equal(fr$aterm_values$time$cens, sd$cens, tol = 1e-12)
  expect_vector_equal(fr$y$time, sd$Y, tol = 1e-12)
})

test_that("multivariate models match brms per-response designs", {
  skip_unless_brms()
  dd <- brms_agree_data()

  bform <- brms::bf(brms::mvbind(yg, x) ~ z + fac) + brms::set_rescor(TRUE)
  sd <- brms_standata(bform, data = dd, family = gaussian())
  fr <- frm(bf(mvbind(yg, x) ~ z + fac) + gaussian() + set_rescor(TRUE),
            data = dd, dry_run = "frame")
  expect_design_equal(fr$linpreds[["yg.mu"]]$X, sd$X_yg)
  expect_design_equal(fr$linpreds[["x.mu"]]$X, sd$X_x)
  # one residual correlation for two responses in both packages
  expect_length(fr$par_template$thetar, as.integer(sd$nrescor))
  expect_identical(length(fr$spec$responses), as.integer(sd$nresp))
})

test_that("nonlinear formulas match brms nlpar designs and covariates", {
  skip_unless_brms()
  dd <- brms_agree_data()

  bform <- brms::bf(yg ~ a * exp(-b * x), a ~ 1, b ~ z, nl = TRUE)
  sd <- brms_standata(bform, data = dd, family = gaussian())
  fr <- frm(bf(yg ~ a * exp(-b * x), a ~ 1, b ~ z, nl = TRUE) + gaussian(),
            data = dd, dry_run = "frame")
  expect_design_equal(fr$linpreds[["yg.a"]]$X, sd$X_a)
  expect_design_equal(fr$linpreds[["yg.b"]]$X, sd$X_b)
  # the covariate of the nonlinear body is plain data in both packages
  # (brms ships it as C_1)
  expect_vector_equal(fr$linpreds[["yg.mu"]]$data_list$x, sd$C_1, tol = 1e-12)
  bt <- brms::brmsterms(bform)
  expect_identical(names(bt$nlpars), fr$spec$responses$yg$nlpars)
})

test_that("mixture dpar naming follows the brms convention", {
  skip_unless_brms()
  dd <- brms_agree_data()

  bfam <- suppressMessages(brms::mixture(gaussian(), gaussian()))
  sd <- brms_standata(brms::bf(yg ~ x), data = dd, family = bfam)
  fr <- frm(bf(yg ~ x) + mixture(gaussian(), gaussian()), data = dd,
            dry_run = "frame")
  expect_design_equal(fr$linpreds[["yg.mu1"]]$X, sd$X_mu1)
  expect_design_equal(fr$linpreds[["yg.mu2"]]$X, sd$X_mu2)
  expect_true(all(c("yg.sigma1", "yg.sigma2") %in% names(fr$linpreds)))

  # DIVERGENCE (convention): brms carries theta1..thetaK under a
  # Dirichlet with a sum-to-one constraint (con_theta has K entries);
  # we carry the K - 1 free logits of the same simplex, so theta2 has
  # no separate parameter here.
  expect_identical(length(sd$con_theta), 2L)
  expect_true("yg.theta1" %in% names(fr$linpreds))
  expect_false("yg.theta2" %in% names(fr$linpreds))
})

test_that("brms get_prior agrees on which special terms exist", {
  skip_unless_brms()
  dd <- brms_agree_data()

  # get_prior's classes are a compact statement of the term types brms
  # found; check that our parse found the same specials.
  classes <- function(bform, family = gaussian()) {
    unique(as.data.frame(brms::get_prior(bform, data = dd,
                                         family = family))$class)
  }
  expect_true("simo" %in% classes(brms::bf(yg ~ mo(inc) + z)))
  expect_length(frm(bf(yg ~ mo(inc) + z) + gaussian(), data = dd,
                    dry_run = "frame")$linpreds[["yg.mu"]]$mo, 1)

  expect_true(all(c("sdgp", "lscale") %in% classes(brms::bf(yg ~ gp(x, k = 8)))))
  expect_length(frm(bf(yg ~ gp(x, k = 8)) + gaussian(), data = dd,
                    dry_run = "frame")$linpreds[["yg.mu"]]$gps, 1)

  expect_true("sds" %in% classes(brms::bf(yg ~ s(x))))
  expect_length(frm(bf(yg ~ s(x)) + gaussian(), data = dd,
                    dry_run = "frame")$linpreds[["yg.mu"]]$smooths, 1)

  cl <- classes(brms::bf(yg ~ x + (x | g), sigma ~ z))
  expect_true(all(c("sd", "cor") %in% cl))
  fr <- frm(bf(yg ~ x + (x | g), sigma ~ z) + gaussian(), data = dd,
            dry_run = "frame")
  expect_length(fr$re_blocks[[1]]$theta_idx, 3)   # 2 sds + 1 correlation
})

# ---------------------------------------------------------------------
# Tier 2: numeric agreement. Stan compiles here (minutes per model), so
# these are opt-in:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true")
#
# brms posterior means are not maximum-likelihood estimates in general.
# Where a model has no latent effects and every prior is flat, the
# posterior mode IS the MLE, so those models are compared through
# rstan::optimizing() on brms's own generated Stan program and agree to
# optimizer precision.
#
# The comparison against posterior means that used to live here is
# gone: the two estimands differ, so the only tolerance that let it
# pass was also wide enough to hide a real divergence. Densities are
# compared at a point in test-brms-likelihood.R instead.
# ---------------------------------------------------------------------

test_that("distributional gaussian ML matches the brms posterior mode", {
  skip_unless_brms_fit()

  set.seed(11)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- 1 + 0.8 * dd$x - 0.4 * dd$z +
    rnorm(n, 0, exp(0.2 + 0.3 * dd$x))

  # empty priors are brms's spelling for improper flat, which makes the
  # posterior mode the MLE
  flat <- c(brms::prior("", class = "Intercept"),
            brms::prior("", class = "b"),
            brms::prior("", class = "Intercept", dpar = "sigma"),
            brms::prior("", class = "b", dpar = "sigma"))
  bform <- brms::bf(y ~ x + z, sigma ~ x)
  code <- brms::make_stancode(bform, data = dd, family = gaussian(),
                              prior = flat)
  sdat <- brms::make_standata(bform, data = dd, family = gaussian(),
                              prior = flat)
  mod <- rstan::stan_model(model_code = code)
  op <- rstan::optimizing(mod, data = sdat, hessian = FALSE,
                          as_vector = TRUE, seed = 1)

  fit <- frm(bf(y ~ x + z, sigma ~ x) + gaussian(), data = dd)
  fe <- fixef(fit)
  expect_vector_equal(op$par[["b_Intercept"]], fe$mu[["(Intercept)"]],
                      tol = 1e-4)
  expect_vector_equal(op$par[["b[1]"]], fe$mu[["x"]], tol = 1e-4)
  expect_vector_equal(op$par[["b[2]"]], fe$mu[["z"]], tol = 1e-4)
  expect_vector_equal(op$par[["b_sigma_Intercept"]],
                      fe$sigma[["(Intercept)"]], tol = 1e-4)
  expect_vector_equal(op$par[["b_sigma[1]"]], fe$sigma[["x"]], tol = 1e-4)
  # with flat priors brms's log posterior IS the log likelihood
  expect_lt(abs(op$value - as.numeric(logLik(fit))), 1e-4)
})

test_that("mo() ML matches brms's monotonic likelihood (vignette model)", {
  skip_unless_brms_fit()

  set.seed(3)
  dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
  dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)

  flat <- c(brms::prior("", class = "Intercept"),
            brms::prior("", class = "b"),
            brms::prior("", class = "sigma"))
  bform <- brms::bf(y ~ mo(inc) + z)
  code <- brms::make_stancode(bform, data = dm, family = gaussian(),
                              prior = flat)
  sdat <- brms::make_standata(bform, data = dm, family = gaussian(),
                              prior = flat)
  mod <- rstan::stan_model(model_code = code)
  sf <- rstan::sampling(mod, data = sdat, chains = 0)

  fit <- frm(bf(y ~ mo(inc) + z) + gaussian(), data = dm)
  fe <- fixef(fit)$mu
  simplex <- exp(c(0, fit$estimates$zeta1))
  simplex <- simplex / sum(simplex)
  # brms centers X but not Xmo, so only the intercept needs translating
  bpars <- list(
    b = array(fe[["z"]], 1),
    Intercept = as.numeric(fe[["(Intercept)"]] + mean(dm$z) * fe[["z"]]),
    bsp = array(fe[["moinc"]], 1),
    simo_1 = simplex,
    sigma = exp(fixef(fit)$sigma[[1]])
  )
  lp <- rstan::log_prob(sf, rstan::unconstrain_pars(sf, bpars),
                        adjust_transform = FALSE, gradient = FALSE)
  # the only prior left is the flat Dirichlet on the simplex, whose
  # normalizing constant is lgamma(D)
  expect_lt(abs(lp - as.numeric(logLik(fit)) - lgamma(length(simplex))), 1e-6)

  # and our estimate is brms's optimum, not merely a point on its surface
  op <- rstan::optimizing(mod, data = sdat, init = bpars, hessian = FALSE,
                          as_vector = TRUE, seed = 1)
  expect_lt(op$value - lp, 1e-4)
  expect_vector_equal(op$par[grep("^simo_1", names(op$par))], simplex,
                      tol = 1e-4)
  expect_vector_equal(op$par[["bsp[1]"]], fe[["moinc"]], tol = 1e-4)
  expect_vector_equal(op$par[["b_Intercept"]], fe[["(Intercept)"]],
                      tol = 1e-4)
})

test_that("R-side autocorrelation matches brms's time-series indexing", {
  skip_unless_brms()
  set.seed(4)
  d <- expand.grid(week = 1:4, subj = factor(1:6))
  d$x <- rnorm(nrow(d))
  d$y <- rnorm(nrow(d))
  # brms sorts by (gr, time) and reports contiguous blocks; our frame
  # keeps the data order and carries one row-index vector per pattern,
  # so the comparable quantities are the group count and the group sizes
  for (tm in c("ar(week, subj, cov = TRUE)", "ma(week, subj, cov = TRUE)",
               "arma(week, subj, cov = TRUE)", "cosy(week, subj)",
               "unstr(week, subj)")) {
    form <- stats::as.formula(paste("y ~ x +", tm))
    sd_b <- brms_standata(brms::bf(form), data = d, family = gaussian())
    ac <- frm(form, data = d, family = gaussian(),
              dry_run = "frame")$autocor[[1L]]
    expect_equal(ac$n_groups, as.integer(sd_b$N_tg), info = tm)
    sizes <- sort(unlist(lapply(ac$patterns, function(p) {
      rep(p$k, p$G)
    })))
    expect_equal(sizes, sort(as.integer(sd_b$nobs_tg)), info = tm)
  }
  # unstr is the one structure brms indexes by TIME LEVEL rather than by
  # position, so its time bookkeeping is directly comparable
  sd_u <- brms_standata(brms::bf(y ~ x + unstr(week, subj)), data = d,
                        family = gaussian())
  fr <- frm(y ~ x + unstr(week, subj), data = d, family = gaussian(),
            dry_run = "frame")
  ac <- fr$autocor[[1L]]
  expect_equal(ac$d, as.integer(sd_u$n_unique_t))
  expect_equal(length(fr$par_template$thetaac),
               as.integer(sd_u$n_unique_cortime))
  # brms's Jtime_tg is one row per group holding that group's time
  # levels; ours is the same set, read out of the pattern gather
  ours <- sort(unique(as.vector(vapply(ac$patterns, function(p) {
    sort(unique((p$gather - 1L) %% ac$d + 1L))
  }, integer(ac$d)))))
  expect_equal(ours, sort(unique(as.vector(sd_u$Jtime_tg))))

  # brms refuses the same aterm combinations at the Stan-code stage
  expect_error(
    suppressMessages(brms::make_stancode(
      brms::bf(y | weights(x) ~ 1 + ar(week, subj, cov = TRUE)),
      data = transform(d, x = abs(d$x) + 1), family = gaussian())),
    "Invalid addition arguments")
  expect_error(
    frm(bf(y | weights(w) ~ 1 + ar(week, subj, cov = TRUE)) + gaussian(),
        data = transform(d, w = abs(d$x) + 1)),
    "cannot be combined with a residual correlation term")
})
