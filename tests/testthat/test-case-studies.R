# Pins the numeric cross-checks of vignette("case-studies"). Each test
# reproduces one section's fit and compares it with the same reference
# the vignette compares it with, so the vignette's claims cannot rot
# silently.

# tabular method for the additive relationship matrix of a pedigree
# sorted parents-first
build_A <- function(ped) {
  n <- nrow(ped)
  A <- diag(n)
  for (i in seq_len(n)) {
    s <- ped$sire[i]
    d <- ped$dam[i]
    if (!is.na(s)) {
      A[i, i] <- 1 + 0.5 * A[s, d]
      for (j in seq_len(i - 1)) {
        A[i, j] <- A[j, i] <- 0.5 * (A[j, s] + A[j, d])
      }
    }
  }
  dimnames(A) <- list(as.character(ped$id), as.character(ped$id))
  A
}

build_pedigree <- function(nsire, ndam_per, noff) {
  nfound <- nsire + nsire * ndam_per
  n <- nfound + nsire * ndam_per * noff
  ped <- data.frame(id = seq_len(n), sire = NA_integer_, dam = NA_integer_)
  k <- nfound
  for (s in seq_len(nsire)) {
    for (j in seq_len(ndam_per)) {
      dam <- nsire + ndam_per * (s - 1L) + j
      for (o in seq_len(noff)) {
        k <- k + 1L
        ped$sire[k] <- s
        ped$dam[k] <- dam
      }
    }
  }
  ped
}

test_that("the animal model matches the closed-form REML likelihood", {
  skip_on_cran()
  ped <- build_pedigree(12, 2, 6)
  n <- nrow(ped)
  A <- build_A(ped)
  set.seed(9)
  dat <- data.frame(id = factor(ped$id, levels = ped$id),
                    sex = factor(rep(c("f", "m"), length.out = n)))
  dat$phen <- 5 + 0.4 * (dat$sex == "m") +
    as.numeric(t(chol(A)) %*% rnorm(n)) + rnorm(n, 0, 1)

  fit <- frm(bf(phen ~ sex + (1 | gr(id, cov = A))) + gaussian(),
             data = dat, data2 = list(A = A), REML = TRUE,
             control = frmtmb_control(check_olre = "ignore"))

  # y ~ N(X beta, sd_a^2 A + sigma^2 I): the restricted likelihood in
  # closed form, maximized without any part of frmtmb
  X <- stats::model.matrix(~ sex, dat)
  y <- dat$phen
  nll <- function(p) {
    V <- exp(p[1]) * A + exp(p[2]) * diag(n)
    Vi <- solve(V)
    M <- crossprod(X, Vi) %*% X
    b <- solve(M, crossprod(X, Vi) %*% y)
    r <- y - X %*% b
    as.numeric(0.5 * ((n - ncol(X)) * log(2 * pi) +
                        determinant(V)$modulus + determinant(M)$modulus +
                        crossprod(r, Vi %*% r)))
  }
  op <- stats::optim(c(0, 0), nll, method = "BFGS",
                     control = list(reltol = 1e-13))

  sd_a <- sqrt(VarCorr(fit)[[1]][1, 1])
  expect_equal(sd_a, sqrt(exp(op$par[1])), tolerance = 1e-5)
  expect_equal(sigma(fit), sqrt(exp(op$par[2])), tolerance = 1e-5)
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-5)
})

test_that("the multi-trait animal model reads the pedigree", {
  skip_on_cran()
  ped <- build_pedigree(12, 2, 6)
  n <- nrow(ped)
  A <- build_A(ped)
  set.seed(4)
  G <- matrix(c(1.0, 0.6, 0.6, 0.8), 2, 2)
  U <- t(chol(A)) %*% matrix(rnorm(n * 2), n, 2) %*% chol(G)
  long <- data.frame(
    id = factor(rep(ped$id, times = 2), levels = ped$id),
    trait = factor(rep(c("y1", "y2"), each = n)),
    value = c(3 + U[, 1] + rnorm(n, 0, 0.7),
              1 + U[, 2] + rnorm(n, 0, 0.9)))

  # the long format is one multi-trait spelling: one random coefficient
  # per trait, so the block covariance IS the genetic matrix G and
  # cov = A makes the whole covariance G %x% A. Since v0.32 the
  # mvbf() + |ID| + gr(cov =) spelling gives the same joint density
  # (test-id-kron.R holds the equivalence).
  fmv <- frm(bf(value ~ 0 + trait + (0 + trait | gr(id, cov = A)),
                sigma ~ 0 + trait) + gaussian(),
             data = long, data2 = list(A = A))
  Gh <- VarCorr(fmv)[[1]]
  expect_equal(cov2cor(Gh)[1, 2], 0.6 / sqrt(1.0 * 0.8), tolerance = 0.05)
  expect_equal(unname(sqrt(diag(Gh))), c(1, sqrt(0.8)), tolerance = 0.2)
  expect_equal(unname(exp(fixef(fmv)$sigma)), c(0.7, 0.9), tolerance = 0.15)

  # a known covariance can be accepted, ignored, and still converge.
  # Refitting with the identity has to move the log likelihood, or A
  # never entered the likelihood at all.
  I <- diag(n)
  dimnames(I) <- dimnames(A)
  fmv_id <- frm(bf(value ~ 0 + trait + (0 + trait | gr(id, cov = I)),
                   sigma ~ 0 + trait) + gaussian(),
                data = long, data2 = list(I = I))
  expect_gt(as.numeric(logLik(fmv)) - as.numeric(logLik(fmv_id)), 5)
})

test_that("the phylogenetic mixed model matches PGLS with Pagel's lambda", {
  skip_on_cran()
  skip_if_not_installed("ape")
  skip_if_not_installed("nlme")
  set.seed(7)
  tree <- ape::rcoal(60)
  tree$tip.label <- paste0("sp", 1:60)
  tree$edge.length <- tree$edge.length /
    max(ape::node.depth.edgelength(tree))
  A_phy <- ape::vcv(tree)
  d <- data.frame(sp = factor(tree$tip.label, levels = tree$tip.label),
                  x = rnorm(60))
  d$y <- 0.5 + 0.8 * d$x +
    as.numeric(t(chol(A_phy)) %*% rnorm(60)) + rnorm(60, 0, 0.6)

  fphy <- frm(bf(y ~ x + (1 | gr(sp, cov = A_phy))) + gaussian(),
              data = d, data2 = list(A_phy = A_phy),
              control = frmtmb_control(check_olre = "ignore"))
  g <- nlme::gls(y ~ x, data = d, method = "ML",
                 correlation = ape::corPagel(0.5, phy = tree, form = ~sp))

  sd_p <- sqrt(VarCorr(fphy)[[1]][1, 1])
  lam <- sd_p^2 / (sd_p^2 + sigma(fphy)^2)
  lam_gls <- as.numeric(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE))
  expect_equal(unname(fixef(fphy)$mu), unname(coef(g)), tolerance = 1e-4)
  expect_equal(lam, lam_gls, tolerance = 1e-4)
  expect_equal(sd_p^2 + sigma(fphy)^2, g$sigma^2, tolerance = 1e-4)
  expect_lt(abs(as.numeric(logLik(fphy)) - as.numeric(logLik(g))), 1e-4)
})

test_that("se() meta-analysis matches metafor::rma", {
  skip_on_cran()
  skip_if_not_installed("metafor")
  # Colditz et al. (1994) BCG vaccine trials
  bcg <- data.frame(
    tpos = c(4, 6, 3, 62, 33, 180, 8, 505, 29, 17, 186, 5, 27),
    tneg = c(119, 300, 228, 13536, 5036, 1361, 2537, 87886, 7470, 1699,
             50448, 2493, 16886),
    cpos = c(11, 29, 11, 248, 47, 372, 10, 499, 45, 65, 141, 3, 29),
    cneg = c(128, 274, 209, 12619, 5761, 1079, 619, 87892, 7232, 1600,
             27197, 2338, 17825),
    ablat = c(44, 55, 42, 52, 13, 44, 19, 13, 27, 42, 18, 33, 33))
  bcg$yi <- log((bcg$tpos / (bcg$tpos + bcg$tneg)) /
                  (bcg$cpos / (bcg$cpos + bcg$cneg)))
  bcg$vi <- 1 / bcg$tpos - 1 / (bcg$tpos + bcg$tneg) +
    1 / bcg$cpos - 1 / (bcg$cpos + bcg$cneg)
  bcg$sei <- sqrt(bcg$vi)
  bcg$study <- factor(seq_len(nrow(bcg)))

  fmeta <- frm(bf(yi | se(sei) ~ 1 + (1 | study)) + gaussian(),
               data = bcg, REML = TRUE)
  rr <- metafor::rma(bcg$yi, bcg$vi, method = "REML")
  expect_equal(unname(fixef(fmeta)$mu), as.numeric(rr$beta),
               tolerance = 1e-5)
  expect_equal(sqrt(VarCorr(fmeta)[[1]][1, 1]), sqrt(rr$tau2),
               tolerance = 1e-5)
  # metafor forms the standard error as (X'WX)^-1 at the REML tau;
  # frmtmb reads it off the joint Hessian, so they are close, not equal
  expect_lt(abs(sqrt(vcov(fmeta)[1, 1]) / rr$se - 1), 0.01)

  freg <- frm(bf(yi | se(sei) ~ ablat + (1 | study)) + gaussian(),
              data = bcg, REML = TRUE)
  rr2 <- metafor::rma(bcg$yi, bcg$vi, mods = ~ bcg$ablat, method = "REML")
  expect_equal(unname(fixef(freg)$mu), unname(coef(rr2)), tolerance = 1e-4)
  expect_equal(sqrt(VarCorr(freg)[[1]][1, 1]), sqrt(rr2$tau2),
               tolerance = 1e-3)
})

test_that("mo() reproduces the saturated factor fit on monotone data", {
  skip_on_cran()
  set.seed(21)
  n <- 400
  L <- 6
  inc <- sample(seq_len(L), n, TRUE)
  shape <- c(0, 0.15, 0.3, 0.75, 0.9, 1)
  dmo <- data.frame(income = factor(inc, ordered = TRUE), z = rnorm(n))
  dmo$ls <- 5 + 2 * shape[inc] + 0.4 * dmo$z + rnorm(n, 0, 0.8)

  fmo <- frm(bf(ls ~ mo(income) + z) + gaussian(), data = dmo)
  fsat <- frm(bf(ls ~ income + z) + gaussian(), data = dmo)
  nd <- data.frame(income = factor(seq_len(L), ordered = TRUE), z = 0)
  p <- predict(fmo, newdata = nd)

  expect_equal(unname(p), unname(predict(fsat, newdata = nd)),
               tolerance = 1e-4)
  expect_equal(as.numeric(logLik(fmo)), as.numeric(logLik(fsat)),
               tolerance = 1e-6)
  expect_identical(attr(logLik(fmo), "df"), attr(logLik(fsat), "df"))
  # the recovered shape tracks the simulated step sizes
  expect_lt(max(abs((p - p[1]) / (p[L] - p[1]) - shape)), 0.12)
})

test_that("location-scale smooths agree with mgcv's gaulss", {
  skip_on_cran()
  set.seed(22)
  dls <- data.frame(x = runif(400))
  dls$y <- rnorm(400, sin(2 * pi * dls$x),
                 exp(-1 + 1.2 * cos(2 * pi * dls$x)))
  fls <- frm(bf(y ~ s(x, k = 10), sigma ~ s(x, k = 10)) + gaussian(),
             data = dls)
  gm <- mgcv::gam(list(y ~ s(x, k = 10), ~ s(x, k = 10)), data = dls,
                  family = mgcv::gaulss(b = 0), method = "ML")
  nd <- data.frame(x = seq(0.05, 0.95, length.out = 6))
  pg <- stats::predict(gm, newdata = nd, type = "response")
  sg <- predict(fls, newdata = nd, dpar = "sigma", type = "response")
  # the two packages pick smoothing parameters by different criteria,
  # so the curves are close rather than identical
  expect_lt(max(abs(predict(fls, newdata = nd) - pg[, 1])), 0.02)
  expect_lt(max(abs(sg * pg[, 2] - 1)), 0.03)
})

test_that("cs() under sratio equals the binomial regression decomposition", {
  skip_on_cran()
  set.seed(51)
  n <- 500
  xo <- rnorm(n)
  th <- c(-1.6, 0, 1.6)
  eff <- c(0.3, 0.6, 0.9)
  Y <- integer(n)
  for (i in seq_len(n)) {
    cp <- c(stats::plogis(th - eff * xo[i]), 1)
    Y[i] <- sample(1:4, 1, prob = diff(c(0, cp)))
  }
  dord <- data.frame(y = factor(Y, ordered = TRUE), x = xo)

  fcs <- frm(bf(y ~ cs(x)) + sratio(), data = dord)
  glm_fits <- lapply(1:3, function(k) {
    sub <- dord[as.integer(dord$y) >= k, ]
    stats::glm(as.integer(as.integer(sub$y) == k) ~ x, data = sub,
               family = stats::binomial())
  })
  b_glm <- vapply(glm_fits, function(g) stats::coef(g)[["x"]], 0)
  ll_glm <- sum(vapply(glm_fits, function(g) as.numeric(logLik(g)), 0))

  # the sequential model factorizes into one binary fit per threshold;
  # frmtmb's cs coefficients carry the opposite sign convention
  expect_equal(unname(confint(fcs)[paste0("bcs2_", 1:3), "est"]), -b_glm,
               tolerance = 1e-4)
  expect_lt(abs(as.numeric(logLik(fcs)) - ll_glm), 1e-6)
})

test_that("the growth mixture recovers class-specific slopes", {
  skip_on_cran()
  set.seed(31)
  nid <- 80
  nt <- 5
  cls <- rbinom(nid, 1, 0.4)
  id <- rep(seq_len(nid), each = nt)
  dg <- data.frame(id = factor(id), time = rep(0:(nt - 1), nid))
  b0 <- 2 + rnorm(nid, 0, 0.5)
  b1 <- c(0.2, 1.2)[cls + 1]
  dg$y <- b0[id] + b1[id] * dg$time + rnorm(nid * nt, 0, 0.6)

  fgmm <- frm(bf(y ~ time + (1 | id)) +
                mixture(gaussian(), gaussian(), groups = ~id),
              data = dg)
  slopes <- sort(c(fixef(fgmm)$mu1[["time"]], fixef(fgmm)$mu2[["time"]]))
  expect_equal(slopes, c(0.2, 1.2), tolerance = 0.15)

  assigned <- max.col(mixture_probs(fgmm)) - 1L
  # label switching is inherent to a mixture, so score both labelings
  expect_gt(max(mean(assigned == cls), mean(assigned != cls)), 0.95)
  expect_lt(AIC(fgmm), AIC(frm(bf(y ~ time + (1 | id)) + gaussian(),
                               data = dg)))
})

test_that("mi(sdx) corrects attenuation from measurement error", {
  skip_on_cran()
  set.seed(3)
  n <- 300
  z <- rnorm(n)
  x_true <- rnorm(n, 0.5 * z, 1)
  su <- 0.6
  dme <- data.frame(x = x_true + rnorm(n, 0, su), z = z, su = su)
  dme$y <- 1 + 1 * x_true + 0.3 * z + rnorm(n, 0, 0.7)

  fme <- frm(bf(y ~ mi(x) + z) + gaussian() +
               bf(x | mi(su) ~ z) + gaussian(), data = dme)
  naive <- stats::coef(stats::lm(y ~ x + z, data = dme))[["x"]]
  corrected <- fixef(fme)$y_mu[["mix"]]

  expect_lt(naive, 0.9)
  expect_equal(corrected, 1, tolerance = 0.12)
  # the correction costs precision: the latent values are estimated
  expect_gt(sqrt(diag(vcov(fme)))[["y_mix"]],
            summary(stats::lm(y ~ x + z, data = dme))$coefficients["x", 2])
})
