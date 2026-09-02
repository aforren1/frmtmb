# Cluster-robust (sandwich) covariance: vcov_cluster(), cluster_scores(),
# and the vcov= feed-through in confint()/hypothesis()/summary().

test_that("per-cluster scores add up to the gradient at the optimum", {
  set.seed(11)
  G <- 25
  dd <- data.frame(g = factor(rep(seq_len(G), each = 8)),
                   x = rnorm(G * 8))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G, 0, 0.8)[dd$g] + rnorm(G * 8)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)

  S <- cluster_scores(fit, ~ g)
  expect_equal(dim(S), c(G, length(fit$opt$par)))
  expect_equal(rownames(S), levels(dd$g))
  expect_equal(colnames(S), frmtmb:::outer_par_names(fit))

  # the scores are the pieces of the gradient the optimizer drove to
  # zero: they must reproduce it exactly, and it must be near zero
  expect_equal(colSums(S), -drop(fit$obj$gr(fit$opt$par)),
               ignore_attr = TRUE, tolerance = 1e-7)
  expect_lt(max(abs(colSums(S))), 1e-3)
})

test_that("one cluster per row reproduces the classical glm sandwich", {
  skip_if_not_installed("sandwich")
  set.seed(12)
  n <- 150
  dd <- data.frame(x = rnorm(n), z = runif(n))
  dd$y <- rpois(n, exp(0.4 + 0.6 * dd$x - 0.3 * dd$z))
  fit <- frm(bf(y ~ x + z) + poisson(), data = dd)
  ref <- glm(y ~ x + z, poisson(), dd)

  expect_equal(unname(fixef(fit)$mu), unname(coef(ref)), tolerance = 1e-6)

  # with no random effects and one cluster per row, our per-cluster
  # scores ARE sandwich::estfun()'s per-observation scores
  S <- cluster_scores(fit, factor(seq_len(n)))
  expect_equal(unname(S), unname(sandwich::estfun(ref)),
               tolerance = 1e-5)

  V <- vcov_cluster(fit, factor(seq_len(n)), type = "CR0")
  expect_equal(unname(V), unname(sandwich::vcovHC(ref, type = "HC0")),
               ignore_attr = TRUE, tolerance = 1e-5)
  # HC1 is the same n / (n - p) rescaling CR1p applies here
  V1 <- vcov_cluster(fit, factor(seq_len(n)), type = "CR1p")
  expect_equal(unname(V1), unname(sandwich::vcovHC(ref, type = "HC1")),
               ignore_attr = TRUE, tolerance = 1e-5)
})

test_that("clustered LMM matches clubSandwich on the matched lme4 fit", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("clubSandwich")
  set.seed(13)
  G <- 30
  m <- 7
  dd <- data.frame(g = factor(rep(seq_len(G), each = m)),
                   x = rnorm(G * m))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G, 0, 0.9)[dd$g] +
    rnorm(G * m, 0, rep(runif(G, 0.4, 2), each = m))
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)
  ref <- lme4::lmer(y ~ x + (1 | g), dd, REML = FALSE)

  # same estimates first, or nothing below means anything
  expect_equal(unname(fixef(fit)$mu), unname(lme4::fixef(ref)),
               tolerance = 1e-5)

  cs0 <- as.matrix(clubSandwich::vcovCR(ref, cluster = dd$g,
                                        type = "CR0"))
  cs1 <- as.matrix(clubSandwich::vcovCR(ref, cluster = dd$g,
                                        type = "CR1"))

  # clubSandwich CONDITIONS on the variance parameters: its bread is the
  # mean-model information alone. Rebuild that form from our pieces and
  # it agrees exactly.
  S <- cluster_scores(fit, ~ g)
  info <- solve(vcov(fit, full = TRUE))
  bi <- seq_along(fixef(fit)$mu)
  A <- solve(info[bi, bi, drop = FALSE])
  cond0 <- A %*% crossprod(S[, bi, drop = FALSE]) %*% A
  # relative, entry by entry: what is left is the two optimizers'
  # disagreement about theta, not the estimator
  expect_lt(max(abs(unname(cond0) / unname(cs0) - 1)), 1e-3)
  expect_lt(max(abs(G / (G - 1) * unname(cond0) / unname(cs1) - 1)),
            1e-3)

  # the shipped estimator sandwiches the WHOLE outer parameter vector,
  # so it differs by the cost of estimating theta: same ballpark, not
  # the same number (documented on the vcov_cluster() page). Its
  # fixed-effect block carries log(sigma) too, as vcov() does.
  V0 <- vcov_cluster(fit, ~ g, type = "CR0")
  expect_equal(dim(V0), dim(vcov(fit)))
  expect_lt(max(abs(sqrt(diag(V0)[bi]) / sqrt(diag(cs0)) - 1)), 0.05)
  expect_equal(attr(V0, "nclusters"), G)
  expect_equal(attr(V0, "df"), G - 1L)
})

test_that("small-sample factors and the surface behave", {
  set.seed(14)
  G <- 20
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   x = rnorm(G * 6))
  dd$y <- 1 + 0.4 * dd$x + rnorm(G, 0, 0.6)[dd$g] + rnorm(G * 6)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)
  N <- nrow(dd)
  p <- length(fit$opt$par)

  V0 <- vcov_cluster(fit, ~ g, type = "CR0")
  expect_equal(vcov_cluster(fit, ~ g, type = "CR1"),
               V0 * G / (G - 1), ignore_attr = TRUE)
  expect_equal(vcov_cluster(fit, ~ g, type = "CR1p"),
               V0 * G / (G - p), ignore_attr = TRUE)
  expect_equal(vcov_cluster(fit, ~ g, type = "CR1S"),
               V0 * G * (N - 1) / ((G - 1) * (N - p)),
               ignore_attr = TRUE)

  # symmetric, positive semidefinite, and named like vcov()
  expect_equal(V0, t(V0), ignore_attr = TRUE)
  expect_gte(min(eigen(V0, only.values = TRUE)$values), -1e-10)
  expect_equal(rownames(V0), rownames(vcov(fit)))
  expect_equal(rownames(vcov_cluster(fit, ~ g, full = TRUE)),
               frmtmb:::outer_par_names(fit))

  # vcov(cluster =) is the sandwich::vcovCL spelling and forwards here
  expect_equal(vcov(fit, cluster = ~ g, type = "CR1"),
               vcov_cluster(fit, ~ g, type = "CR1"))
  # a plain factor and a variable name resolve the same way
  expect_equal(vcov_cluster(fit, dd$g), V0)
  expect_equal(vcov_cluster(fit, "g"), V0)

  # CR2 is refused, not faked
  expect_error(vcov_cluster(fit, ~ g, type = "CR2"), "Bell-McCaffrey")
  expect_error(vcov_cluster(fit, ~ g, type = "CR3"), "Bell-McCaffrey")
})

test_that("confint, hypothesis and summary accept the robust matrix", {
  set.seed(15)
  G <- 24
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   x = rnorm(G * 6))
  dd$y <- 1 + 0.4 * dd$x + rnorm(G, 0, 0.6)[dd$g] +
    rnorm(G * 6, 0, rep(runif(G, 0.3, 2), each = 6))
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)
  V <- vcov_cluster(fit, ~ g, type = "CR1", full = TRUE)

  ci_m <- confint(fit, parm = "x")
  ci_r <- confint(fit, parm = "x", vcov = V)
  expect_equal(ci_m[, "est"], ci_r[, "est"])
  # robust intervals here are wider: the generator has cluster-specific
  # error scales the model does not describe
  expect_gt(diff(ci_r[1, 1:2]), diff(ci_m[1, 1:2]))
  # a t(G - 1) reference, not a normal one
  se <- sqrt(diag(V))[match("x", frmtmb:::outer_par_names(fit))]
  expect_equal(unname(ci_r[1, "upr"] - ci_r[1, "est"]),
               unname(qt(0.975, G - 1) * se))

  # a function of the fit is accepted as well
  expect_equal(confint(fit, parm = "x",
                       vcov = function(f) {
                         vcov_cluster(f, ~ g, "CR1", full = TRUE)
                       }),
               ci_r)

  hy_m <- hypothesis(fit, "x = 0")
  hy_r <- hypothesis(fit, "x = 0", vcov = V)
  expect_equal(hy_m$estimate, hy_r$estimate)
  expect_gt(hy_r$se, hy_m$se)

  sm <- summary(fit, vcov = V)
  expect_true("t value" %in% colnames(sm$coefficients$mu))
  expect_equal(unname(sm$coefficients$mu[, "Std. Error"]),
               unname(sqrt(diag(V))[1:2]))

  # the fixed-effect block alone is not enough, and says so
  expect_error(confint(fit, vcov = vcov_cluster(fit, ~ g)),
               "full = TRUE")
  expect_error(confint(fit, parm = "x", method = "profile", vcov = V),
               "method = 'wald' only")
})

test_that("guards refuse every structure the likelihood does not factor over", {
  set.seed(16)
  dd <- data.frame(a = factor(rep(1:6, each = 12)),
                   b = factor(rep(1:12, 6)),
                   x = rnorm(72))
  dd$y <- rnorm(72, 1 + 0.3 * dd$x + rnorm(6, 0, .5)[dd$a] +
                  rnorm(12, 0, .5)[dd$b])

  # crossed random effects: b spans every level of a
  f_cross <- frm(bf(y ~ x + (1 | a) + (1 | b)) + gaussian(), data = dd,
                 REML = FALSE)
  expect_error(vcov_cluster(f_cross, ~ a), "crosses `cluster`")

  # nested is fine, and coarser than the random effect is fine too
  dd2 <- data.frame(school = factor(rep(1:12, each = 8)))
  dd2$class <- factor(paste0(dd2$school, ".",
                             rep(rep(1:2, each = 4), 12)))
  dd2$x <- rnorm(96)
  dd2$y <- rnorm(96, 1 + 0.3 * dd2$x + rnorm(12, 0, .6)[dd2$school] +
                   rnorm(24, 0, .4)[dd2$class])
  f_nest <- frm(bf(y ~ x + (1 | school) + (1 | class)) + gaussian(),
                data = dd2, REML = FALSE)
  expect_silent(V <- vcov_cluster(f_nest, ~ school))
  expect_equal(attr(V, "nclusters"), 12L)
  # the finer factor does NOT contain the school effect
  expect_error(vcov_cluster(f_nest, ~ class), "crosses `cluster`")

  # a global smooth is one random effect over every row
  skip_if_not_installed("mgcv")
  f_sm <- frm(bf(y ~ s(x) + (1 | school)) + gaussian(), data = dd2,
              REML = FALSE)
  expect_error(vcov_cluster(f_sm, ~ school), "crosses `cluster`")

  # REML, priors and a one-level cluster
  f_reml <- frm(bf(y ~ x + (1 | school)) + gaussian(), data = dd2,
                REML = TRUE)
  expect_error(vcov_cluster(f_reml, ~ school), "maximum-likelihood fit")
  f_pri <- frm(bf(y ~ x + (1 | school)) + gaussian(), data = dd2,
               REML = FALSE,
               prior = set_prior("normal(0, 1)", class = "b",
                                 coef = "x"))
  expect_error(vcov_cluster(f_pri, ~ school), "made with priors")
  expect_error(vcov_cluster(f_nest, rep(1, nrow(dd2))),
               "at least 2")
  expect_error(vcov_cluster(f_nest, dd2$school[1:10]),
               "fitted to 96 rows")
  expect_error(vcov_cluster(f_nest, ~ school + class),
               "more than one variable")
  bad <- as.character(dd2$school)
  bad[3] <- NA
  expect_error(vcov_cluster(f_nest, bad), "missing values")
})

test_that("group-level mixtures factor at their own grouping", {
  set.seed(21)
  G <- 30
  m <- 5
  gv <- factor(rep(seq_len(G), each = m))
  cls <- rbinom(G, 1, 0.4)
  dd <- data.frame(y = rnorm(G * m, c(0, 3)[cls[gv] + 1L], 1), g = gv)
  fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~ g),
             data = dd)

  # a masked-out mixture group contributes log(sum_k pi_k) = 0 exactly,
  # so the identity holds here as it does without a mixture
  S <- cluster_scores(fit, ~ g)
  expect_equal(nrow(S), G)
  expect_equal(colSums(S), -drop(fit$obj$gr(fit$opt$par)),
               ignore_attr = TRUE, tolerance = 1e-7)

  # coarser than the mixture grouping is still a factorization
  pairs <- factor(((as.integer(gv) - 1L) %/% 2L) + 1L)
  expect_equal(nrow(cluster_scores(fit, pairs)), G / 2)

  # finer is not: one class draw covers the whole group
  expect_error(cluster_scores(fit, factor(seq_len(G * m))),
               "mixture group spans more than one cluster")
})

test_that("mm() pooled levels are refused unless nested", {
  set.seed(22)
  n <- 120
  dd <- data.frame(x = rnorm(n),
                   g1 = factor(sample(letters[1:5], n, TRUE)),
                   g2 = factor(sample(letters[3:8], n, TRUE)),
                   cl = factor(rep(1:10, each = 12)))
  dd$y <- rnorm(n, 1 + 0.4 * dd$x)
  fit <- frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), data = dd,
             REML = FALSE)
  # every pooled level is loaded by rows from several clusters
  expect_error(vcov_cluster(fit, ~ cl), "crosses `cluster`")
})

test_that("robust SEs recover nominal coverage where model-based ones do not", {
  skip_on_cran()
  set.seed(17)
  G <- 40
  m <- 12
  n <- G * m
  gvar <- factor(rep(seq_len(G), each = m))
  beta_x <- 0.5
  nsim <- 60
  cover_m <- logical(nsim)
  cover_r <- logical(nsim)
  for (s in seq_len(nsim)) {
    # heteroskedastic clusters AND an unmodeled random slope on a
    # predictor that varies mostly BETWEEN clusters: the fitted
    # (1 | g) model is misspecified in exactly the way cluster-robust
    # inference is meant to survive. See dev/sandwich/probe-coverage.R
    # for the tuning - plain heteroskedasticity alone barely dents the
    # model-based interval, because sigma absorbs its average.
    x <- rep(rnorm(G, 0, 2), each = m) + rnorm(n)
    sd_g <- rep(runif(G, 0.2, 2.5), each = m)
    slope_g <- rep(rnorm(G, 0, 1), each = m)
    y <- 1 + (beta_x + slope_g) * x + rnorm(G, 0, 0.7)[gvar] +
      rnorm(n, 0, sd_g)
    dd <- data.frame(y = y, x = x, g = gvar)
    # one replicate in 60 stops just short of the gradient tolerance;
    # a loose-tolerance coverage count does not care
    fit <- suppressWarnings(try(
      frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE),
      silent = TRUE))
    if (inherits(fit, "try-error")) next
    ci_m <- confint(fit, parm = "x")
    V <- vcov_cluster(fit, ~ g, type = "CR1", full = TRUE)
    ci_r <- confint(fit, parm = "x", vcov = V)
    cover_m[s] <- ci_m[1, "lwr"] <= beta_x && beta_x <= ci_m[1, "upr"]
    cover_r[s] <- ci_r[1, "lwr"] <= beta_x && beta_x <= ci_r[1, "upr"]
  }
  # loose, seeded thresholds: the point is the direction and the size of
  # the gap, not a precise coverage number from 60 replicates
  expect_lt(mean(cover_m), 0.80)
  expect_gt(mean(cover_r), 0.88)
  expect_gt(mean(cover_r) - mean(cover_m), 0.15)
})
