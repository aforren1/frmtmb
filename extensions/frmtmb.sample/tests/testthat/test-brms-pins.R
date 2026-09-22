# Pins on the CONTENT behind brms's names, none of which routes through
# brms_par_labels(): each value under a name is compared with what the
# model itself computes at the same draw, by position. A mislabeled
# r_ level or coefficient, a b_ coefficient under another's name, a
# correlation matrix read in the wrong order or a sigma column left on
# the log scale fails here. dev/brmsnames-mutants.R records each of
# those mutations failing this file.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

bp_rel <- function(got, ref) {
  # the same floating-point operations reach both sides, in another
  # order at most; element by element, so a small value is not judged at
  # the scale of a large one
  expect_true(all(abs(got - ref) <= 64 * .Machine$double.eps * abs(ref)),
              info = paste("worst relative difference",
                           max(abs(got - ref) / abs(ref))))
}

bp_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      # K = 4, so brms's correlation order and R's column-major lower
      # triangle differ (they agree up to K = 3)
      set.seed(5)
      G <- 10
      n <- 300
      d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                      w = stats::rnorm(n), g = factor(rep(1:G, 30)))
      S <- matrix(c(1, .5, -.3, .2, .5, 1, .3, -.2, -.3, .3, 1, .1,
                    .2, -.2, .1, 1), 4)
      L <- t(chol(S * 0.4))
      U <- t(L %*% matrix(stats::rnorm(4 * G), 4))
      d$y <- 1 + d$x + U[d$g, 1] + U[d$g, 2] * d$x + U[d$g, 3] * d$z +
        U[d$g, 4] * d$w + stats::rnorm(n)
      fit <- suppressWarnings(frm(bf(y ~ x + (1 + x + z + w | g)),
                                  family = gaussian(), data = d))
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 8)))
      cache <<- list(d = d, fit = fit, ds = ds)
    }
    cache
  }
})

test_that("r_, b_ and sigma columns hold what the model reads there", {
  cs <- bp_case()
  ds <- cs$ds
  idx <- frmtmb.sample:::draws_par_index(cs$fit)
  re <- ranef(ds, summary = FALSE)$g
  raw <- as.matrix(ds$stanfit)
  for (i in c(1L, 25L, 50L)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
    # the ML accessors read the estimate vector by position
    M <- unclass(ranef(sh))[[1L]]
    cn <- gsub("[()]", "", colnames(M))
    for (lv in rownames(M)) {
      for (k in seq_along(cn)) {
        col <- paste0("r_g[", lv, ",", cn[k], "]")
        expect_identical(unname(ds$draws[i, col]), unname(M[lv, k]),
                         info = col)
        expect_identical(unname(re[i, lv, cn[k]]), unname(M[lv, k]),
                         info = col)
      }
    }
    fe <- fixef_by_dpar(sh)$mu
    for (nm in names(fe)) {
      col <- paste0("b_", gsub("[()]", "", nm))
      expect_identical(unname(ds$draws[i, col]), unname(fe[[nm]]),
                       info = col)
    }
    # sigma is brms's natural-scale column, the exponential of the
    # log sigma the sampler itself stored
    bp_rel(unname(ds$draws[i, "sigma"]), exp(unname(raw[i, "betad"])))
    bp_rel(unname(ds$draws[i, "sigma"]), sigma(sh))
  }
  expect_false("sigma" %in% rownames(fixef(ds)))
  expect_false("b_sigma_Intercept" %in% variables(ds))
})

test_that("VarCorr() correlations are in brms's order at K = 4", {
  cs <- bp_case()
  ds <- cs$ds
  idx <- frmtmb.sample:::draws_par_index(cs$fit)
  vc <- VarCorr(ds, summary = FALSE)$g$cor
  nm <- c("Intercept", "x", "z", "w")
  for (i in c(1L, 25L, 50L)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
    C <- stats::cov2cor(varcorr_matrices(sh)[[1L]])
    for (a in 1:4) {
      for (b in 1:4) {
        bp_rel(vc[i, nm[a], nm[b]], C[a, b])
      }
    }
  }
})

test_that("log_lik() refuses brms's pointwise and point-estimate forms", {
  cs <- bp_case()
  expect_error(log_lik(cs$ds, pointwise = TRUE),
               "log_lik(pointwise = TRUE) is not supported", fixed = TRUE)
  expect_error(log_lik(cs$ds, add_point_estimate = TRUE),
               "log_lik(add_point_estimate = TRUE) is not supported",
               fixed = TRUE)
  # the defaults, spelled out, are the ordinary matrix
  expect_identical(dim(log_lik(cs$ds, pointwise = FALSE)),
                   c(nrow(cs$ds$draws), nrow(cs$d)))
})

test_that("ranef() and coef() compute a reduced-rank block per draw", {
  skip_sampler()
  set.seed(21)
  G <- 12
  n <- 240
  d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                  g = factor(rep(seq_len(G), length.out = n)))
  f <- stats::rnorm(G)
  d$y <- 1 + d$x1 + f[d$g] * (0.8 * d$x1 - 0.5 * d$x2) + stats::rnorm(n)
  fit <- suppressWarnings(frm(bf(y ~ x1 + rr(x1 + x2 | g, d = 1)),
                              family = gaussian(), data = d))
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 4)))
  re <- ranef(ds, summary = FALSE)$g
  expect_false(anyNA(re))
  expect_false(anyNA(coef(ds)$g))
  idx <- frmtmb.sample:::draws_par_index(fit)
  for (i in c(1L, 50L)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
    M <- unclass(ranef(sh))[[1L]]
    colnames(M) <- gsub("[()]", "", colnames(M))
    bp_rel(re[i, rownames(M), colnames(M)], M)
  }
})

test_that("ranef() and coef() refuse draws with no group-level draws", {
  # laplace-shaped draws: the random effects are in the model and no r_
  # column is in the draws. Real frm_sample(laplace = TRUE) draws give the
  # same refusal (dev/brmsnames-log/probe-laplace.txt); without it the
  # fill dies on "subscript out of bounds"
  cs <- bp_case()
  ld <- cs$ds
  ld$draws <- ld$draws[, !startsWith(colnames(ld$draws), "r_"),
                       drop = FALSE]
  expect_error(ranef(ld), "no draws of the group-level coefficients of 'g'",
               fixed = TRUE)
  expect_error(coef(ld), "frm_sample(laplace = TRUE)", fixed = TRUE)
})

test_that("multivariate VarCorr() and bayes_R2() have a row per response", {
  skip_sampler()
  set.seed(3)
  n <- 150
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:10, 15)))
  d$y_a <- 1 + 0.5 * d$x + stats::rnorm(10)[d$g] + stats::rnorm(n)
  d$y2 <- 0.4 * d$x + stats::rnorm(n)
  fit <- suppressWarnings(frm(mvbf(bf(y_a ~ x + (1 | g)), bf(y2 ~ x),
                                   rescor = FALSE),
                              family = gaussian(), data = d))
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 3)))
  vc <- VarCorr(ds)
  expect_named(vc, c("g", "residual__"))
  expect_identical(rownames(vc$residual__$sd), c("ya", "y2"))
  expect_true(all(c("sigma_ya", "sigma_y2") %in% variables(ds)))
  r2 <- bayes_R2(ds)
  expect_identical(rownames(r2), c("R2ya", "R2y2"))
  expect_false(isTRUE(all.equal(r2[1L, "Estimate"], r2[2L, "Estimate"])))
  expect_identical(rownames(bayes_R2(ds, resp = "y2")), "R2y2")
  expect_error(bayes_R2(ds, resp = "y_a"), "Valid response variables")
})

test_that("a mixture's weights are brms's theta1 and theta2 on draws", {
  skip_sampler()
  # dev/brmsnames-rev2-mixture.R, data seed 64: a share near 0.88, whose
  # log ratio near 2 cannot pass for it
  set.seed(64)
  n <- 400
  z <- stats::runif(n) < 0.9
  d <- data.frame(x = stats::rnorm(n))
  d$y <- ifelse(z, -2 + 0.3 * d$x + stats::rnorm(n, 0, 0.8),
                2.5 + 0.3 * d$x + stats::rnorm(n, 0, 1.2))
  fit <- suppressWarnings(frm(bf(y ~ x),
                              family = mixture(gaussian(), gaussian()),
                              data = d))
  set.seed(5)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3)))
  cn <- colnames(ds$draws)
  expect_true(all(c("theta1", "theta2") %in% cn))
  expect_identical(cn[length(cn)], "lp__")
  # the sampler's own log ratio, read from the stanfit by position
  raw <- as.matrix(ds$stanfit)
  j <- grep("^betad", colnames(raw))[3L]
  p1 <- stats::plogis(unname(raw[, j]))
  expect_gt(mean(p1), 0.8)
  bp_rel(unname(ds$draws[, "theta1"]), p1)
  bp_rel(unname(ds$draws[, "theta2"]), 1 - p1)
  h <- hypothesis(ds, "theta1 = 0.5", class = NULL)$hypothesis
  bp_rel(h$Estimate, mean(p1) - 0.5)
  # every reader that hands a draw back to the model gets the log ratio
  idx <- frmtmb.sample:::draws_par_index(fit)
  for (i in c(1L, 50L, 100L)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
    bp_rel(unname(sh$estimates$betad[3L]), unname(raw[i, j]))
  }
  expect_true("theta2" %in% rownames(summary(ds)))
})

test_that("r_ labels a level repeats are suffixed as brms suffixes them", {
  skip_sampler()
  # dev/brmsnames-rev2-collide.R C4, data seed 52: `lvl 1` and `lvl.1`
  # are both r_gd[lvl.1,Intercept] in brms, which suffixes the later one
  set.seed(52)
  n <- 240
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(paste("lvl", rep(1:12, length.out = n))))
  d$y <- 1 + 0.4 * d$x + stats::rnorm(12, 0, 0.8)[as.integer(d$g)] +
    stats::rnorm(n, 0, 0.8)
  gd <- as.character(d$g)
  i1 <- which(gd == "lvl 1")
  gd[i1[c(TRUE, FALSE)]] <- "lvl.1"
  d$gd <- factor(gd)
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | gd)), family = gaussian(),
                              data = d))
  set.seed(3)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 3)))
  cn <- colnames(ds$draws)
  expect_false(anyDuplicated(cn) > 0L)
  expect_identical(dim(posterior::as_draws_df(ds))[2L], length(cn) + 3L)
  # each name holds its own level: the suffixed one is `lvl.1`
  re <- ranef(ds, summary = FALSE)$gd
  bp_rel(unname(ds$draws[, "r_gd[lvl.1,Intercept]"]),
         unname(re[, "lvl 1", "Intercept"]))
  bp_rel(unname(ds$draws[, "r_gd[lvl.1,Intercept]__1"]),
         unname(re[, "lvl.1", "Intercept"]))
  # the unsuffixed name reads its own column, the first level
  h <- hypothesis(ds, "r_gd[lvl.1,Intercept] = 0", class = NULL)
  bp_rel(h$hypothesis$Estimate, mean(re[, "lvl 1", "Intercept"]))
  # draws that carry a name twice are refused, not read at the first
  dup <- ds
  colnames(dup$draws)[colnames(dup$draws) ==
                        "r_gd[lvl.1,Intercept]__1"] <- "r_gd[lvl.1,Intercept]"
  expect_error(hypothesis(dup, "r_gd[lvl.1,Intercept] = 0", class = NULL),
               "more than one column")
})
