# The Dirichlet draw that dirichlet(), brmsfamily("dirichlet2") and
# dirichlet_multinomial() share (simplex_rdirichlet()). Drawn as
# normalized gamma variates on the probability scale, every variate of
# a row underflows to 0 at a small concentration and the row is 0 / 0;
# the draw is on the log scale. The references are the analytic
# Dirichlet moments, judged against Monte Carlo standard errors the run
# measures.

rd_cats <- c("a", "b", "c")
rd_eta <- c(0.3, -0.2)
rd_p <- exp(c(0, rd_eta)) / sum(exp(c(0, rd_eta)))

rd_dirichlet <- function(phi, n) {
  fam <- frmtmb:::fam_dirichlet_impl(rd_cats, 1L)
  dp <- c(stats::setNames(lapply(rd_eta, rep, n), fam$primary_dpars),
          list(phi = rep(phi, n)))
  fam$sim(dp, list(), n)
}

rd_dirichlet2 <- function(phi, n) {
  fam <- frmtmb:::fam_dirichlet2_impl(rd_cats)
  dp <- stats::setNames(lapply(phi * rd_p, rep, n), fam$primary_dpars)
  fam$sim(dp, list(), n)
}

test_that("a small concentration leaves every Dirichlet row a composition", {
  set.seed(5)
  n <- 20000
  for (phi in c(0.01, 1e-3)) {
    for (draw in list(rd_dirichlet, rd_dirichlet2)) {
      x <- draw(phi, n)
      # rellib-r10: 18 rows NaN at phi 0.01 and 9570 at 1e-3
      expect_identical(sum(!is.finite(rowSums(x))), 0L,
                       label = paste("phi", phi))
      expect_lt(max(abs(rowSums(x) - 1)), 64 * .Machine$double.eps)
    }
  }
  fam <- dirichlet_multinomial(K = 3)
  x <- fam$sim(list(mu2 = rep(rd_eta[1], n), mu3 = rep(rd_eta[2], n),
                    phi = rep(1e-3, n)), list(trials = rep(12, n)), n)
  expect_true(all(rowSums(x) == 12))
})

test_that("the draws keep the Dirichlet's means and covariances", {
  set.seed(6)
  n <- 1e5
  z <- stats::qnorm(1 - 1e-6 / 2)
  for (phi in c(0.01, 0.5, 20)) {
    for (draw in list(rd_dirichlet, rd_dirichlet2)) {
      x <- draw(phi, n)
      lab <- paste("phi", phi)
      for (k in 1:3) {
        expect_lt(abs(mean(x[, k]) - rd_p[k]),
                  z * stats::sd(x[, k]) / sqrt(n),
                  label = paste(lab, "mean", k))
        v <- rd_p[k] * (1 - rd_p[k]) / (phi + 1)
        d2 <- (x[, k] - rd_p[k])^2
        expect_lt(abs(mean(d2) - v), z * stats::sd(d2) / sqrt(n),
                  label = paste(lab, "variance", k))
      }
      cv <- -rd_p[2] * rd_p[3] / (phi + 1)
      pr <- (x[, 2] - rd_p[2]) * (x[, 3] - rd_p[3])
      expect_lt(abs(mean(pr) - cv), z * stats::sd(pr) / sqrt(n),
                label = paste(lab, "covariance"))
    }
  }
})

test_that("predict() runs on a dirichlet_multinomial fit with a small phi", {
  # the fit of the review's case: phi near 0.02, where predict()'s
  # parameter draws reach smaller values; it stopped with "NA in
  # probability vector"
  set.seed(41)
  n <- 200L
  dd <- data.frame(x = stats::rnorm(n), N = sample(5:30, n, TRUE))
  E <- cbind(0, 0.3 + 0.5 * dd$x, -0.2 - 0.4 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  Y <- t(vapply(seq_len(n), function(i) {
    lg <- log(stats::rgamma(3, P[i, ] * 0.02 + 1)) +
      log(stats::runif(3)) / (P[i, ] * 0.02)
    p <- exp(lg - max(lg))
    as.numeric(stats::rmultinom(1, dd$N[i], p / sum(p)))
  }, numeric(3)))
  colnames(Y) <- rd_cats
  dd$Y <- Y
  fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = dd)
  phi <- exp(fit$estimates$betad[[1]])
  expect_lt(phi, 0.1)
  set.seed(7)
  pr <- predict(fit)
  expect_identical(dim(pr), c(n, 4L, 3L))
  expect_true(all(is.finite(pr)))
  s <- simulate(fit, nsim = 3, seed = 8)
  expect_true(all(vapply(s, function(m) all(rowSums(m) == dd$N), NA)))
})
