# Lane wt-correct: the data every probe uses, so brms and frmtmb are
# measured on the same rows. Each generator fixes its own seed.

correct_data_gauss <- function() {
  set.seed(101)
  n <- 60
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(letters[1:6], length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(6, 0, 0.7)[d$g] + rnorm(n)
  d
}

correct_data_ord <- function() {
  set.seed(102)
  n <- 120
  d <- data.frame(x = rnorm(n))
  eta <- 0.8 * d$x
  u <- stats::rlogis(n)
  d$ord <- factor(cut(eta + u, c(-Inf, -0.5, 0.7, Inf), labels = FALSE),
                  ordered = TRUE)
  d$cat <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  d$y <- d$x + rnorm(n)
  d
}

correct_data_multinom <- function() {
  set.seed(103)
  n <- 60
  d <- data.frame(x = rnorm(n), n = 20L)
  p <- cbind(1, exp(0.3 + 0.5 * d$x), exp(-0.2 - 0.4 * d$x))
  p <- p / rowSums(p)
  d$Y <- t(vapply(seq_len(n), function(i) {
    as.vector(stats::rmultinom(1, 20, p[i, ]))
  }, numeric(3)))
  colnames(d$Y) <- c("a", "b", "c")
  d
}

correct_data_beta <- function() {
  set.seed(104)
  n <- 120
  d <- data.frame(x = rnorm(n), g = factor(rep(1:8, length.out = n)))
  mu <- stats::plogis(0.2 + 0.4 * d$x + rnorm(8, 0, 0.4)[d$g])
  phi <- exp(2 + rnorm(8, 0, 0.3)[d$g])
  d$yb <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  d
}

correct_data_mv <- function() {
  set.seed(105)
  n <- 80
  d <- data.frame(x = rnorm(n), g = factor(rep(1:8, length.out = n)))
  d$y1 <- 1 + d$x + rnorm(8)[d$g] + rnorm(n)
  d$y2 <- -1 + 0.5 * d$x + rnorm(8)[d$g] + rnorm(n)
  d
}

# brms's tests.priors.R "default_prior returns correct priors for
# multivariate models", seed 1
correct_data_p15 <- function() {
  set.seed(1)
  data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)), x = rnorm(10),
             g = rep(1:2, 5))
}

correct_data_offset <- function() {
  set.seed(106)
  n <- 50
  d <- data.frame(x = rnorm(n), z = rnorm(n), o = rnorm(n, 3, 0.5),
                  o2 = rnorm(n, 0.4, 0.1))
  d$y <- 2 + d$o + d$x + rnorm(n)
  d
}

correct_data_offset_pois <- function() {
  set.seed(107)
  n <- 50
  d <- data.frame(x = rnorm(n), expo = runif(n, 5, 50))
  d$cnt <- stats::rpois(n, d$expo * exp(-1 + 0.3 * d$x))
  d
}

correct_data_mi <- function() {
  set.seed(108)
  n <- 80
  d <- data.frame(z = rnorm(n))
  d$xm <- 0.5 * d$z + rnorm(n)
  d$y <- 1 + 0.7 * d$xm - 0.3 * d$z + rnorm(n)
  d$xm[sample.int(n, 12)] <- NA
  d
}

correct_data_mix <- function() {
  set.seed(109)
  n <- 400
  d <- data.frame(x = rnorm(n))
  # component 2 (mean 2) grows with x; component 1 (mean -1) is the rest
  p2 <- stats::plogis(-0.3 + 0.9 * d$x)
  k2 <- stats::rbinom(n, 1, p2) == 1
  d$y <- ifelse(k2, rnorm(n, 2, 0.7), rnorm(n, -1, 0.7))
  d
}

correct_data_counts <- function() {
  set.seed(110)
  n <- 60
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(rep(letters[1:6], length.out = n)))
  eta <- 0.3 + 0.5 * d$x + rnorm(6, 0, 0.5)[d$g]
  d$cnt <- stats::rpois(n, exp(eta))
  d$bin <- stats::rbinom(n, 1, stats::plogis(eta - 0.3))
  d
}
