# One Student-t latent scale per grouping factor, as brms 2.23.0 builds
# it: every Student-t term of a factor reads the same dfm_1, in mu and
# in sigma alike, and a term on another factor gets its own
# (dev/prilabel-log/brms-ref.txt, brms-ref2.txt; dev/prilabel-brms-ref.R).

st_data <- function(n = 360, seed = 21) {
  set.seed(seed)
  g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
  h <- factor(sample(sprintf("h%d", 1:8), n, TRUE))
  x <- rnorm(n)
  z <- rnorm(n)
  w <- sqrt(5 / rchisq(24, 5))
  u0 <- (w * rnorm(24, 0, 0.8))[as.integer(g)]
  u1 <- (w * rnorm(24, 0, 0.4))[as.integer(g)]
  y <- 1 + 0.5 * x + 0.3 * z + u0 + u1 * x + rnorm(8, 0, 0.5)[as.integer(h)] +
    rnorm(n, 0, 0.6)
  data.frame(y, x, z, g, h)
}

test_that("two Student-t terms of one factor are one multivariate t", {
  dd <- st_data()
  two <- frm(bf(y ~ x + (1 | gr(g, dist = "student")) +
                  (0 + x | gr(g, dist = "student"))),
             family = gaussian(), data = dd)
  one <- frm(bf(y ~ x + (1 + x || gr(g, dist = "student"))),
             family = gaussian(), data = dd)
  # the same model, the shared scale written as one diag() block: the
  # same maximum and the same scales
  expect_equal(as.numeric(logLik(two)), as.numeric(logLik(one)))
  expect_equal(unname(two$estimates$theta), unname(one$estimates$theta),
               tolerance = 1e-5)
  expect_identical(two$frame$re_blocks[[1L]]$t_share, 1:2)
})

test_that("the shared density is the mixture over one scale per level", {
  dd <- st_data()
  fit <- frm(bf(y ~ x + (1 + x | gr(g, dist = "student", dist_nu = 4)) +
                  (0 + z | gr(g, dist = "student", dist_nu = 4))),
             family = gaussian(), data = dd)
  bks <- fit$frame$re_blocks
  expect_identical(bks[[2L]]$t_share, 1:2)
  b <- fit$estimates$b
  th <- fit$estimates$theta
  S1 <- frmtmb:::covstruct_registry$us_t$vcov(th[bks[[1L]]$theta_idx],
                                             bks[[1L]])
  S2 <- frmtmb:::covstruct_registry$us_t$vcov(th[bks[[2L]]$theta_idx],
                                             bks[[2L]])
  S <- matrix(0, 3, 3)
  S[1:2, 1:2] <- S1
  S[3, 3] <- S2
  B <- rbind(matrix(b[bks[[1L]]$b_idx], 2), b[bks[[2L]]$b_idx])
  # the independent reference: per level, the gaussian with scale S / l
  # averaged over l ~ gamma(nu / 2, rate nu / 2) by numerical integration
  nu <- 4
  Si <- solve(S)
  ldS <- as.numeric(determinant(S)$modulus)
  ref <- sum(apply(B, 2L, function(bj) {
    q <- drop(crossprod(bj, Si %*% bj))
    f <- function(l) {
      exp(1.5 * log(l) - 1.5 * log(2 * pi) - ldS / 2 - l * q / 2 +
            stats::dgamma(l, nu / 2, nu / 2, log = TRUE))
    }
    log(stats::integrate(f, 0, Inf, rel.tol = 1e-12)$value)
  }))
  got <- frmtmb:::student_shared_lpdf(bks, b, th)
  expect_equal(got, ref, tolerance = 1e-9)
  # and the fit's tape is that density plus the gaussian likelihood:
  # obj$env$f() at the joint mode is the joint negative log density
  lp <- Filter(function(l) identical(l$dpar, "mu"),
               fit$frame$linpreds)[[1L]]
  beta <- fit$estimates$beta
  eta <- as.vector(lp$X %*% beta[lp$idx]) + as.vector(lp$Z %*% b)
  sig <- exp(fit$estimates$betad[[1L]])
  joint <- -sum(stats::dnorm(dd$y, eta, sig, log = TRUE)) - ref
  par <- fit$obj$env$last.par.best
  expect_equal(as.numeric(fit$obj$env$f(par)), joint, tolerance = 1e-9)
})

test_that("a term in sigma shares the scale of its factor, as in brms", {
  dd <- st_data()
  fit <- frm(bf(y ~ x + (1 | gr(g, dist = "student")),
                sigma ~ (1 | gr(g, dist = "student"))),
             family = gaussian(), data = dd)
  expect_identical(lapply(fit$frame$re_blocks, `[[`, "t_share"),
                   list(1:2, 1:2))
})

test_that("only the Student-t terms of one factor share a scale", {
  dd <- st_data()
  a <- frm(bf(y ~ x + (1 | gr(g, dist = "student")) +
                (1 | gr(h, dist = "student")) + (0 + x | g)),
           family = gaussian(), data = dd, dry_run = "frame")
  expect_true(all(vapply(a$re_blocks, function(bk) is.null(bk$t_share),
                         NA)))
  # brms has one df per factor: two dist_nu on one factor are refused
  expect_error(frm(bf(y ~ x + (1 | gr(g, dist = "student", dist_nu = 3)) +
                        (0 + x | gr(g, dist = "student", dist_nu = 6))),
                   family = gaussian(), data = dd),
               "take one dist_nu. Got 3 and 6", fixed = TRUE)
})

test_that("draws of shared terms read one mixing variable per level", {
  dd <- st_data()
  fit <- frm(bf(y ~ x + (1 | gr(g, dist = "student")) +
                  (0 + x | gr(g, dist = "student"))),
             family = gaussian(), data = dd)
  bks <- fit$frame$re_blocks
  n <- bks[[1L]]$n_levels
  sc <- exp(fit$estimates$theta)
  set.seed(9)
  bd <- frmtmb:::draw_b(fit)
  # the first block draws its normals and the level's mixing variable,
  # the second its normals and the same mixing variable
  set.seed(9)
  z1 <- stats::rnorm(n)
  w <- sqrt(5 / stats::rchisq(n, df = 5))
  z2 <- stats::rnorm(n)
  expect_equal(bd[bks[[1L]]$b_idx], z1 * sc[1L] * w)
  expect_equal(bd[bks[[2L]]$b_idx], z2 * sc[2L] * w)
  # so a simulation that drops one and keeps the other is refused, and
  # one that drops both runs
  expect_error(simulate(fit, nsim = 1, seed = 1,
                        re_formula = ~ (1 | gr(g, dist = "student"))),
               "keeps some of the Student-t terms", fixed = TRUE)
  expect_length(simulate(fit, nsim = 2, seed = 1, re_formula = NA), 2L)
  # quadrature integrates one effect at a time
  expect_error(frm(bf(y ~ x + (1 | gr(g, dist = "student")) +
                        (0 + x | gr(g, dist = "student"))),
                   family = gaussian(), data = dd, quadrature = TRUE),
               "share one latent scale", fixed = TRUE)
})

# lane prilabel's review, finding 1 (dev/prilabel-rev-tby.R, data seed
# 21): frmtmb splits gr(g, by = f) into one block per by-level, each
# over its own levels of g, and the first sharing refused that.
st_by_data <- function() {
  set.seed(21)
  n <- 360
  g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
  d <- data.frame(g, x = rnorm(n), z = rnorm(n),
                  hg = factor(ifelse(as.integer(g) <= 12L, "u", "v")))
  d$y <- 1 + 0.5 * d$x + rnorm(24, 0, 0.8)[d$g] + rnorm(n, 0, 0.6)
  d
}

test_that("a Student-t by-term alone is the model of 0.71.0", {
  d <- st_by_data()
  # rellib-r9 (0.71.0), dev/prilabel-rev-log/tby-base.txt: one block per
  # by-level, each its own multivariate t; they hold different levels
  # of g, so there is nothing to share
  ll <- c(-361.6239067, -364.5278394)
  forms <- list(y ~ x + (1 | gr(g, by = hg, dist = "student")),
                y ~ x + (1 + x | gr(g, by = hg, dist = "student")))
  for (k in 1:2) {
    # nlminb can stop with "singular convergence (7)" here; frmtmb's
    # Newton check finds the stop at the maximum and says so
    fit <- allow_warnings(frm(bf(forms[[k]]), family = gaussian(),
                              data = d),
                          "Optimizer did not report convergence")
    expect_equal(as.numeric(logLik(fit)), ll[k], tolerance = 1e-8)
    expect_true(all(vapply(fit$frame$re_blocks,
                           function(bk) is.null(bk$t_share), NA)))
  }
  # and a scalar by-term still takes quadrature, as on 0.71.0
  q <- frm(bf(forms[[1]]), family = gaussian(), data = d,
           quadrature = TRUE)
  expect_true(is.finite(as.numeric(logLik(q))))
})

test_that("a by-term and a plain term share one scale per level of g", {
  d <- st_by_data()
  # brms 2.23.0 writes dfm_1 into both terms, over every level of g
  # (dev/prilabel-rev-log/misc.txt)
  fit <- frm(bf(y ~ x + (1 | gr(g, by = hg, dist = "student")) +
                  (0 + z | gr(g, dist = "student"))),
             family = gaussian(), data = d)
  bks <- fit$frame$re_blocks
  expect_identical(lapply(bks, `[[`, "t_share"), rep(list(1:3), 3))
  # the independent reference: per level of g, a bivariate t with nu 5
  # over (the by-block's intercept, the plain block's slope), scale
  # diag(s_by^2, s_z^2), written from the textbook density
  nu <- 5
  lp <- Filter(function(l) identical(l$dpar, "mu"),
               fit$frame$linpreds)[[1L]]
  ref_nll <- function(par) {
    pv <- split(unname(par), names(par))
    b <- pv$b
    th <- pv$theta
    eta <- as.vector(lp$X %*% pv$beta[lp$idx]) + as.vector(lp$Z %*% b)
    ll <- sum(stats::dnorm(d$y, eta, exp(pv$betad[1L]), log = TRUE))
    lev <- bks[[3L]]$levels
    for (j in seq_along(lev)) {
      k <- if (lev[j] %in% bks[[1L]]$levels) 1L else 2L
      at <- match(lev[j], bks[[k]]$levels)
      v <- c(b[bks[[k]]$b_idx][at], b[bks[[3L]]$b_idx][j])
      s <- exp(c(th[bks[[k]]$theta_idx], th[bks[[3L]]$theta_idx]))
      q <- sum((v / s)^2)
      ll <- ll + lgamma((nu + 2) / 2) - lgamma(nu / 2) - log(nu * pi) -
        sum(log(s)) - (nu + 2) / 2 * log(1 + q / nu)
    }
    -ll
  }
  par0 <- fit$obj$env$last.par.best
  set.seed(7)
  for (r in 0:3) {
    par <- par0
    if (r) par <- par + stats::rnorm(length(par), 0, 0.3)
    expect_equal(as.numeric(fit$obj$env$f(par)), ref_nll(par),
                 tolerance = 1e-10)
  }
  # draws: one mixing variable per level of g for all three blocks
  n1 <- bks[[1L]]$n_levels
  n2 <- bks[[2L]]$n_levels
  sc <- exp(fit$estimates$theta)
  set.seed(9)
  bd <- frmtmb:::draw_b(fit)
  set.seed(9)
  z1 <- stats::rnorm(n1)
  w <- sqrt(nu / stats::rchisq(n1 + n2, df = nu))
  z2 <- stats::rnorm(n2)
  z3 <- stats::rnorm(n1 + n2)
  expect_equal(bd[bks[[1L]]$b_idx], z1 * sc[1L] * w[seq_len(n1)])
  expect_equal(bd[bks[[2L]]$b_idx], z2 * sc[2L] * w[n1 + seq_len(n2)])
  expect_equal(bd[bks[[3L]]$b_idx], z3 * sc[3L] * w)
  # dropping the plain term keeps the by-term at levels they share
  expect_error(simulate(fit, nsim = 1, seed = 1,
                        re_formula = ~ (1 | gr(g, by = hg))),
               "keeps some of the Student-t terms", fixed = TRUE)
})
