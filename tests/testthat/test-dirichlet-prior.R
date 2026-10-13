# brms's simo prior, dirichlet(alpha) on a mo() simplex
# (dev/moparity-findings.md, item 2). On 0.69.0 (rellib-r7) frm()
# refused every class = "simo" row ("brms's \"simo\" is the Dirichlet on
# a mo() simplex. frmtmb holds that simplex as its free coordinates and
# puts no density on it at all") and set_prior() had no Dirichlet:
# "Unsupported prior distribution 'dirichlet'".

dir_data <- function() {
  set.seed(3)
  n <- 120
  dm <- data.frame(inc = sample(0:3, n, TRUE), z = stats::rnorm(n))
  dm$y <- 1 + 0.6 * c(0, 1, 1.4, 1.5)[dm$inc + 1] + 0.3 * dm$z +
    stats::rnorm(n)
  dm
}

test_that("set_prior() reads brms's Dirichlet and its rules", {
  p <- set_prior("dirichlet(c(5, 1, 1))", class = "simo", coef = "mox1")
  expect_identical(unclass(p)[[1L]]$dist$alpha, c(5, 1, 1))
  expect_identical(unclass(p)[[1L]]$prior, "dirichlet(c(5, 1, 1))")
  p1 <- set_prior("dirichlet(2)", class = "simo", coef = "mox1")
  expect_identical(unclass(p1)[[1L]]$dist$alpha, 2)
  # brms's three refusals (dev/moparity-brms-simo.R)
  expect_error(set_prior("normal(0, 1)", class = "simo", coef = "mox1"),
               "only valid prior for simplex")
  expect_error(set_prior("dirichlet(2)", class = "simo"),
               "addresses one simplex")
  expect_error(set_prior("dirichlet(2)", class = "simo", coef = "mox1",
                         lb = 0), "no lb/ub")
  expect_error(set_prior("dirichlet(2)", class = "b"), "class = \"simo\"")
  expect_error(set_prior("dirichlet(c(1, x))", class = "simo",
                         coef = "mox1"), "Cannot parse")
  expect_error(prior_dirichlet(c(1, -1)), "positive")
})

test_that("a simo row resolves onto its simplex, by brms's coef", {
  dm <- dir_data()
  pr <- set_prior("dirichlet(c(2, 1.5, 3))", class = "simo",
                  coef = "moinc:z1")
  fit <- frm(bf(y ~ mo(inc) * z), family = gaussian(), data = dm,
             prior = pr)
  ent <- frmtmb:::resolve_prior_input(fit, fit$prior)$entries
  expect_length(ent, 1L)
  # brms numbers mo(inc) before mo(inc):z, and the interaction is zeta2
  expect_identical(ent[[1L]]$comp, "zeta2")
  expect_identical(ent[[1L]]$dist$alpha, c(2, 1.5, 3))
  expect_error(frm(bf(y ~ mo(inc) * z), family = gaussian(), data = dm,
                   prior = set_prior("dirichlet(c(1, 1))", class = "simo",
                                     coef = "moinc1")),
               "Expected input of length 3")
  expect_error(frm(bf(y ~ mo(inc) * z), family = gaussian(), data = dm,
                   prior = set_prior("dirichlet(1)", class = "simo",
                                     coef = "mox1")),
               "moinc1")
})

test_that("a MAP fit under a simo prior is the hand log posterior's mode", {
  dm <- dir_data()
  a <- c(2, 1.5, 3)
  pr <- set_prior("dirichlet(c(2, 1.5, 3))", class = "simo",
                  coef = "moinc1")
  fit <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm,
             prior = pr)
  est <- fit$estimates
  hand <- function(p) {
    beta <- p[1:3]
    zeta <- p[5:6]
    w <- frmtmb:::mo_simplex(zeta)
    mo <- 3 * c(0, cumsum(w))[dm$inc + 1]
    mu <- beta[1] + beta[2] * dm$z + beta[3] * mo
    ll <- sum(stats::dnorm(dm$y, mu, exp(p[4]), log = TRUE))
    # the Dirichlet on w, times the stereographic chart's Jacobian: the
    # sphere point's coordinates, the conformal factor, and the 2^D
    # sheets the chart covers
    r2 <- sum(zeta^2)
    u <- (2 * as.vector(frmtmb:::mo_sphere_basis(3) %*% zeta) -
            (r2 - 1) / sqrt(3)) / (r2 + 1)
    ll + lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(w)) +
      sum(log(abs(u))) + 2 * log(4) - 2 * log(1 + r2) - 3 * log(2)
  }
  p_fit <- c(est$beta, est$betad, est$zeta1)
  h <- hand(p_fit)
  expect_lt(abs(-fit$obj$fn(fit$opt$par) - h) / abs(h), 1e-12)
  # and it is the mode: the hand function from another start comes back
  op <- stats::optim(p_fit + c(0.2, -0.1, 0.2, 0.1, 0.3, -0.2),
                     function(p) -hand(p), method = "BFGS",
                     control = list(reltol = 1e-15, maxit = 5000))
  expect_lt(-op$value - h, 1e-9 * abs(h))
  # the inverse: the ML fit's weights are elsewhere, far beyond how
  # closely two searches for the MAP weights agree
  ml <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm)
  w_map <- frmtmb:::mo_simplex(est$zeta1)
  w_op <- frmtmb:::mo_simplex(op$par[5:6])
  w_ml <- frmtmb:::mo_simplex(ml$estimates$zeta1)
  expect_gt(max(abs(w_map - w_ml)), 1e3 * max(abs(w_map - w_op)))
})

test_that("a concentration below 1/2 has no MAP and is refused", {
  dm <- dir_data()
  expect_error(frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm,
                   prior = set_prior("dirichlet(0.4)", class = "simo",
                                     coef = "moinc1")),
               "has no maximum")
  # the guard absent: 1/2 is the smallest concentration with a mode
  fit <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm,
             prior = set_prior("dirichlet(0.5)", class = "simo",
                               coef = "moinc1"))
  expect_true(all(is.finite(fit$opt$par)))
})

test_that("the chart Jacobians match finite differences", {
  set.seed(2)
  fd <- function(z, chart) {
    f <- function(z) frmtmb:::mo_simplex(z, chart)[-1L]
    h <- 1e-6
    J <- vapply(seq_along(z), function(j) {
      e <- numeric(length(z))
      e[j] <- h
      (f(z + e) - f(z - e)) / (2 * h)
    }, numeric(length(z)))
    log(abs(det(J)))
  }
  gaps <- c()
  rel <- c()
  for (D in 3:4) {
    a <- stats::runif(D, 0.5, 3)
    for (chart in list(NULL, "softmax")) {
      for (i in 1:5) {
        z <- stats::rnorm(D - 1)
        w <- frmtmb:::mo_simplex(z, chart)
        terms <- c(lgamma(sum(a)), -lgamma(a), (a - 1) * log(w),
                   fd(z, chart), if (is.null(chart)) -D * log(2))
        want <- sum(terms)
        got <- frmtmb:::dirichlet_logdens(z,
                                          frmtmb:::dirichlet_dist(a, D, chart))
        gaps <- c(gaps, abs(got - want))
        # relative to the size of the terms summed, which is what the
        # rounding of the sum and of the differences scales with; the
        # sum itself can be near 0
        rel <- c(rel, abs(got - want) / sum(abs(terms)))
      }
    }
  }
  # central differences at h = 1e-6 are good to about 1e-9 here; the
  # inverse below shows the check would see a missing Jacobian
  expect_lt(max(rel), 1e-7)
  z <- c(0.3, -0.4)
  no_jac <- lgamma(4) - 3 * lgamma(1.5) +
    sum(0.5 * log(frmtmb:::mo_simplex(z)))
  expect_gt(abs(no_jac - frmtmb:::dirichlet_logdens(
    z, frmtmb:::dirichlet_dist(1.5, 3L))), 1e4 * max(gaps))
})

test_that("get_prior() lists a simo row per simplex, as brms does", {
  dm <- dir_data()
  tab <- get_prior(bf(y ~ mo(inc) * z, sigma ~ mo(inc)), data = dm,
                   family = gaussian())
  s <- tab[tab$class == "simo", ]
  expect_setequal(paste0(s$coef, "|", s$dpar),
                  c("moinc1|", "moinc:z1|", "moinc1|sigma"))
  # the fit route applies no density there
  expect_true(all(s$prior == "(flat)"))
})

test_that("brms's own default_prior() table translates", {
  skip_unless_brms()
  dm <- dir_data()
  bp <- suppressMessages(brms::default_prior(y ~ mo(inc), data = dm))
  bp <- bp[bp$class == "simo", ]
  bp$prior <- "dirichlet(c(3, 1, 1))"
  fit <- frm(bf(y ~ mo(inc)), family = gaussian(), data = dm, prior = bp)
  ent <- frmtmb:::resolve_prior_input(fit, fit$prior)$entries
  expect_identical(ent[[1L]]$dist$alpha, c(3, 1, 1))
})
