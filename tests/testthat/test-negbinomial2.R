# brmsfamily("negbinomial2"): brms 2.23.0's negative binomial in
# sigma = 1 / shape. The family is negbinomial()'s under fam_reparam(),
# so most of what is asserted here is an identity with negbinomial(),
# checked bit for bit where the two run the same arithmetic. The
# references that are not negbinomial() itself are stats::dnbinom() and
# the analytic mean and variance. brms's own program is the gated
# test-nb2dm-brms.R.

nb2_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                   g = gl(15, n / 15))
  dd$e <- exp(stats::runif(n))
  dd$y <- stats::rnbinom(n, mu = exp(1 + 0.4 * dd$x), size = 2.5)
  dd
}

# a relative difference in units of the larger value's rounding, so a
# bound reads as a count of ulps rather than as a tolerance
nb2_ulps <- function(a, b) {
  abs(a - b) / (.Machine$double.eps * pmax(abs(a), abs(b)))
}

test_that("brmsfamily('negbinomial2') is brms's family", {
  fam <- brmsfamily("negbinomial2")
  expect_s3_class(fam, "frmtmb_family")
  expect_identical(fam$family, "negbinomial2")
  expect_identical(fam$dpars, c("mu", "sigma"))
  expect_identical(fam$link, "log")
  expect_identical(fam$link_sigma, "log")
  expect_null(fam$link_shape)
  expect_identical(brmsfamily("negbinomial2", "sqrt",
                              link_sigma = "softplus")$link_sigma,
                   "softplus")
  expect_identical(brmsfamily("negbinomial2", "sqrt")$link, "sqrt")
  # brms's sentence, and brms's set: brms:::.family_negbinomial2()$links
  expect_error(brmsfamily("negbinomial2", "logit"),
               "'logit' is not a supported link for family 'negbinomial2'",
               fixed = TRUE)
  expect_error(brmsfamily("negbinomial2", link_sigma = "logit"),
               "not a supported link for parameter 'sigma'", fixed = TRUE)
  # brms does not export a constructor, and neither does frmtmb
  expect_false("negbinomial2" %in% getNamespaceExports("frmtmb"))
})

test_that("a brms negbinomial2 family object is read with its links", {
  skip_if_not_installed("brms")
  f <- as_frmtmb_family(brms::brmsfamily("negbinomial2",
                                         link_sigma = "identity"))
  expect_identical(f$family, "negbinomial2")
  expect_identical(f$link_sigma, "identity")
})

test_that("the density is negbinomial's at shape = 1 / sigma", {
  nb <- negbinomial()
  nb2 <- brmsfamily("negbinomial2")
  g <- expand.grid(y = c(0, 1, 3, 12, 80), eta = c(-3, 0, 2.5),
                   sigma = c(1e-6, 0.05, 0.7, 4, 300))
  dbl2 <- row_lpdf(nb2, g$y, g$y, list(mu = exp(g$eta), sigma = g$sigma),
                   list(), NULL)
  dbl1 <- row_lpdf(nb, g$y, g$y, list(mu = exp(g$eta), shape = 1 / g$sigma),
                   list(), NULL)
  expect_identical(dbl2, dbl1)
  # on the tape, where the objective hands each dpar's linear predictor
  # along, against negbinomial() at the mirrored shape, 1 / sigma
  t2 <- RTMB::MakeTape(function(p) {
    n <- length(g$y)
    row_lpdf(nb2, g$y, g$y,
             list(mu = exp(p[1:n]), sigma = exp(p[n + 1:n]),
                  .eta_mu = p[1:n], .eta_sigma = p[n + 1:n]), list(), NULL)
  }, c(g$eta, log(g$sigma)))
  t1 <- RTMB::MakeTape(function(p) {
    n <- length(g$y)
    row_lpdf(nb, g$y, g$y,
             list(mu = exp(p[1:n]), shape = 1 / exp(p[n + 1:n]),
                  .eta_mu = p[1:n]), list(), NULL)
  }, c(g$eta, log(g$sigma)))
  p0 <- c(g$eta, log(g$sigma))
  expect_identical(t2(p0), t1(p0))
  expect_identical(t2$jacobian(p0), t1$jacobian(p0))
  # and against stats::dnbinom where that is accurate (sigma >= 1e-4:
  # below it dnbinom loses digits to its Poisson branch,
  # dev/nb2dm-findings.md)
  ok <- g$sigma >= 1e-4
  ref <- stats::dnbinom(g$y, size = 1 / g$sigma, mu = exp(g$eta), log = TRUE)
  expect_lt(max(nb2_ulps(as.numeric(t2(p0))[ok], ref[ok])), 4096)
})

test_that("a response that is not a count is refused in negbinomial2's name", {
  dd <- nb2_data(30)
  dd$w <- dd$y + 0.5
  expect_error(frm(w ~ x, family = brmsfamily("negbinomial2"), data = dd),
               "negbinomial2: response must be non-negative integers",
               fixed = TRUE)
})

test_that("a fit is negbinomial's fit with sigma = 1 / shape", {
  dd <- nb2_data(31)
  f1 <- frm(y ~ x, family = negbinomial(), data = dd)
  f2 <- frm(y ~ x, family = brmsfamily("negbinomial2"), data = dd)
  l1 <- as.numeric(logLik(f1))
  l2 <- as.numeric(logLik(f2))
  # two optimizer runs on one function: equal to the optimizer's stop
  expect_lt(abs(l1 - l2), sqrt(.Machine$double.eps) * abs(l1))
  b1 <- fixef_by_dpar(f1)
  b2 <- fixef_by_dpar(f2)
  expect_lt(max(abs(b2$mu - b1$mu) / abs(b1$mu)), 1e-5)
  expect_lt(abs(b2$sigma + b1$shape) / abs(b1$shape), 1e-5)
  s1 <- sqrt(diag(vcov(f1)))
  s2 <- sqrt(diag(vcov(f2)))
  expect_lt(max(abs(s2 / s1 - 1)), 1e-4)
  # sigma(), as brms reports it
  expect_lt(abs(sigma(f2) * exp(b1$shape[[1]]) - 1), 1e-5)
  # sigma ~ z is shape ~ -z
  f3 <- frm(bf(y ~ x, sigma ~ z), family = brmsfamily("negbinomial2"),
            data = dd)
  f4 <- frm(bf(y ~ x, shape ~ z), family = negbinomial(), data = dd)
  expect_lt(abs(as.numeric(logLik(f3)) - as.numeric(logLik(f4))),
            sqrt(.Machine$double.eps) * abs(l1))
  expect_lt(max(abs(fixef_by_dpar(f3)$sigma + fixef_by_dpar(f4)$shape) /
                  abs(fixef_by_dpar(f4)$shape)), 1e-4)
})

test_that("rate() and a group term carry over from negbinomial()", {
  dd <- nb2_data(32)
  f1 <- frm(y | rate(e) ~ x, family = negbinomial(), data = dd)
  f2 <- frm(y | rate(e) ~ x, family = brmsfamily("negbinomial2"), data = dd)
  l1 <- as.numeric(logLik(f1))
  expect_lt(abs(as.numeric(logLik(f2)) - l1),
            sqrt(.Machine$double.eps) * abs(l1))
  # brms multiplies the mean AND the shape 1 / sigma by the exposure
  fm <- fitted(f2)[, "Estimate"]
  expect_lt(max(abs(fm / fitted(f1)[, "Estimate"] - 1)), 1e-5)
  r1 <- frm(y ~ x + (1 | g), family = negbinomial(), data = dd)
  r2 <- frm(y ~ x + (1 | g), family = brmsfamily("negbinomial2"), data = dd)
  lr <- as.numeric(logLik(r1))
  expect_lt(abs(as.numeric(logLik(r2)) - lr), sqrt(.Machine$double.eps) *
              abs(lr))
})

test_that("the post-fit paths read negbinomial's mean and variance", {
  dd <- nb2_data(33)
  f1 <- frm(y ~ x, family = negbinomial(), data = dd)
  f2 <- frm(y ~ x, family = brmsfamily("negbinomial2"), data = dd)
  # the estimates agree to the optimizer's stop, about 1e-8 relative, so
  # every derived quantity agrees to that and no closer
  near <- function(a, b) max(abs(a - b)) <= 1e-5 * max(abs(b))
  expect_true(near(fitted(f2), fitted(f1)))
  expect_true(near(residuals(f2, type = "pearson")[, "Estimate"],
                   residuals(f1, type = "pearson")[, "Estimate"]))
  expect_true(near(frm_linpred(f2, dpar = "sigma"),
                   -frm_linpred(f1, dpar = "shape")))
  nd <- data.frame(x = c(-1, 0, 2))
  expect_true(near(fitted(f2, newdata = nd), fitted(f1, newdata = nd)))
  s <- simulate(f2, nsim = 3, seed = 4)
  expect_identical(dim(s), c(nrow(dd), 3L))
  expect_true(all(s >= 0 & s == round(s)))
  # the deviance is negbinomial's: its residuals are defined, not refused
  expect_true(near(residuals(f2, type = "deviance")[, "Estimate"],
                   residuals(f1, type = "deviance")[, "Estimate"]))
})

test_that("the simulator's moments are mu and mu + sigma mu^2", {
  nb2 <- brmsfamily("negbinomial2")
  set.seed(77)
  n <- 2e5
  z <- stats::qnorm(1 - 1e-6 / 2)
  for (pt in list(c(mu = 3, sigma = 0.5), c(mu = 0.4, sigma = 2.5),
                  c(mu = 25, sigma = 0.02))) {
    x <- nb2$sim(list(mu = rep(pt[["mu"]], n), sigma = rep(pt[["sigma"]], n)),
                 list(), n)
    m <- pt[["mu"]]
    v <- m + pt[["sigma"]] * m^2
    lab <- paste0("mu ", m, ", sigma ", pt[["sigma"]])
    expect_lt(abs(mean(x) - m), z * sqrt(v / n), label = paste(lab, "mean"))
    # the standard error of a sample variance, from the sample itself
    se_v <- stats::sd((x - mean(x))^2) / sqrt(n)
    expect_lt(abs(stats::var(x) - v), z * se_v, label = paste(lab, "var"))
    # and the family's own variance function says the same
    expect_identical(nb2$post$var_fn(list(mu = m, sigma = pt[["sigma"]]),
                                     list()),
                     m + m^2 / (1 / pt[["sigma"]]))
  }
})

test_that("every negbinomial slot reaches negbinomial2, a later one too", {
  nb <- negbinomial()
  nb2 <- brmsfamily("negbinomial2")
  takes_dpars <- function(f) {
    is.function(f) && "dpars" %in% names(formals(f))
  }
  for (nm in names(nb)) {
    v <- nb[[nm]]
    if (takes_dpars(v)) {
      expect_true(takes_dpars(nb2[[nm]]), label = nm)
      expect_identical(body(attr(nb2[[nm]], "frmtmb_inner")), body(v),
                       label = nm)
    }
  }
  for (k in names(nb$post)) {
    expect_identical(body(attr(nb2$post[[k]], "frmtmb_inner")),
                     body(nb$post[[k]]),
                     label = k)
  }
  expect_identical(nb2$accepts_aterms, nb$accepts_aterms)
  expect_identical(nb2$type, nb$type)
  # a slot negbinomial() gains (a CDF, say) is negbinomial2's at once,
  # read at shape = 1 / sigma, with the term it accepts
  real <- frmtmb:::fam_negbinomial
  local_mocked_bindings(fam_negbinomial = function(link = "log",
                                                   link_shape = "log") {
    f <- real(link, link_shape)
    # the slot contract: probabilities, P(Y <= q)
    f[["lcdf"]] <- function(q, dpars, aterms) {
      stats::pnbinom(q, size = dpars[["shape"]], mu = dpars[["mu"]])
    }
    f[["accepts_aterms"]] <- c(f[["accepts_aterms"]], "cens")
    f
  }, .package = "frmtmb")
  later <- brmsfamily("negbinomial2")
  expect_true("cens" %in% later$accepts_aterms)
  expect_identical(later$lcdf(4, list(mu = 2, sigma = 0.25), list()),
                   stats::pnbinom(4, size = 4, mu = 2))
  expect_length(formals(later$lcdf), 3L)
})

test_that("it starts where negbinomial() starts", {
  dd <- nb2_data(34)
  f1 <- frm(y ~ x, family = negbinomial(), data = dd)
  f2 <- frm(y ~ x, family = brmsfamily("negbinomial2"), data = dd)
  p1 <- frmtmb:::start_ls_plan(f1$frame)
  p2 <- frmtmb:::start_ls_plan(f2$frame)
  expect_length(p2, 1L)
  expect_identical(p2[[1]]$beta, p1[[1]]$beta)
  expect_identical(p2[[1]]$scale$sigma, 1 / p1[[1]]$scale$shape)
})

test_that("the compatibility table answers for negbinomial2", {
  tab <- frm_compat("negbinomial2")
  row <- tab[tab$feature_a == "rate()" | tab$feature_b == "rate()", ]
  expect_identical(row$status, "works")
  # what negbinomial() refuses, this refuses with it
  cen <- tab[tab$feature_a == "cens()" | tab$feature_b == "cens()", ]
  nb <- frm_compat("negbinomial")
  cen1 <- nb[nb$feature_a == "cens()" | nb$feature_b == "cens()", ]
  expect_identical(cen$status, cen1$status)
  expect_true("negbinomial2" %in%
                frmtmb:::frmtmb_compat_groups_lst$discrete)
})
