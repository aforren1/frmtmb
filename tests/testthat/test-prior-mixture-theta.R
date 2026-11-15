# A mixture's weights are brms's class "theta": one simplex with a
# dirichlet() prior, listed by default_prior() and taken by set_prior(),
# where no weight has a formula or a fixed value (brms 2.23.0,
# dev/priors2-brms.R, cases M1 to M30). frmtmb's internal covariance
# class, "theta" through 0.82.0, is "thetacov".

pmt_data <- function(seed = 3) {
  set.seed(seed)
  data.frame(y = c(stats::rnorm(48, 0, 1), stats::rnorm(112, 4, 1)),
             x = stats::rnorm(160))
}

# the weights of a fit whose last component is the reference
pmt_weights <- function(fit) {
  b <- fit$estimates[["betad"]]
  v <- c(b[grepl("^theta[0-9]+_", names(b))], 0)
  exp(v) / sum(exp(v))
}

# an independent MAP of a two- or three-component gaussian mixture under
# dirichlet(a) on the softmax chart: the density of the log ratios is
# the Dirichlet times the Jacobian prod(w), so sum(a * log(w)) up to a
# constant
pmt_reference <- function(y, a, start) {
  K <- length(a)
  nlp <- function(p) {
    mu <- p[seq_len(K)]
    s <- exp(p[K + seq_len(K)])
    v <- c(p[2 * K + seq_len(K - 1L)], 0)
    lw <- v - log(sum(exp(v)))
    dens <- vapply(seq_len(K), function(k) {
      lw[k] + stats::dnorm(y, mu[k], s[k], log = TRUE)
    }, numeric(length(y)))
    m <- apply(dens, 1L, max)
    -sum(m + log(rowSums(exp(dens - m)))) - sum(a * lw)
  }
  o <- stats::optim(start, nlp, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 1000))
  v <- c(o$par[2 * K + seq_len(K - 1L)], 0)
  # `value` leaves out the Dirichlet's normalizing constant
  list(w = exp(v) / sum(exp(v)), value = o$value)
}

# a fit's internal estimates in pmt_reference()'s order: the component
# means, their log sds, the log ratios
pmt_start <- function(fit) {
  unname(c(fit$estimates[["beta"]], fit$estimates[["betad"]]))
}

test_that("default_prior() lists one theta row per weight simplex", {
  d <- pmt_data()
  m2 <- mixture(gaussian(), gaussian())
  tab <- default_prior(y ~ 1, data = d, family = m2)
  th <- tab[tab$class == "theta", ]
  expect_identical(nrow(th), 1L)
  expect_identical(th$prior, "(flat)")
  expect_true(all(c(th$coef, th$group, th$dpar, th$nlpar, th$resp) == ""))
  # brms's place in the table: after the class rows of no dpar, before
  # the components' rows (brms 2.23.0, M1)
  expect_identical(tab$class[seq_len(3L)], c("sigma1", "sigma2", "theta"))
  # three components, an ordinal mixture: still one row
  expect_identical(sum(default_prior(y ~ 1, data = d,
                                     family = mixture(gaussian(), gaussian(),
                                                      gaussian()))$class ==
                         "theta"), 1L)
  d$o <- cut(d$y, 3, labels = FALSE)
  expect_identical(sum(default_prior(o ~ x, data = d,
                                     family = mixture(cumulative(),
                                                      cumulative()))$class ==
                         "theta"), 1L)
  # one per response, under its resp (M6)
  d$y2 <- d$y + stats::rnorm(nrow(d))
  tmv <- default_prior(mvbind(y, y2) ~ 1, data = d, family = m2)
  expect_identical(tmv$resp[tmv$class == "theta"], c("y", "y2"))
  # the guard absent: a weight with a formula or a fixed value removes
  # the simplex (M3, M4), and a model with no mixture has none
  for (f in list(bf(y ~ 1, theta1 ~ x), bf(y ~ 1, theta1 = 0.4))) {
    expect_false(any(startsWith(default_prior(f, data = d,
                                              family = m2)$class, "theta")))
  }
  expect_false(any(default_prior(y ~ x, data = d)$class == "theta"))
})

test_that("class theta takes dirichlet() alone, as brms does", {
  d <- pmt_data()
  m2 <- mixture(gaussian(), gaussian())
  # brms's words (M11), and where frmtmb's old class went
  expect_error(set_prior("normal(0, 1)", class = "theta"),
               "Currently 'dirichlet' is the only valid prior for simplex",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(set_prior("normal(0, 1)", class = "theta"),
               "class \"thetacov\"", fixed = TRUE)
  # a narrowing names no parameter (M14 to M17, M30)
  for (a in list(list(coef = "a"), list(group = "g"),
                 list(dpar = "theta1"), list(nlpar = "a"))) {
    expect_error(do.call(set_prior, c(list("dirichlet(2)", class = "theta"),
                                      a)),
                 paste0("parameter: theta_", a[[1L]], " ~ dirichlet(2)"),
                 fixed = TRUE, class = "frmtmb_error")
  }
  expect_error(set_prior("dirichlet(2)", class = "theta", lb = 0),
               "takes no lb/ub", class = "frmtmb_error")
  # one weight is no class (M12, M13)
  expect_error(set_prior("dirichlet(2)", class = "theta1"),
               "one mixture proportion", class = "frmtmb_error")
  # the model decides the rest, in brms's words: a weight with a
  # formula (M19), no mixture (M20), the wrong length (M10), no resp in
  # a multivariate model (M22)
  expect_error(frm(bf(y ~ 1, theta1 ~ x), data = d, family = m2,
                   prior = set_prior("dirichlet(2)", class = "theta")),
               "do not correspond to any model parameter: theta ~",
               fixed = TRUE, class = "frmtmb_error")
  expect_error(frm(y ~ x, data = d,
                   prior = set_prior("dirichlet(2)", class = "theta")),
               "this model has no mixture()", fixed = TRUE)
  expect_error(frm(y ~ 1, data = d,
                   family = mixture(gaussian(), gaussian(), gaussian()),
                   prior = set_prior("dirichlet(c(1, 2))", class = "theta")),
               "Expected input of length 3", fixed = TRUE)
  d$y2 <- d$y + stats::rnorm(nrow(d))
  expect_error(validate_prior(set_prior("dirichlet(5)", class = "theta"),
                              mvbind(y, y2) ~ 1, data = d, family = m2),
               "resp = \"y2\"", fixed = TRUE)
  vp <- validate_prior(set_prior("dirichlet(5)", class = "theta",
                                 resp = "y2"),
                       mvbind(y, y2) ~ 1, data = d, family = m2)
  expect_identical(vp$prior[vp$class == "theta"], c("(flat)", "dirichlet(5)"))
})

test_that("frm() is the mode of the density frm_sample() samples", {
  d <- pmt_data()
  m2 <- mixture(gaussian(), gaussian())
  ml <- frm(y ~ 1, data = d, family = m2)
  # each fit is a mode of the independent density: started a step away
  # from it, the reference optimizer comes back to within a thousandth
  # of the move the prior made, and the objectives agree up to the
  # Dirichlet's constant. A mixture has several modes, so the reference
  # starts near the fit, and from the ML fit it finds none lower
  check <- function(fit, ml, y, a) {
    lab <- paste(a, collapse = ",")
    w <- pmt_weights(fit)
    move <- max(abs(w - pmt_weights(ml)))
    expect_gt(move, 0)
    near <- pmt_reference(y, a, pmt_start(fit) + 0.05)
    expect_lt(max(abs(w - near$w)) / move, 1e-3, label = lab)
    lognorm <- lgamma(sum(a)) - sum(lgamma(a))
    expect_equal(fit$opt$objective, near$value - lognorm, tolerance = 1e-8,
                 label = lab)
    far <- pmt_reference(y, a, pmt_start(ml))
    expect_lte(fit$opt$objective,
               far$value - lognorm + 1e-8 * abs(far$value))
    # the log prior is the Dirichlet of the weights times prod(w)
    lp <- lognorm + sum(a * log(w))
    expect_equal(fit$opt$objective,
                 -as.numeric(logLik(fit)) - lp, tolerance = 1e-8,
                 label = lab)
  }
  for (a in list(c(2, 2), c(20, 1), c(0.5, 0.5))) {
    pr <- set_prior(paste0("dirichlet(c(", a[1], ", ", a[2], "))"),
                    class = "theta")
    check(frm(y ~ 1, data = d, family = m2, prior = pr), ml, d$y, a)
  }
  # three components and unequal concentrations: alpha_k goes with
  # component k, and the last is the reference
  m3 <- mixture(gaussian(), gaussian(), gaussian())
  set.seed(8)
  d3 <- data.frame(y = c(stats::rnorm(50, 0), stats::rnorm(70, 4),
                         stats::rnorm(40, 8)))
  ml3 <- frm(y ~ 1, data = d3, family = m3)
  f3 <- frm(y ~ 1, data = d3, family = m3,
            prior = set_prior("dirichlet(c(9, 1, 3))", class = "theta"))
  check(f3, ml3, d3$y, c(9, 1, 3))
})

test_that("a theta row passed back is the same model", {
  d <- pmt_data()
  m2 <- mixture(gaussian(), gaussian())
  ml <- frm(y ~ 1, data = d, family = m2)
  # the default table's (flat) theta row applies nothing
  back <- frm(y ~ 1, data = d, family = m2, prior = default_prior(ml))
  expect_identical(back$opt$par, ml$opt$par)
  pr <- set_prior("dirichlet(c(3, 2))", class = "theta")
  fit <- frm(y ~ 1, data = d, family = m2, prior = pr)
  for (tab in list(validate_prior(pr, y ~ 1, data = d, family = m2),
                   prior_summary(fit))) {
    again <- frm(y ~ 1, data = d, family = m2, prior = tab)
    expect_identical(again$opt$par, fit$opt$par)
  }
  # brms's own spelling translates (brms 2.23.0, M24)
  skip_unless_brms()
  bf2 <- frm(y ~ 1, data = d, family = m2,
             prior = brms::prior(dirichlet(3, 2), class = theta))
  expect_identical(bf2$opt$par, fit$opt$par)
})

test_that("frm_simulate() refuses a draw from a dirichlet() prior", {
  d <- pmt_data()
  m2 <- mixture(gaussian(), gaussian())
  np <- list(b_mu1_Intercept = 0, b_mu2_Intercept = 4, sigma1 = 1,
             sigma2 = 1)
  expect_error(frm_simulate(y ~ 1, data = d, family = m2, nsim = 1,
                            newparams = np,
                            prior = set_prior("dirichlet(2)",
                                              class = "theta")),
               "is a density on a whole simplex", fixed = TRUE,
               class = "frmtmb_error")
})

test_that("the covariance class is thetacov, and class theta is not it", {
  set.seed(2)
  d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:8, 10)),
                  t = factor(rep(1:4, 20)))
  d$y <- stats::rnorm(80)
  f <- y ~ toep(0 + t | g)
  tab <- default_prior(f, data = d)
  expect_true(any(tab$class == "thetacov"))
  expect_false(any(tab$class == "theta"))
  cf <- tab$coef[tab$class == "thetacov" & nzchar(tab$coef)][1L]
  vp <- validate_prior(set_prior("normal(0, 1)", class = "thetacov",
                                 coef = cf), f, data = d)
  expect_identical(vp$prior[vp$class == "thetacov" & vp$coef == cf],
                   "normal(0, 1)")
  # the old spelling, with a density, a bound or neither, points at the
  # new class (review of lane priors2, L1)
  for (args in list(list("normal(0, 1)", coef = cf),
                    list("", coef = "theta_4", lb = -1, ub = 1),
                    list("", coef = "thetaac_1"))) {
    expect_error(do.call(set_prior, c(args, class = "theta")),
                 paste0("Write class = \"thetacov\", coef = \"",
                        args$coef, "\""), fixed = TRUE)
  }
  # an empty density on a narrowed weights row prints as brms's (flat)
  expect_error(set_prior("", class = "theta", coef = "a"),
               "theta_a ~ (flat).", fixed = TRUE)
})
