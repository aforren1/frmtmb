# Equating one distributional parameter to another, brms's
# bf(y ~ x, sigma1 = "sigma2"). The rules and their messages are
# brms 2.23.0's (brms:::brmsformula(), read in dev/formula2-brms-src/);
# the names are those of a brms fit of the same model
# (dev/formula2-brms-equate-fit.R and its log).

equate_data <- function() {
  set.seed(11)
  n <- 300
  d <- data.frame(x = rnorm(n))
  k <- rbinom(n, 1, 0.4)
  d$y <- ifelse(k == 1, rnorm(n, 3 + 0.5 * d$x, 1), rnorm(n, -1, 1))
  d
}

test_that("bf() records an equation and applies brms's rules to it", {
  f <- bf(y ~ x, sigma1 = "sigma2")
  expect_identical(f$pfix, list(sigma1 = "sigma2"))
  # the five rules of brms's bf(), with brms's messages
  expect_error(bf(y ~ x, sigma1 = "sigmaa2"),
               "Can only equate parameters of the same class",
               class = "frmtmb_error")
  expect_error(bf(y ~ x, mu3 = "mu2"),
               "Equating parameters of class 'mu' is not allowed",
               class = "frmtmb_error")
  expect_error(bf(y ~ x, sigma1 = "sigma1"),
               "Equating 'sigma1' with itself is not meaningful",
               class = "frmtmb_error")
  expect_error(bf(y ~ x, shape1 ~ x, shape2 = "shape1"),
               "Cannot use predicted parameters on the right-hand side",
               class = "frmtmb_error")
  expect_error(bf(y ~ x, shape1 = "shape3", shape2 = "shape1"),
               "Cannot use fixed parameters on the right-hand side",
               class = "frmtmb_error")
  # brms stops in an internal sum on this one; frmtmb says why
  expect_error(bf(y ~ x, theta1 = "theta2"),
               "Equating mixing proportions", class = "frmtmb_error")
  # the rules hold on every route an equation or a formula arrives by
  expect_error(bf(bf(y ~ x, sigma1 = "sigma2"), sigma2 ~ x),
               "Cannot use predicted parameters", class = "frmtmb_error")
  expect_error(bf(y ~ x, sigma1 = "sigma2") + lf(sigma2 ~ x),
               "Cannot use predicted parameters", class = "frmtmb_error")
  expect_error(lf(sigma2 ~ x, sigma1 = "sigma2"),
               "Cannot use predicted parameters", class = "frmtmb_error")
  # a string that names nothing is still refused, not read as a formula
  expect_error(bf(y ~ x, "sigma2"), "Cannot interpret bf[(][)] argument",
               class = "frmtmb_error")
  expect_error(bf(y ~ x, sigma1 = c("sigma2", "sigma3")),
               "Cannot interpret bf[(][)] argument 'sigma1'",
               class = "frmtmb_error")
})

test_that("an equation naming no parameter of the family is refused", {
  d <- equate_data()
  fam <- mixture(gaussian(), gaussian())
  expect_error(frm(bf(y ~ x, sigma1 = "sigma3"), family = fam, data = d),
               "Parameter 'sigma3' cannot be found", class = "frmtmb_error")
  expect_error(frm(bf(y ~ x, sigma3 = "sigma1"), family = fam, data = d),
               "not available for family", class = "frmtmb_error")
})

test_that("an equated mixture is the model with one shared sigma", {
  d <- equate_data()
  fam <- mixture(gaussian(), gaussian())
  fit <- frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)

  # the same model written by hand: one log sigma for both components,
  # theta2 the reference as in frmtmb
  X <- cbind(1, d$x)
  y <- d$y
  nll <- function(p) {
    mu1 <- as.vector(X %*% p$b1)
    mu2 <- as.vector(X %*% p$b2)
    s <- exp(p$ls)
    a <- p$t1 - log(1 + exp(p$t1)) + RTMB::dnorm(y, mu1, s, log = TRUE)
    b <- -log(1 + exp(p$t1)) + RTMB::dnorm(y, mu2, s, log = TRUE)
    m <- 0.5 * (a + b + abs(a - b))
    -sum(m + log(exp(a - m) + exp(b - m)))
  }
  est <- fit$estimates
  p_frm <- list(b1 = unname(est$beta[1:2]), b2 = unname(est$beta[3:4]),
                ls = unname(est$betad[1]), t1 = unname(est$betad[2]))
  obj <- RTMB::MakeADFun(nll, p_frm, silent = TRUE)
  ll <- as.numeric(logLik(fit))
  # one parameter fewer than the mixture with two sigmas
  expect_identical(attr(logLik(fit), "df"), 6L)
  expect_length(est$betad, 2L)
  # the same function at the same point
  expect_lt(abs(obj$fn(obj$par) + ll), 1e-12 * abs(ll))
  # and the same optimum, from the hand-written model's own start
  st <- obj$par
  st[] <- 0
  st[c(1, 3)] <- c(-1, 3)
  opt <- stats::nlminb(st, obj$fn, obj$gr)
  expect_lt(abs(opt$objective + ll), 1e-9 * abs(ll))
  se <- sqrt(diag(solve(obj$he(opt$par))))
  expect_lt(max(abs(opt$par - obj$par) / se), 1e-3)

  # the lf() spelling is the same model
  fit_lf <- frm(bf(y ~ x) + lf(sigma1 = "sigma2"), family = fam, data = d)
  expect_identical(as.numeric(logLik(fit_lf)), ll)
  # sigma1 reads sigma2's value at every row
  expect_identical(frm_linpred(fit, dpar = "sigma1"),
                   frm_linpred(fit, dpar = "sigma2"))
})

test_that("an equated dpar is named as brms names it", {
  d <- equate_data()
  fam <- mixture(gaussian(), gaussian())
  fit <- frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)
  # brms 2.23.0's variables() of this model, less Intercept_mu1,
  # Intercept_mu2, lprior and lp__, which an ML fit does not have
  # (dev/formula2-brms-equate-fit.log)
  brms_vars <- c("b_mu1_Intercept", "b_mu2_Intercept", "b_mu1_x",
                 "b_mu2_x", "sigma1", "sigma2", "theta1", "theta2")
  expect_setequal(variables(fit), brms_vars)
  h <- hypothesis(fit, "sigma1 = sigma2", class = NULL)
  expect_identical(h$hypothesis$Estimate, 0)
  # brms's summary lists both under Further Distributional Parameters
  sp <- summary(fit)$spec_pars
  expect_identical(rownames(sp)[1:2], c("sigma1", "sigma2"))
  expect_identical(sp["sigma1", ], sp["sigma2", ],
                   ignore_attr = TRUE)
  # sigma2 is the one parameter: the prior table and fixef() hold it
  # alone, and a prior on sigma1 is refused as addressing nothing
  pr <- get_prior(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)
  expect_true("sigma2" %in% pr$class)
  expect_false("sigma1" %in% pr$class)
  expect_error(frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d,
                   prior = set_prior("normal(1, 1)", class = "sigma1")),
               "sigma1 is equated to sigma2", class = "frmtmb_error")
  expect_named(par_template(fit)$betad,
               c("sigma2_(Intercept)", "theta1_(Intercept)"))
})

test_that("a prediction reads the equated value", {
  d <- equate_data()
  fam <- mixture(gaussian(), gaussian())
  fit <- frm(bf(y ~ x, sigma1 = "sigma2"), family = fam, data = d)
  nd <- d[1:5, ]
  expect_equal(fitted(fit, newdata = nd), fitted(fit)[1:5, ],
               ignore_attr = TRUE)
  sim <- simulate(fit, nsim = 2, seed = 1)
  expect_identical(dim(as.matrix(sim)), c(nrow(d), 2L))
})

test_that("coef() has no block for an equated dpar", {
  d <- equate_data()
  fit <- frm(bf(y ~ x, sigma1 = "sigma2"),
             family = mixture(gaussian(), gaussian()), data = d)
  cf <- coef(fit)
  expect_false("sigma1" %in% names(cf))
  expect_true("sigma2" %in% names(cf))
})

test_that("an equation onto a link-scale dpar lists no natural value", {
  # an hmm() transition cell (frmtmb.latent) is a link-scale dpar: brms
  # has no such parameter, so it stays a coefficient, and the cell
  # equated to it is no parameter at all. The family's declaration is
  # set by hand here, which is the one field that rule reads
  # (dev/formula2-rev-hmm.R has the hmm() case itself).
  d <- equate_data()
  fit <- frm(bf(y ~ x, sigma1 = "sigma2"),
             family = mixture(gaussian(), gaussian()), data = d)
  fit$spec$responses$y$family$link_scale_dpars <- c("sigma1", "sigma2")
  v <- variables(fit)
  expect_false("sigma1" %in% v)
  expect_true("b_sigma2_Intercept" %in% v)
  expect_false("sigma1" %in% rownames(summary(fit)$spec_pars))
})
