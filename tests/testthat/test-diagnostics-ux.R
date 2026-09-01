# The diagnostics / UX cluster of dev/test-backlog.md: the papercuts a
# user porting from glm, lme4, glmmTMB or brms meets before anything
# statistical goes wrong. Each test names the backlog entry it closes.

# --- cbind(successes, failures) responses [glmmTMB#1319, #1325] -------

make_binom <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 5), each = 5)))
  d$trials <- rep(10L, n)
  d$succ <- rbinom(n, 10, plogis(0.3 + 0.5 * d$x))
  d$fail <- d$trials - d$succ
  d
}

test_that("cbind(successes, failures) matches glm() exactly (glmmTMB#1319)", {
  d <- make_binom()
  f <- frm(cbind(succ, fail) ~ x, data = d, family = binomial())
  g <- glm(cbind(succ, fail) ~ x, data = d, family = binomial())
  expect_equal(unname(fixef(f)$mu), unname(coef(g)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
               tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(f)))[1:2]),
               unname(summary(g)$coefficients[, 2]), tolerance = 1e-5)
})

test_that("cbind() is exactly the trials() spelling (glmmTMB#1325)", {
  d <- make_binom()
  f1 <- frm(cbind(succ, fail) ~ x + (1 | g), data = d, family = binomial())
  f2 <- frm(succ | trials(trials) ~ x + (1 | g), data = d,
            family = binomial())
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
  expect_equal(fixef(f1)$mu, fixef(f2)$mu)
  # the response really is the success count, and the aterm the total
  expect_equal(f1$frame$y[[1]], d$succ)
  expect_equal(f1$frame$aterm_values[[1]]$trials, as.numeric(d$trials))
  expect_equal(fitted(f1), fitted(f2))
})

test_that("cbind() responses are refused where they are ambiguous", {
  d <- make_binom()
  expect_error(
    frm(cbind(succ, fail) | trials(trials) ~ x, data = d,
        family = binomial()),
    "already carries the number of trials"
  )
  expect_error(
    frm(cbind(succ, fail, x) ~ x, data = d, family = binomial()),
    "exactly two unnamed columns"
  )
  d2 <- d
  d2$fail <- d2$fail + 0.5
  expect_error(
    frm(cbind(succ, fail) ~ x, data = d2, family = binomial()),
    "non-negative integer counts"
  )
})

# --- degenerate fits: no free outer parameters [glmmTMB#1325, #1317] --

test_that("a model with no free parameters fits degenerately (glmmTMB#1317)", {
  d <- make_binom()
  f <- frm(succ | trials(trials) ~ 0, data = d, family = binomial())
  expect_length(f$opt$par, 0L)
  expect_equal(f$opt$convergence, 0L)
  expect_match(f$opt$message, "no free parameters")
  # the likelihood is the template's, evaluated once: mu = 0.5 on the
  # logit scale with no coefficients to move it
  expect_equal(as.numeric(logLik(f)),
               sum(dbinom(d$succ, d$trials, 0.5, log = TRUE)))
  expect_equal(attr(logLik(f), "df"), 0L)
  # the post-fit surface still answers rather than dying in nlminb
  expect_equal(unname(fitted(f)), rep(5, nrow(d)))
  expect_s3_class(summary(f), "summary.frmtmb_fit")
  expect_type(diagnose(f, quiet = TRUE), "list")
})

# --- grouping-factor structure checks [lme4 lmerControl] --------------

make_gauss <- function(n = 120, seed = 4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 6), each = 6)))
  d$y <- 1 + d$x + rnorm(n)
  d$one <- factor(rep("a", n))
  d$obs <- factor(seq_len(n))
  d$cnt <- rpois(n, 3)
  d
}

test_that("a one-level grouping factor is reported (lme4 lmerControl)", {
  d <- make_gauss()
  expect_warning(frm(y ~ x + (1 | one), data = d, family = gaussian()),
                 "single level")
  expect_error(
    frm(y ~ x + (1 | one), data = d, family = gaussian(),
        control = frmtmb_control(check_nlev_1 = "stop")),
    "single level"
  )
  expect_silent(
    f <- frm(y ~ x + (1 | one), data = d, family = gaussian(),
             control = frmtmb_control(check_nlev_1 = "ignore"))
  )
  # ignoring it still fits: the variance simply collapses to zero
  expect_lt(as.data.frame(VarCorr(f))$sdcor[1], 1e-3)
})

test_that("gaussian OLRE warns about confounding with sigma (lme4)", {
  d <- make_gauss()
  expect_warning(frm(y ~ x + (1 | obs), data = d, family = gaussian()),
                 "confounded with the residual sd")
  expect_silent(
    frm(y ~ x + (1 | obs), data = d, family = gaussian(),
        control = frmtmb_control(check_olre = "ignore"))
  )
  # overdispersion for a discrete family is the legitimate use
  expect_silent(frm(cnt ~ x + (1 | obs), data = d, family = poisson()))
  # se() pins the residual sd row by row, so the split is identified:
  # this is the random-effects meta-analysis
  dm <- data.frame(yi = rnorm(12), sei = runif(12, 0.2, 0.6),
                   obs = factor(1:12))
  expect_silent(
    frm(bf(yi | se(sei) ~ 1 + (1 | obs)) + gaussian(), data = dm)
  )
  # and an ordinary grouping factor is never flagged
  expect_silent(frm(y ~ x + (1 | g), data = d, family = gaussian()))
})

# --- diagnose() upgrades [glmmTMB diagnose(), lme4 isSingular] --------

test_that("diagnose() names complete separation (glmmTMB diagnose())", {
  set.seed(7)
  ds <- data.frame(z = rep(0:1, each = 25))
  ds$xx <- c(rnorm(25, -3), rnorm(25, 3))
  f <- suppressWarnings(frm(z ~ xx, data = ds, family = bernoulli()))
  dg <- diagnose(f, quiet = TRUE)
  expect_s3_class(dg$separation, "data.frame")
  expect_true(any(grepl("xx", dg$separation$parameter)))
  expect_gt(abs(dg$separation$estimate[1]), 10)
  expect_output(diagnose(f), "complete separation")
  # a well-behaved binomial fit is not flagged
  d <- make_binom()
  expect_null(diagnose(frm(cbind(succ, fail) ~ x, data = d,
                           family = binomial()),
                       quiet = TRUE)$separation)
})

test_that("diagnose() points badly scaled predictors at autoscale", {
  set.seed(8)
  dp <- data.frame(y = rnorm(60))
  dp$xbig <- rnorm(60) * 1e5
  dg <- diagnose(frm(y ~ xbig, data = dp, family = gaussian()),
                 quiet = TRUE)
  expect_s3_class(dg$predictor_scale, "data.frame")
  expect_true(any(grepl("xbig", dg$predictor_scale$column)))
  expect_output(diagnose(frm(y ~ xbig, data = dp, family = gaussian())),
                "autoscale")
  dp$xok <- rnorm(60)
  expect_null(diagnose(frm(y ~ xok, data = dp, family = gaussian()),
                       quiet = TRUE)$predictor_scale)
})

test_that("diagnose() gives an isSingular verdict (lme4 isSingular)", {
  d <- make_gauss()
  f <- frm(y ~ x + (1 | one), data = d, family = gaussian(),
           control = frmtmb_control(check_nlev_1 = "ignore"))
  dg <- diagnose(f, quiet = TRUE)
  # the verdict is read off the estimates, so it stands even though the
  # Hessian is positive definite and the gradient is tiny
  expect_true(dg$pdHess)
  expect_s3_class(dg$singular, "data.frame")
  expect_lt(dg$singular$value[1], 1e-4)
  expect_output(diagnose(f), "Singular fit")
})

test_that("diagnose() works on a fit with no random effects", {
  # theta is NULL there, and abs(NULL) is an error, not an empty result
  d <- make_gauss()
  expect_type(diagnose(frm(y ~ x, data = d, family = gaussian()),
                       quiet = TRUE), "list")
  expect_output(diagnose(frm(y ~ x, data = d, family = gaussian())),
                "No convergence problems detected")
})

# --- simulate() return types [glmmTMB test-simulate.R; lme4#737] ------

make_ordinal <- function(n = 300, seed = 5) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  eta <- 1.2 * d$x
  P <- cbind(plogis(-0.5 - eta),
             plogis(1.0 - eta) - plogis(-0.5 - eta),
             1 - plogis(1.0 - eta))
  k <- 1L + rowSums(t(apply(P, 1, cumsum)) < runif(n))
  d$o <- factor(c("lo", "mid", "hi")[k], levels = c("lo", "mid", "hi"),
                ordered = TRUE)
  d
}

test_that("ordinal simulate() returns ordered factors (lme4#737)", {
  d <- make_ordinal()
  for (fam in list(cumulative(), sratio(), cratio(), acat())) {
    f <- frm(o ~ x, data = d, family = fam)
    s <- simulate(f, nsim = 2)
    expect_true(is.ordered(s[[1]]), info = fam$family)
    expect_equal(levels(s[[1]]), levels(d$o), info = fam$family)
    expect_equal(nrow(s), nrow(d), info = fam$family)
  }
})

test_that("ordinal category probabilities reproduce the lpdf exactly", {
  # the simulator needs the whole category distribution; the lpdf only
  # scores the observed one, so pin them against each other
  d <- make_ordinal()
  n <- nrow(d)
  for (fam in c("cumulative", "sratio", "cratio", "acat")) {
    f <- frm(o ~ x, data = d, family = get(fam, envir = asNamespace(
      "frmtmb"))())
    ordered <- fam %in% c("cumulative", "sratio")
    tau <- frmtmb:::ord_tau_from_raw(f$estimates$tau_raw, ordered)
    eta <- as.numeric(model.matrix(~ x - 1, d) %*% fixef(f)$mu)
    P <- frmtmb:::ord_cat_probs(fam, eta, tau, NULL, "logit")
    expect_equal(unname(rowSums(P)), rep(1, n), tolerance = 1e-12,
                 info = fam)
    yv <- f$frame$y[[1]]
    lp <- f$spec$responses[[1]]$family$lpdf(
      yv, list(mu = eta), list(), list(tau_raw = f$estimates$tau_raw))
    expect_equal(log(P[cbind(seq_len(n), yv)]), as.numeric(lp),
                 tolerance = 1e-10, info = fam)
  }
})

test_that("multinomial simulate() returns matrices (glmmTMB test-simulate.R)", {
  set.seed(6)
  m <- data.frame(x = rnorm(150))
  m$tr <- rep(6L, 150)
  mm <- t(sapply(seq_len(150), function(i) rmultinom(1, 6, c(.3, .4, .3))))
  colnames(mm) <- c("c1", "c2", "c3")
  m$Y <- mm
  f <- frm(Y | trials(tr) ~ x, data = m, family = multinomial(K = 3))
  s <- simulate(f, nsim = 2)
  expect_true(is.matrix(s[[1]]))
  expect_equal(dim(s[[1]]), c(150L, 3L))
  expect_equal(colnames(s[[1]]), c("c1", "c2", "c3"))
  expect_equal(unique(rowSums(s[[1]])), 6)
})

test_that("simulate() respects na.exclude padding (glmmTMB test-simulate.R)", {
  d <- make_ordinal(n = 200)
  d$x[c(3, 9)] <- NA
  f <- frm(o ~ x, data = d, family = cumulative(), na.action = na.exclude)
  s <- simulate(f, nsim = 2)
  expect_equal(nrow(s), nrow(d))
  expect_true(all(is.na(s[[1]][c(3, 9)])))
  expect_false(anyNA(s[[1]][-c(3, 9)]))
  # na.omit keeps the old fitted-row shape
  f2 <- frm(o ~ x, data = d, family = cumulative())
  expect_equal(nrow(simulate(f2, nsim = 1)), nrow(d) - 2L)
  # the internal consumers work in fitted-row space and must unpad
  d$yg <- rnorm(nrow(d))
  fg <- frm(yg ~ x, data = d, family = gaussian(), na.action = na.exclude)
  expect_equal(nrow(simulate(fg, nsim = 1)), nrow(d))
  # 3 draws x (mu intercept, mu slope, sigma intercept)
  expect_equal(dim(frm_bootstrap(fg, nsim = 3, seed = 1)$t), c(3L, 3L))
})

# --- nl parameter names and starting values [brms#391, #734] ----------

test_that("an nlpar colliding with a data column is refused (brms#391)", {
  set.seed(9)
  d <- data.frame(x = rnorm(80), a = rnorm(80), b = rnorm(80))
  d$yn <- 2 * exp(0.5 * d$x) + rnorm(80, sd = 0.1)
  # the nonlinear body resolved 'a' and 'b' to the PARAMETERS and
  # silently ignored the identically named columns
  expect_error(
    frm(bf(yn ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d,
        family = gaussian()),
    "also name columns of the data"
  )
  # renaming them fits
  f <- frm(bf(yn ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE), data = d,
           family = gaussian())
  expect_equal(unname(fixef(f)$A), 2, tolerance = 0.1)
  expect_equal(unname(fixef(f)$B), 0.5, tolerance = 0.1)
})

test_that("a failed nl fit names start= (brms#734)", {
  set.seed(10)
  dd <- data.frame(t = seq(0.5, 5, length.out = 40))
  dd$y <- 3 * exp(dd$t / 2) + rnorm(40, sd = 0.1)
  # B = 0 is the default start and t / 0 is not finite there
  expect_error(
    frm(bf(y ~ A * exp(t / B), A ~ 1, B ~ 1, nl = TRUE), data = dd,
        family = gaussian()),
    "start"
  )
  # with a start in the right region the same model fits
  f <- frm(bf(y ~ A * exp(t / B), A ~ 1, B ~ 1, nl = TRUE), data = dd,
           family = gaussian(), start = list(beta = c(3, 2)))
  expect_equal(unname(fixef(f)$A), 3, tolerance = 0.2)
  expect_equal(unname(fixef(f)$B), 2, tolerance = 0.2)
})

# --- REML anova() [glmmTMB#776] ---------------------------------------

test_that("REML anova() allows equal fixed-effect designs (glmmTMB#776)", {
  d <- make_gauss()
  d$z <- rnorm(nrow(d))
  f1 <- frm(y ~ x + (1 | g), data = d, family = gaussian(), REML = TRUE)
  f2 <- frm(y ~ x, data = d, family = gaussian(), REML = TRUE)
  # same X, different random-effect structure: the REML likelihoods are
  # for the same error contrasts, so the LRT is meaningful
  tab <- anova(f1, f2)
  expect_s3_class(tab, "anova")
  expect_equal(nrow(tab), 2L)
  expect_true(is.finite(tab$Chisq[2]))
  # a reordered but equivalent design spans the same column space
  g1 <- frm(y ~ x + z + (1 | g), data = d, family = gaussian(),
            REML = TRUE)
  g2 <- frm(y ~ z + x + (1 | g), data = d, family = gaussian(),
            REML = TRUE)
  expect_s3_class(anova(g1, g2), "anova")
  expect_equal(as.numeric(logLik(g1)), as.numeric(logLik(g2)),
               tolerance = 1e-6)
})

test_that("REML anova() still refuses different fixed effects (glmmTMB#776)", {
  d <- make_gauss()
  f1 <- frm(y ~ x + (1 | g), data = d, family = gaussian(), REML = TRUE)
  f2 <- frm(y ~ 1 + (1 | g), data = d, family = gaussian(), REML = TRUE)
  expect_error(anova(f1, f2), "same column space")
  # and refuses to mix the two methods
  f3 <- frm(y ~ x + (1 | g), data = d, family = gaussian())
  expect_error(anova(f1, f3), "cannot mix REML and ML")
  # ML comparisons are untouched
  f4 <- frm(y ~ 1 + (1 | g), data = d, family = gaussian())
  expect_s3_class(anova(f3, f4), "anova")
})

# --- offset() in dpar formulas [glmmTMB test-offset.R] -----------------

test_that("offset() in a dpar formula reaches the likelihood", {
  # this already worked; the test pins it down so it cannot regress into
  # a silent drop, which is the failure mode glmmTMB#625 describes
  set.seed(11)
  n <- 120
  d <- data.frame(x = rnorm(n), off = runif(n, 0.5, 1.5))
  d$y <- rnorm(n, 1 + d$x, sd = exp(0.2 + 0.3 * d$x))
  f <- frm(bf(y ~ x, sigma ~ x + offset(off)), data = d,
           family = gaussian())
  mu <- as.numeric(model.matrix(~x, d) %*% fixef(f)$mu)
  sg <- exp(as.numeric(model.matrix(~x, d) %*% fixef(f)$sigma) + d$off)
  expect_equal(as.numeric(logLik(f)),
               sum(dnorm(d$y, mu, sg, log = TRUE)), tolerance = 1e-8)
  # the offset moves the estimate: dropping it is not a no-op
  f0 <- frm(bf(y ~ x, sigma ~ x), data = d, family = gaussian())
  expect_gt(abs(fixef(f)$sigma[1] - fixef(f0)$sigma[1]), 0.5)
  # and it follows through to prediction, in sample and on newdata
  expect_equal(as.numeric(predict(f, dpar = "sigma", type = "link")),
               log(sg))
  nd <- data.frame(x = c(-1, 0, 1), off = c(0, 1, 2))
  expect_equal(
    as.numeric(predict(f, newdata = nd, dpar = "sigma", type = "link")),
    as.numeric(model.matrix(~x, nd) %*% fixef(f)$sigma) + nd$off
  )
})
