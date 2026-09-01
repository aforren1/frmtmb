# Regressions harvested from the OPEN issue trackers of brms, lme4 and
# glmmTMB. Each test names the upstream issue it came from: either a
# defect frmtmb shared and has now fixed, or a trap frmtmb dodges that
# is worth pinning down so a refactor cannot reintroduce it.

test_that("random-effect terms cannot be crossed with '*' or ':' (lme4#196)", {
  data(sleepstudy, package = "lme4")
  # reformulas hoists the bar out and silently fits the '+' model, so
  # lme4 accepts this typo and reports an interaction it never fit
  expect_error(
    frm(bf(Reaction ~ Days * (1 | Subject)) + gaussian(), data = sleepstudy),
    "cannot be crossed"
  )
  expect_error(
    frm(bf(Reaction ~ Days:(1 | Subject)) + gaussian(), data = sleepstudy),
    "cannot be crossed"
  )
  expect_error(
    frm(bf(Reaction ~ Days * (Days || Subject)) + gaussian(),
        data = sleepstudy),
    "cannot be crossed"
  )
  # the guard must not touch legitimate interactions or nested/crossed
  # grouping expressions, where ':' and '/' live inside the bar
  ss <- sleepstudy
  ss$b <- factor(ss$Days %% 3)
  ss$z <- as.numeric(ss$Days) / 10
  expect_s3_class(
    frm(bf(Reaction ~ Days:z + (1 | Subject)) + gaussian(), data = ss),
    "frmtmb_fit"
  )
  expect_s3_class(
    frm(bf(Reaction ~ Days + (1 | Subject:b)) + gaussian(), data = ss),
    "frmtmb_fit"
  )
  expect_s3_class(
    frm(bf(Reaction ~ Days + (1 | Subject/b)) + gaussian(), data = ss),
    "frmtmb_fit"
  )
})

test_that("mo()/mi() interaction multipliers must be numeric (brms#1828)", {
  set.seed(123)
  lev <- c("a1", "a2", "a3", "a4")
  income <- factor(sample(lev, 120, TRUE), levels = lev, ordered = TRUE)
  dat <- data.frame(
    income = income,
    ls = c(30, 60, 70, 75)[income] + rnorm(120, sd = 7),
    z = rnorm(120),
    xc = sample(c("a", "b"), 120, TRUE)
  )
  dat$xf <- factor(dat$xc)
  dat$xl <- dat$xc == "a"

  # a character multiplier used to pass the type gate (as.numeric() on a
  # character vector is still numeric, just all NA) and only surfaced as
  # "NA/NaN gradient evaluation" from the optimizer
  expect_error(frm(bf(ls ~ mo(income) * xc) + gaussian(), data = dat),
               "numeric multipliers only")
  expect_error(frm(bf(ls ~ mo(income) * xf) + gaussian(), data = dat),
               "numeric multipliers only")

  # numeric and logical multipliers keep working
  fit <- frm(bf(ls ~ mo(income) * z) + gaussian(), data = dat)
  expect_true("moincome:z" %in% names(fixef(fit)$mu))
  fit_l <- frm(bf(ls ~ mo(income) * xl) + gaussian(), data = dat)
  expect_true("moincome:xl" %in% names(fixef(fit_l)$mu))
})

test_that("anova() rejects fits with different numbers of observations (lme4#622)", {
  data(sleepstudy, package = "lme4")
  full <- frm(bf(Reaction ~ Days + (1 | Subject)) + gaussian(),
              data = sleepstudy)
  null <- frm(bf(Reaction ~ 1 + (1 | Subject)) + gaussian(),
              data = sleepstudy)
  short <- frm(bf(Reaction ~ 1 + (1 | Subject)) + gaussian(),
               data = sleepstudy[1:150, ])
  # without the guard this returns a negative Chisq and p = 1
  expect_error(anova(short, full), "same number of observations")
  tab <- anova(null, full)
  expect_true(tab[["Chisq"]][2] > 0)
})

test_that("zero prior weights equal subsetting (lme4#880)", {
  set.seed(1)
  gm <- rnorm(10, sd = 5) + c(rep(0, 8), -50, 50)
  dw <- data.frame(group = factor(rep(LETTERS[1:10], each = 3)),
                   x = rnorm(30))
  dw$y <- gm[as.integer(dw$group)] + 10 * dw$x + rnorm(30)
  dw$w <- as.integer(dw$group %in% LETTERS[1:8])

  sub <- frm(bf(y ~ x + (1 | group)) + gaussian(), data = subset(dw, w == 1))
  wt <- frm(bf(y | weights(w) ~ x + (1 | group)) + gaussian(), data = dw)
  # lme4 reports 4.24 vs 3.18 for these two; the weighted REML/ML
  # criterion there keeps a contribution from the zero-weight rows
  expect_equal(as.numeric(logLik(wt)), as.numeric(logLik(sub)),
               tolerance = 1e-6)
  expect_equal(unname(fixef(wt)$mu), unname(fixef(sub)$mu),
               tolerance = 1e-6)
  expect_equal(VarCorr(wt)[[1]][1, 1], VarCorr(sub)[[1]][1, 1],
               tolerance = 1e-6)
})

test_that("crossed grouping (1 | a * b) expands like a + b + a:b (lme4#234)", {
  set.seed(4)
  dd <- expand.grid(site = factor(1:8), year = factor(1:6), rep = 1:3)
  dd$y <- rnorm(nrow(dd))
  star <- frm(bf(y ~ 1 + (1 | site * year)) + gaussian(), data = dd)
  hand <- frm(bf(y ~ 1 + (1 | site) + (1 | year) + (1 | site:year)) +
                gaussian(), data = dd)
  expect_equal(as.numeric(logLik(star)), as.numeric(logLik(hand)),
               tolerance = 1e-6)
  expect_length(ranef(star), 3L)
})

test_that("grouping levels match across numeric and character newdata (lme4#616)", {
  set.seed(510)
  fd <- data.frame(state = sample(1:40, 200, replace = TRUE))
  fd$y <- rbinom(200, 1, fd$state / 40)
  fit <- frm(bf(y ~ (1 | state)) + bernoulli(), data = fd)
  # lme4 silently returns the wrong level's BLUP when newdata restates
  # an integer grouping variable as character
  expect_equal(predict(fit, newdata = fd),
               predict(fit, newdata = transform(fd,
                                                state = as.character(state))),
               tolerance = 1e-10)
})

test_that("REML predictions agree with fixef (glmmTMB#1143, #983)", {
  set.seed(7)
  td <- data.frame(cat = factor(rep(paste0("cat", 1:4), each = 50)),
                   year = factor(sample(1:4, 200, TRUE)))
  lk <- c(-0.01, 1.17, 0.95, 0.31)[as.integer(td$cat)] + rnorm(200)
  td$count <- rpois(200, exp(lk))
  for (reml in c(FALSE, TRUE)) {
    fit <- frm(bf(count ~ 0 + cat + (1 | year)) + poisson(), data = td,
               REML = reml)
    pr <- predict(fit, newdata = data.frame(cat = "cat1"), re.form = NA,
                  type = "link", allow_new_levels = TRUE)
    expect_equal(unname(pr), unname(fixef(fit)$mu[1]), tolerance = 1e-8)
  }
})

test_that("non-integer responses are rejected, not silently fit (lme4#682, #180)", {
  set.seed(3)
  bd <- data.frame(g = factor(rep(1:20, 5)))
  bd$Y <- ifelse(rbinom(100, 1, 0.4), 0.9, 0.1)
  # lme4 warns but fits, and shrinks the random-effect SD toward zero
  expect_error(frm(bf(Y ~ 1 + (1 | g)) + bernoulli(), data = bd),
               "must be 0/1")
  pd <- data.frame(g = bd$g, y = rpois(100, 3) + 0.5)
  expect_error(frm(bf(y ~ 1 + (1 | g)) + poisson(), data = pd),
               "non-negative integers")
})

test_that("discrete truncation normalizes with F(lb - 1) (brms#1903, #1923)", {
  set.seed(1)
  y <- rpois(6000, 3)
  y <- y[y >= 2 & y <= 6][1:400]
  fit <- frm(bf(y | trunc(lb = 2, ub = 6) ~ 1) + poisson(),
             data = data.frame(y = y))
  mu <- exp(fixef(fit)$mu[[1]])
  # the inclusive lower bound must keep its mass: P(2 <= Y <= 6) uses
  # ppois(1), not ppois(2)
  ll_ok <- sum(dpois(y, mu, log = TRUE)) -
    length(y) * log(ppois(6, mu) - ppois(1, mu))
  expect_equal(as.numeric(logLik(fit)), ll_ok, tolerance = 1e-6)
})

test_that("nonlinear fixed-effect SEs match nlme, not nlmer (lme4#819, #164)", {
  skip_if_not_installed("nlme")
  Soy <- nlme::Soybean
  fit <- frm(bf(weight ~ Asym / (1 + exp((xmid - Time) / scal)),
                Asym ~ 1 + (1 | Plot), xmid ~ 1, scal ~ 1, nl = TRUE) +
               gaussian(),
             data = Soy, start = list(beta = c(19, 55, 8)))
  se <- sqrt(diag(vcov(fit)))[1:3]
  # nlmer reports about 0.0063/0.0065/0.0065 here - two orders of
  # magnitude too small; nlme gets 0.629/0.458/0.290
  expect_equal(unname(se), c(0.629, 0.458, 0.290), tolerance = 0.05)
  # nlmer's predict() errors as soon as newdata is supplied (lme4#164)
  expect_equal(predict(fit, newdata = Soy), predict(fit), tolerance = 1e-8)
})

test_that("tensor-product smooths are supported or refused clearly (glmmTMB#1082)", {
  set.seed(101)
  dd <- data.frame(x = runif(200), z = runif(200))
  dd$y <- sin(2 * dd$x) + dd$z^2 + rnorm(200, sd = 0.2)
  # te()/ti() have no random-effect representation; say so rather than
  # dropping the term or erroring inside the model frame
  expect_error(frm(bf(y ~ te(x, z)) + gaussian(), data = dd),
               "not supported")
  fit <- frm(bf(y ~ t2(x, z)) + gaussian(), data = dd)
  ref <- mgcv::gam(y ~ t2(x, z), data = dd, method = "ML")
  expect_equal(as.numeric(fitted(fit)), as.numeric(fitted(ref)),
               tolerance = 1e-4)
})
