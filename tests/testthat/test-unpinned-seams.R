# Seams the simulate re-check left unpinned: a refusal and two internal
# call sites that no test exercised. Each block says what would have to
# change for it to fail, because a pin whose failure mode is not
# obvious is a pin nobody trusts.

test_that("method = 'predict' refuses a family that draws the response whole", {
  skip_on_cran()

  # mixture_mvn() carries a sim_ctx and no sim: its draw is an n by D
  # matrix for the WHOLE response, which simulate() consumes and the
  # per-row quantiles of a prediction band cannot. So simulate() works
  # on this fit and conditional_effects(method = "predict") refuses,
  # and the refusal has to arrive as a message naming the reason rather
  # than as a length error from inside the loop - which is what it was
  # until the family checks moved out of the grid loop.
  set.seed(42)
  n <- 200
  cl <- rbinom(n, 1, 0.4)
  x <- rnorm(n)
  Y <- cbind(ifelse(cl == 1, 0, 3) + 0.5 * x + rnorm(n),
             ifelse(cl == 1, 0, 4) + rnorm(n))
  dd <- data.frame(x = x)
  dd$Y <- Y
  fit <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd)

  # simulate() accepts it, through the sim_ctx route
  sm <- simulate(fit, nsim = 2, seed = 3)
  expect_identical(dim(sm[[1]]), dim(Y))

  # conditional_effects() does not, by name
  expect_error(conditional_effects(fit, method = "predict"),
               "draws the response whole")
  # and the refusal is raised before any grid work, so a family with no
  # simulator at all reaches its own message rather than this one
  expect_error(conditional_effects(fit, method = "predict",
                                   effects = "nosuchvariable"),
               "draws the response whole")

  # the expected-response display is refused too, and for the other
  # reason: this family's mean is an n by D matrix, so there is no one
  # curve to draw or to put a delta-method band on
  expect_error(conditional_effects(fit),
               "not one number per observation")
  # what IS drawable is a named parameter, and the mixing weight comes
  # back on the softmax scale like any other mixture's
  ce <- conditional_effects(fit, dpar = "theta1")
  expect_identical(names(ce), "x")
  expect_true(all(ce$x$estimate__ > 0 & ce$x$estimate__ < 1))
})

test_that("a finalized family reaches the spec par_template and get_prior read", {
  skip_on_cran()

  # R/par-template.R and R/priors.R both call
  # carry_finalized_responses() before building their table, so that a
  # family which rewrote itself against the response is priced under the
  # vocabulary it ENDED with rather than the one it was written with.
  #
  # Neither call site is observable from outside today, and this block
  # is the measurement rather than a claim: with a family whose
  # family_finalize() replaces a link, and again with one that REORDERS
  # family$dpars, get_prior() and par_template() return identical
  # output whether or not the two lines are there (checked against a
  # scratch install with both removed, on both routes of get_prior()
  # and on the two par_template() paths that reach
  # resolve_prior_input()). A finalize cannot swap which dpar is
  # PRIMARY - the response spec's own primary_dpars is never refreshed,
  # which is half the reason the seam is unobservable; the other half
  # is that neither function reads a dpar's link. So what is pinned
  # here is the helper's contract, which is what those call sites
  # depend on.
  fam <- frmtmb_family(
    "bounded_sigma",
    dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      stats::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
    },
    family_finalize = function(fam, y, aterms) {
      ub <- 4 * stats::sd(y)
      fam$links$sigma <- list(
        name = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
        linkfun = function(mu) log(mu / (ub - mu)),
        linkinv = function(eta) ub / (1 + exp(-eta)),
        mu_eta = function(eta) {
          p <- 1 / (1 + exp(-eta))
          ub * p * (1 - p)
        }
      )
      fam
    }
  )

  set.seed(1)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- 2 + dd$x + stats::rnorm(80)

  spec <- parse_spec(as_bform(bf(y ~ x, sigma ~ x) + fam, NULL))
  frame <- assemble_frame(spec, dd)
  # written with a log link
  expect_identical(spec$responses$y$family$links$sigma$name, "log")
  # and finalized to the scaled logit the data decided
  expect_match(frame$spec$responses$y$family$links$sigma$name,
               "^scaled_logit")

  carried <- carry_finalized_responses(spec, frame)
  expect_identical(carried$responses$y$family$links$sigma$name,
                   frame$spec$responses$y$family$links$sigma$name)
  expect_identical(carried$responses$y$dpars[[2]]$link$name,
                   frame$spec$responses$y$family$links$sigma$name)
  # a family with no family_finalize() is passed through untouched, so
  # the helper cannot cost anything on the models that do not need it
  spec2 <- parse_spec(as_bform(bf(y ~ x) + gaussian(), NULL))
  frame2 <- assemble_frame(spec2, dd)
  expect_identical(carry_finalized_responses(spec2, frame2), spec2)

  # and both call sites still answer on such a family, which is the
  # regression the lines were written against
  gp <- get_prior(bf(y ~ x, sigma ~ x) + fam, data = dd)
  expect_true(all(c("Intercept", "b") %in% gp$class))
  expect_true("sigma" %in% gp$dpar)
  pt <- par_template(bf(y ~ x, sigma ~ x) + fam, data = dd)
  expect_named(unclass(pt), c("beta", "betad"), ignore.order = TRUE)
})

test_that("a reported dpar scale survives the rows na.action dropped", {
  skip_on_cran()

  # A mixture's mixing weight reports on the softmax scale, which is a
  # function of ALL the dpars and so is read through predict() rather
  # than from this predictor's own eta. Those reads go through the same
  # na.action bookkeeping as any other prediction, and the value has to
  # come back on the same rows as the standard error it is paired with:
  # a padding applied twice, or once on one of the pair only, is a
  # length mismatch rather than a wrong number.
  set.seed(4)
  n <- 300
  dx <- data.frame(x = stats::rnorm(n))
  p1 <- stats::plogis(0.8 + 0.5 * dx$x)
  dx$y <- ifelse(stats::rbinom(n, 1, p1), stats::rnorm(n, -1, 0.6),
                 stats::rnorm(n, 2, 0.9))
  dx$x[c(3, 17, 100)] <- NA          # three rows na.omit drops
  fit <- suppressMessages(
    frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()),
        data = dx))
  expect_identical(fit$frame$n_obs, 297L)

  rv <- as.numeric(predict(fit, type = "response", dpar = "theta1"))
  eta <- as.numeric(predict(fit, type = "link", dpar = "theta1"))
  # one value per FITTED row, as every other prediction gives, and not
  # one per row of the data frame the rows were dropped from
  expect_length(rv, 297L)
  expect_length(eta, 297L)
  expect_false(anyNA(rv))
  expect_true(all(rv > 0 & rv < 1))
  expect_equal(rv, stats::plogis(eta), tolerance = 1e-12)

  se <- predict(fit, type = "response", dpar = "theta1", se.fit = TRUE)
  el <- predict(fit, type = "link", dpar = "theta1", se.fit = TRUE)
  expect_length(se$fit, 297L)
  expect_length(se$se.fit, 297L)
  # the delta method through the softmax: dp/deta = p (1 - p)
  expect_equal(as.numeric(se$se.fit),
               rv * (1 - rv) * as.numeric(el$se.fit),
               tolerance = 1e-12)
})

test_that("a mixing weight's response-scale SE is the one-predictor rule", {
  skip_on_cran()

  # DOCUMENTED LIMIT, pinned. The softmax response scale reports its
  # standard error through the family's `deriv` hook, p (1 - p), times
  # THIS predictor's own se_eta. At K = 2 that is exact: there is one
  # theta predictor and nothing else the weight depends on. At K >= 3
  # the softmax also moves with the other theta predictors
  # (dp1/dtheta2 = -p1 p2) and those terms are dropped, so the reported
  # standard error is WIDER than the joint delta method - conservative,
  # never optimistic. ?mixture and ?predict.frmtmb_fit say so; this
  # block is the measurement behind that sentence.
  set.seed(6)
  n <- 600
  x <- stats::rnorm(n)
  e1 <- 0.4 + 0.6 * x
  e2 <- -0.2 + 0.3 * x
  den <- exp(e1) + exp(e2) + 1
  pr <- cbind(exp(e1), exp(e2), 1) / den
  cl <- apply(pr, 1, function(q) sample.int(3, 1, prob = q))
  mu <- c(-2, 1, 4)[cl]
  d <- data.frame(x = x, y = stats::rnorm(n, mu, 0.7))
  fit <- frm(bf(y ~ 1, theta1 ~ x, theta2 ~ x) +
               mixture(gaussian(), gaussian(), gaussian()), data = d)

  eta <- cbind(as.numeric(predict(fit, type = "link", dpar = "theta1")),
               as.numeric(predict(fit, type = "link", dpar = "theta2")))
  den2 <- 1 + exp(eta[, 1]) + exp(eta[, 2])
  p1 <- exp(eta[, 1]) / den2
  p2 <- exp(eta[, 2]) / den2
  X <- cbind(1, d$x)
  V <- vcov(fit)
  cn <- c("theta1_(Intercept)", "theta1_x", "theta2_(Intercept)", "theta2_x")
  skip_if_not(all(cn %in% rownames(V)))

  rep_se <- as.numeric(predict(fit, type = "response", dpar = "theta1",
                               se.fit = TRUE)$se.fit)
  # the one-predictor rule, by hand: it IS what is reported
  Gown <- p1 * (1 - p1) * X
  Vown <- V[cn[1:2], cn[1:2], drop = FALSE]
  one <- sqrt(pmax(rowSums((Gown %*% Vown) * Gown), 0))
  expect_lt(max(abs(rep_se / one - 1)), 1e-10)

  # the joint delta over BOTH theta predictors is smaller, everywhere
  G <- cbind(Gown, -p1 * p2 * X)
  Vt <- V[cn, cn, drop = FALSE]
  joint <- sqrt(pmax(rowSums((G %*% Vt) * G), 0))
  expect_true(all(rep_se >= joint))
  # measured 5.5% to 26.1% wide on this fit; the review measured
  # 8.9% to 17.3% on its own. The pin is the DIRECTION and an order of
  # magnitude, not a number that would move with the data.
  expect_lt(max(rep_se / joint - 1), 0.6)
  expect_gt(max(rep_se / joint - 1), 0.01)

  # and at K = 2 the same rule is exact, which is why this is a limit
  # rather than a defect in the mechanism
  set.seed(9)
  d2 <- data.frame(x = stats::rnorm(400))
  p1b <- stats::plogis(0.8 + 0.5 * d2$x)
  d2$y <- ifelse(stats::rbinom(400, 1, p1b), stats::rnorm(400, -1, 0.6),
                 stats::rnorm(400, 2, 0.9))
  f2 <- frm(bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian()),
            data = d2)
  r2 <- as.numeric(predict(f2, type = "response", dpar = "theta1",
                           se.fit = TRUE)$se.fit)
  q <- as.numeric(predict(f2, type = "response", dpar = "theta1"))
  V2 <- vcov(f2)
  c2 <- c("theta1_(Intercept)", "theta1_x")
  skip_if_not(all(c2 %in% rownames(V2)))
  G2 <- q * (1 - q) * cbind(1, d2$x)
  ex <- sqrt(pmax(rowSums((G2 %*% V2[c2, c2, drop = FALSE]) * G2), 0))
  expect_lt(max(abs(r2 / ex - 1)), 1e-10)
})
