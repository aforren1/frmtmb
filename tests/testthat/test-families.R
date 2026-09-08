test_that("Gamma matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(21)
  n <- 400
  x <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  y <- rgamma(n, shape = 2, scale = mu / 2)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + Gamma(link = "log"), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = Gamma(link = "log"), data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
})

test_that("lognormal equals lm on the log scale plus the Jacobian", {
  set.seed(22)
  n <- 300
  x <- rnorm(n)
  y <- exp(rnorm(n, 0.5 + 0.7 * x, 0.4))
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + lognormal(), data = dd)
  ref <- lm(log(y) ~ x, dd)
  expect_lt(abs(as.numeric(logLik(fit)) -
                  (as.numeric(logLik(ref)) - sum(log(y)))), 1e-6)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-5)
})

test_that("student matches MASS::fitdistr", {
  skip_if_not_installed("MASS")
  set.seed(23)
  y <- 3 + 1.5 * rt(800, df = 5)
  dd <- data.frame(y = y)
  fit <- frm(bf(y ~ 1) + student(), data = dd)
  ref <- suppressWarnings(MASS::fitdistr(y, "t"))
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
  expect_lt(abs(fixef(fit)$mu[[1]] - ref$estimate[["m"]]), 1e-3)
  est_sigma <- exp(fixef(fit)$sigma[[1]])
  expect_lt(abs(est_sigma - ref$estimate[["s"]]), 1e-3)
  est_nu <- 1 + exp(fixef(fit)$nu[[1]])
  expect_lt(abs(est_nu - ref$estimate[["df"]]) / ref$estimate[["df"]], 1e-2)
})

# The two tests below are about ONE thing: the student log density has
# to keep its digits when the degrees of freedom run to infinity, which
# is what data with no heavy tails asks of it. `RTMB::dt()` sends a
# double straight to `stats::dt()` and an AD number to a tape that
# subtracts two `lgamma()` values agreeing in every leading digit, so
# the accurate branch is the one no fit uses. dev/remlopt-findings.md
# carries the measurement.

test_that("the student density keeps its digits as nu runs off", {
  set.seed(11)
  z <- stats::rnorm(120)
  fam <- student()
  obj <- RTMB::MakeADFun(
    function(p) {
      -sum(fam$lpdf(z, list(mu = 0, sigma = 1,
                            nu = 1 + exp(p[["e"]])), list()))
    },
    list(e = 5), silent = TRUE)
  for (e in c(0, 2, 5, 10, 15, 20, 25, 30, 35)) {
    # stats::dt is the reference: it holds 1e-13 against a 300-bit
    # Rmpfr reference at every df tried, up to 1e50
    ref <- -sum(stats::dt(z, df = 1 + exp(e), log = TRUE))
    got <- as.numeric(obj$fn(e))
    # the tolerance is the answer's own magnitude in units of the
    # double that carries it, not a number chosen in advance
    expect_lt(abs(got - ref), 64 * .Machine$double.eps * abs(ref))
  }
})

test_that("a student objective has no false optimum in the large-nu tail", {
  set.seed(12)
  n <- 150
  # uniform errors have lighter tails than any student-t, so the
  # likelihood rises all the way to the gaussian limit: there is no
  # interior optimum in nu and no sign change in its gradient
  y <- 1 + (stats::runif(n) - 0.5) * 3.4
  fam <- student()
  obj <- RTMB::MakeADFun(
    function(p) {
      -sum(fam$lpdf(y, list(mu = p[["mu"]], sigma = exp(p[["ls"]]),
                            nu = 1 + exp(p[["e"]])), list()))
    },
    list(mu = mean(y), ls = log(stats::sd(y)), e = 5), silent = TRUE)
  base <- c(mean(y), log(stats::sd(y)))
  es <- seq(4, 34, by = 0.5)
  val <- vapply(es, function(e) as.numeric(obj$fn(c(base, e))), 0)
  grd <- vapply(es, function(e) as.numeric(obj$gr(c(base, e)))[3], 0)
  # every rise in the negative log likelihood is the density having
  # lost digits, so the ceiling is what a double can hold at this
  # objective's magnitude
  noise <- 64 * .Machine$double.eps * max(abs(val))
  expect_lt(max(c(0, diff(val))), noise)
  # The optimizer stops where the gradient turns uphill. Past the point
  # where nu is unidentified to the last bit the sign of a gradient
  # this small is meaningless, and the ladder steps by 0.5, so the same
  # floor bounds the derivative: what may not happen is an uphill
  # gradient BIGGER than the objective's own noise.
  expect_equal(sum(grd > noise), 0L)
})

test_that("negbinomial GLM matches MASS::glm.nb", {
  skip_if_not_installed("MASS")
  set.seed(24)
  n <- 500
  x <- rnorm(n)
  y <- rnbinom(n, size = 1.5, mu = exp(0.5 + 0.4 * x))
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + negbinomial(), data = dd)
  ref <- MASS::glm.nb(y ~ x, data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-5)
  expect_vector_equal(fixef(fit)$mu, coef(ref), tol = 1e-4)
  expect_lt(abs(exp(fixef(fit)$shape[[1]]) - ref$theta), 1e-2)
})

test_that("nbinom1 matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  set.seed(25)
  n <- 500
  x <- rnorm(n)
  mu <- exp(1 + 0.3 * x)
  phi <- 1.5
  y <- rnbinom(n, size = mu / phi, mu = mu)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + nbinom1(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::nbinom1, data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-5)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
})

test_that("beta matches glmmTMB beta_family", {
  skip_if_not_installed("glmmTMB")
  set.seed(26)
  n <- 400
  x <- rnorm(n)
  mu <- plogis(0.3 + 0.6 * x)
  phi <- 8
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + Beta(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::beta_family(),
                          data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-6)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-4)
  expect_lt(abs(exp(fixef(fit)$phi[[1]]) - glmmTMB::sigma(ref)), 1e-2)
})

test_that("tweedie matches glmmTMB", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")
  set.seed(27)
  n <- 400
  x <- rnorm(n)
  mu <- exp(0.8 + 0.3 * x)
  y <- mgcv::rTweedie(mu, p = 1.5, phi = 1.2)
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + tweedie(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x, family = glmmTMB::tweedie(), data = dd)
  expect_loglik_equal(fit, ref, tol = 1e-4)
  expect_vector_equal(fixef(fit)$mu, unname(glmmTMB::fixef(ref)$cond),
                      tol = 1e-3)
})

test_that("compois fits and recovers the mean model", {
  set.seed(28)
  n <- 120
  x <- rnorm(n)
  y <- rpois(n, exp(0.6 + 0.4 * x))   # nu = 1 reduces to poisson
  dd <- data.frame(y, x)
  fit <- frm(bf(y ~ x) + compois(), data = dd)
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_vector_equal(fixef(fit)$mu, c(0.6, 0.4), tol = 0.3)
  # poisson nested in compois: likelihood at optimum can't be worse
  ref <- glm(y ~ x, family = poisson, data = dd)
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(ref)) - 1e-6)
})

test_that("simulate() round-trips through family simulators", {
  set.seed(29)
  dd <- data.frame(x = rnorm(300), g = factor(rep(1:10, 30)))
  dd$y <- rpois(300, exp(0.5 + 0.3 * dd$x + rnorm(10, 0, 0.4)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  s1 <- simulate(fit, nsim = 3, seed = 1)
  s2 <- simulate(fit, nsim = 3, seed = 1)
  expect_identical(s1, s2)
  expect_identical(dim(s1), c(300L, 3L))
  expect_true(all(s1$sim_1 >= 0))

  # marginal simulation redraws random effects
  m1 <- simulate(fit, nsim = 1, seed = 2, re.form = NA)
  c1 <- simulate(fit, nsim = 1, seed = 2)
  expect_false(identical(m1, c1))
})

# --- links on distributional parameters other than the mean ----------

test_that("every constructor takes the dpar links brms allows", {
  # the sets are keyed on the parameter's support, and `nu` means two
  # different parameters: student's degrees of freedom take logm1,
  # compois's dispersion is an ordinary positive scale
  pos <- c("log", "identity", "softplus", "squareplus")
  spec <- list(
    student = list(sigma = pos, nu = c("logm1", "identity")),
    lognormal = list(sigma = pos), negbinomial = list(shape = pos),
    nbinom1 = list(phi = pos), Beta = list(phi = pos),
    tweedie = list(phi = pos), compois = list(nu = pos),
    weibull = list(shape = pos), huber = list(sigma = pos),
    beta_binomial = list(phi = pos), von_mises = list(kappa = pos),
    exgaussian = list(sigma = pos, beta = pos),
    shifted_lognormal = list(sigma = pos, ndt = pos),
    skew_normal = list(sigma = pos, alpha = pos),
    zero_inflated_poisson = list(zi = c("logit", "identity")),
    zero_inflated_binomial = list(zi = c("logit", "identity")),
    hurdle_poisson = list(hu = c("logit", "identity")),
    zero_inflated_negbinomial = list(shape = pos,
                                     zi = c("logit", "identity")),
    hurdle_gamma = list(shape = pos, hu = c("logit", "identity")),
    hurdle_lognormal = list(sigma = pos, hu = c("logit", "identity")),
    zero_inflated_beta = list(phi = pos, zi = c("logit", "identity")),
    asym_laplace = list(sigma = pos, quantile = c("logit", "identity")))
  for (cn in names(spec)) {
    ctor <- get(cn, envir = asNamespace("frmtmb"))
    for (dp in names(spec[[cn]])) {
      arg <- paste0("link_", dp)
      expect_true(arg %in% names(formals(ctor)),
                  label = paste(cn, "has", arg))
      for (l in spec[[cn]][[dp]]) {
        fam <- do.call(ctor, stats::setNames(list(l), arg))
        expect_identical(fam$links[[dp]]$name, l,
                         label = paste(cn, arg, l))
      }
      # a link outside the parameter's range is refused, by name
      expect_error(do.call(ctor, stats::setNames(list("tan_half"), arg)),
                   dp, fixed = TRUE)
    }
  }
})

test_that("frm_family() reaches the links the stats families cannot", {
  # stats::gaussian() has no link_sigma and stats::poisson() refuses
  # softplus, so these are unreachable through the constructors
  expect_identical(
    frm_family("gaussian", link_sigma = "softplus")$links$sigma$name,
    "softplus")
  expect_identical(frm_family("poisson", link = "softplus")$links$mu$name,
                   "softplus")
  expect_identical(
    frm_family("Gamma", link = "inverse", link_shape = "identity")$links$
      shape$name, "identity")
  expect_error(frm_family("gaussian", link_shape = "log"), "link_shape")
  expect_error(frm_family("nope"), "no family called")
  expect_error(frm_family("gaussian", "softplus"), "has to be named")
})

test_that("a dpar link that is not log keeps the density honest", {
  # log_dpar() used to hand the linear predictor back AS the log of the
  # dpar, which is right only on a log link. On a softplus shape that
  # would be a different, wrong, density; the fit is compared with the
  # log-link fit of the SAME model, which has to reach the same place.
  set.seed(4)
  n <- 300
  d <- data.frame(x = rnorm(n))
  d$y <- rnbinom(n, mu = exp(1 + 0.5 * d$x), size = 2)
  a <- frm(bf(y ~ x), family = negbinomial(link_shape = "log"), data = d)
  b <- frm(bf(y ~ x), family = negbinomial(link_shape = "softplus"),
           data = d)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-6)
  expect_equal(unlist(fixef(a)$mu), unlist(fixef(b)$mu), tolerance = 1e-5)
  # the shape itself, read back off each link, agrees
  sa <- exp(unlist(fixef(a))[["shape.(Intercept)"]])
  sb <- log1p(exp(unlist(fixef(b))[["shape.(Intercept)"]]))
  expect_equal(sa, sb, tolerance = 1e-4)
})

test_that("a gate link that is not logit keeps the density honest", {
  # the same hazard for gate_logs(), which assumed the logit
  set.seed(5)
  n <- 400
  d <- data.frame(x = rnorm(n))
  d$y <- ifelse(runif(n) < 0.25, 0L, rpois(n, exp(1 + 0.4 * d$x)))
  a <- frm(bf(y ~ x), family = zero_inflated_poisson(link_zi = "logit"),
           data = d)
  b <- frm(bf(y ~ x), family = zero_inflated_poisson(link_zi = "identity"),
           data = d)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)),
               tolerance = 1e-5)
  za <- plogis(unlist(fixef(a))[["zi.(Intercept)"]])
  zb <- unlist(fixef(b))[["zi.(Intercept)"]]
  expect_equal(za, zb, tolerance = 1e-4)
})

test_that("a dpar's prior goes through that dpar's own link", {
  # set_prior(class = "sigma") is a density on sigma itself, so the
  # objective carries |d sigma / d eta|. That Jacobian has to be the
  # dpar's OWN link: with a log link assumed, the softplus row below
  # would land on the log row instead of on its own reference.
  set.seed(7)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.7 * d$x + rnorm(n, 0, 1.4)
  X <- model.matrix(~ x, d)
  rss <- sum((d$y - X %*% qr.solve(X, d$y))^2)   # profile beta is OLS
  pr <- set_prior("normal(0, 0.5)", class = "sigma")
  got <- ref <- stats::setNames(numeric(4), c("log", "identity",
                                              "softplus", "squareplus"))
  for (l in names(got)) {
    lk <- frmtmb:::get_link(l)
    f <- frm(bf(y ~ x), family = frm_family("gaussian", link_sigma = l),
             data = d, prior = pr)
    got[l] <- lk$linkinv(unname(unlist(fixef(f))[["sigma.(Intercept)"]]))
    obj <- function(eta) {
      s <- lk$linkinv(eta)
      if (s <= 0) return(1e10)
      -(-n * log(s) - rss / (2 * s^2) +
        stats::dnorm(s, 0, 0.5, log = TRUE) + log(abs(lk$mu_eta(eta))))
    }
    ref[l] <- lk$linkinv(stats::optimize(obj, c(-8, 8), tol = 1e-12)$minimum)
  }
  expect_equal(unname(got), unname(ref), tolerance = 1e-5)
  # and the four are genuinely different places, so the test above is
  # not passing on a coincidence
  expect_gt(diff(range(got)), 1e-4)
})

test_that("summary() and print() name the link of every dpar", {
  set.seed(3)
  n <- 80
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + d$x + rnorm(n)
  f <- frm(bf(y ~ x), family = gaussian(), data = d)
  expect_match(paste(utils::capture.output(print(summary(f))),
                     collapse = "\n"),
               "Links: mu = identity; sigma = log", fixed = TRUE)
  expect_match(paste(utils::capture.output(print(f)), collapse = "\n"),
               "Links: mu = identity; sigma = log", fixed = TRUE)
  g <- frm(bf(y ~ x), family = frm_family("gaussian",
                                          link_sigma = "softplus"), data = d)
  expect_match(paste(utils::capture.output(print(summary(g))),
                     collapse = "\n"),
               "sigma = softplus", fixed = TRUE)
})
