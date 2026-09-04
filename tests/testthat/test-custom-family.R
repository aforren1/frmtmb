my_poisson <- function() {
  custom_family(
    "my_poisson",
    dpars = "mu",
    links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
    },
    init_dpars = list(mu = function(y, aterms) mean(y) + 0.1),
    type = "discrete",
    post = list(mean_fn = function(dpars, aterms) dpars$mu),
    sim = function(dpars, aterms, n) stats::rpois(n, dpars$mu)
  )
}

test_that("a hand-written custom family matches the built-in", {
  set.seed(121)
  dd <- data.frame(x = rnorm(300), g = factor(rep(1:15, 20)))
  dd$y <- rpois(300, exp(0.5 + 0.3 * dd$x + rnorm(15, 0, 0.4)[dd$g]))

  fit_c <- frm(bf(y ~ x + (1 | g)) + my_poisson(), data = dd)
  fit_b <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit_c)) - as.numeric(logLik(fit_b))),
            1e-8)
  expect_vector_equal(fixef(fit_c)$mu, fixef(fit_b)$mu, tol = 1e-6)
  # the whole post-processing stack works on the custom family
  expect_equal(fitted(fit_c), fitted(fit_b), tolerance = 1e-6)
  s <- simulate(fit_c, nsim = 1, seed = 1)
  expect_true(all(s$sim_1 >= 0))
})

test_that("check_custom_family passes a correct lpdf", {
  expect_true(check_custom_family(
    my_poisson(), y = rpois(50, 3),
    dpars = list(mu = rep(2.5, 50))
  ))
})

test_that("check_custom_family catches tape-unsafe code", {
  bad <- custom_family(
    "bad",
    dpars = "mu",
    links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      # base::matrix strips the advector class: values silently wrong
      m <- matrix(dpars$mu, ncol = 1)
      y * log(m[, 1]) - m[, 1] - lgamma(y + 1)
    }
  )
  expect_error(check_custom_family(bad, y = rpois(50, 3),
                                   dpars = list(mu = rep(2.5, 50))),
               "tape|differently|gradient")
})

test_that("a custom lpdf needs no ADoverload boilerplate of its own", {
  # numeric-first c() and pad[i] <- on advectors, no ADoverload lines:
  # frmtmb_family() splices the overloads onto the body
  cf <- custom_family(
    "gauss_bare", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      pad <- c(0, dpars$mu)
      pad[1] <- pad[2]
      RTMB::dnorm(y, pad[-1], dpars$sigma, log = TRUE)
    })
  # the wrap is visible on the stored function, and a function that
  # binds the overloads itself is left untouched
  expect_true("ADoverload" %in% all.names(body(cf$lpdf)))
  own <- function(y, dpars, aterms) {
    "c" <- RTMB::ADoverload("c")
    RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE)
  }
  cf2 <- custom_family("gauss_own", dpars = c("mu", "sigma"),
                       links = list(mu = "identity", sigma = "log"),
                       lpdf = own)
  expect_identical(cf2$lpdf, own)
  set.seed(2)
  dd <- data.frame(x = stats::rnorm(120))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x, 1.3)
  fc <- frm(bf(y ~ x), family = cf, data = dd)
  ref <- frm(bf(y ~ x), family = gaussian(), data = dd)
  expect_equal(as.numeric(logLik(fc)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_equal(unname(unlist(fixef(fc)$mu)),
               unname(unlist(fixef(ref)$mu)), tolerance = 1e-4)
})

# --- the extension API's hardening (frmtmb.ddm findings 2, 3, 5) ------
#
# Each of these closes a route by which an extension-written family
# produced a wrong answer or a distant error rather than a refusal.

# A family whose density indexes a vint() payload it never checks for.
# Without a declaration this is the silent case: an absent vint1 is
# NULL, the arithmetic gives numeric(0), and the fit returns.
needs_vint <- function(required = character(0)) {
  custom_family(
    "needs_vint", dpars = "mu", links = list(mu = "logit"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
    },
    required_aterms = required,
    type = "discrete",
    post = list(mean_fn = function(dpars, aterms) {
      dpars$mu * aterms$vint1
    })
  )
}

binom_dat <- function(n = 80) {
  set.seed(407)
  dd <- data.frame(size = 5L, x = stats::rnorm(n))
  dd$y <- stats::rbinom(n, 5L, stats::plogis(0.3 + 0.5 * dd$x))
  dd
}

test_that("required_aterms refuses a missing addition term by name", {
  dd <- binom_dat()
  expect_error(
    frm(bf(y ~ x) + needs_vint("vint1"), data = dd),
    "needs_vint: the density needs `vint1`"
  )
  # the same family with the term supplied fits
  fit <- frm(bf(y | vint(size) ~ x) + needs_vint("vint1"), data = dd)
  expect_s3_class(fit, "frmtmb_fit")
  # the refusal names the term the way the user has to write it
  err <- tryCatch(frm(bf(y ~ x) + needs_vint(c("vint1", "vint2")),
                      data = dd),
                  error = conditionMessage)
  expect_match(err, "vint(<column>) + vint(..., <column>)", fixed = TRUE)
})

test_that("required_aterms is checked before the family's own valid_y", {
  # valid_y reading the absent term would otherwise fail first, with a
  # message about NULL rather than about the missing term
  fam <- custom_family(
    "vint_valid_y", dpars = "mu", links = list(mu = "logit"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dbinom(y, aterms$vint1, dpars$mu, log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y > aterms$vint1)) stop("y exceeds the trial count")
      invisible(NULL)
    },
    required_aterms = "vint1", type = "discrete")
  expect_error(frm(bf(y ~ x) + fam, data = binom_dat()),
               "vint_valid_y: the density needs")
})

test_that("a zero-length log-density is refused even undeclared", {
  # The backstop: no required_aterms, so nothing but the length of the
  # density's own result says anything is wrong. The aterm enters as
  # plain data before it meets an advector, which is the reproduced
  # hazard: NULL times a number is numeric(0), and an advector plus a
  # zero-length vector is zero-length rather than an error.
  undeclared <- custom_family(
    "undeclared_vint", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      off <- aterms$vint1 * 0
      RTMB::dpois(y, dpars$mu, log = TRUE) + off
    },
    type = "discrete")
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(80))
  dd$y <- stats::rpois(80, exp(0.5 + 0.3 * dd$x))
  expect_error(frm(bf(y ~ x) + undeclared, data = dd),
               "returned no values for 80 observations")
})

test_that("required_aterms takes a character vector only", {
  expect_error(
    custom_family("bad_req", dpars = "mu", links = list(mu = "log"),
                  lpdf = function(y, dpars, aterms) y,
                  required_aterms = list("vint1")),
    "names the addition terms"
  )
})

test_that("family_finalize derives a link from the response", {
  # the shifted-family problem: the ndt link's upper bound is min(y),
  # which the family cannot know until frm() has the data
  shifted <- function() {
    custom_family(
      "shifted_exp", dpars = c("mu", "ndt"),
      links = list(mu = "log", ndt = "log"),
      lpdf = function(y, dpars, aterms) {
        RTMB::dexp(y - dpars$ndt, 1 / dpars$mu, log = TRUE)
      },
      init_dpars = list(mu = function(y, aterms) mean(y) - min(y) / 2,
                        ndt = function(y, aterms) min(y) / 2),
      family_finalize = function(fam, y, aterms) {
        ub <- min(y)
        fam$links$ndt <- list(
          name = "ndt_bounded",
          linkfun = function(mu) log(mu / (ub - mu)),
          linkinv = function(eta) ub / (1 + exp(-eta)),
          mu_eta = function(eta) {
            p <- 1 / (1 + exp(-eta))
            ub * p * (1 - p)
          }
        )
        fam
      },
      type = "continuous",
      post = list(mean_fn = function(dpars, aterms) dpars$mu + dpars$ndt)
    )
  }
  set.seed(11)
  dd <- data.frame(x = stats::rnorm(200))
  dd$y <- 0.3 + stats::rexp(200, 1 / exp(0.2 + 0.3 * dd$x))
  fit <- frm(bf(y ~ x) + shifted(), data = dd)
  # the finalized link is what the fit reports and what predict() uses,
  # not the "log" the family was constructed with
  expect_identical(family(fit)$links$ndt$name, "ndt_bounded")
  ndt <- unname(predict(fit, type = "response", dpar = "ndt")[1])
  expect_lt(ndt, min(dd$y))
  expect_gt(ndt, 0)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("family_finalize must return a family", {
  fam <- custom_family(
    "bad_finalize", dpars = "mu", links = list(mu = "identity"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars$mu, 1, log = TRUE)
    },
    family_finalize = function(fam, y, aterms) list(links = fam$links))
  set.seed(3)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- stats::rnorm(40, dd$x)
  expect_error(frm(bf(y ~ x) + fam, data = dd),
               "family_finalize() must return a family object",
               fixed = TRUE)
})

test_that("family_finalize is not a function", {
  expect_error(
    custom_family("nf", dpars = "mu", links = list(mu = "log"),
                  lpdf = function(y, dpars, aterms) y,
                  family_finalize = "later"),
    "must be a function"
  )
})

test_that("a custom link object is validated at family construction", {
  # the old failure was inside predict(se.fit = TRUE), far from the
  # family that caused it
  half <- list(name = "half", linkfun = function(mu) mu,
               linkinv = function(eta) eta)
  expect_error(
    custom_family("no_mu_eta", dpars = c("mu", "sigma"),
                  links = list(mu = "identity", sigma = half),
                  lpdf = function(y, dpars, aterms) y),
    "dpar 'sigma'.*`mu_eta`"
  )
  expect_error(
    custom_family("bad_link_fn", dpars = "mu",
                  links = list(mu = list(name = "x", linkfun = 1,
                                         linkinv = identity,
                                         mu_eta = identity)),
                  lpdf = function(y, dpars, aterms) y),
    "non-function"
  )
  expect_error(
    custom_family("bad_link_name", dpars = "mu",
                  links = list(mu = list(name = 42, linkfun = identity,
                                         linkinv = identity,
                                         mu_eta = identity)),
                  lpdf = function(y, dpars, aterms) y),
    "must name itself with a single string"
  )
  # a list missing several fields names all of them at once
  expect_error(
    custom_family("empty_link", dpars = "mu", links = list(mu = list()),
                  lpdf = function(y, dpars, aterms) y),
    "has no `name`, `linkfun`, `linkinv`, `mu_eta`"
  )
})

test_that("a non-finite init_dpars value is reported, not dropped", {
  fam <- custom_family(
    "init_out_of_range", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars$mu, log = TRUE)
    },
    # a mean of zero is outside the log link's range
    init_dpars = list(mu = function(y, aterms) 0),
    type = "discrete")
  set.seed(5)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- stats::rpois(60, exp(0.5 + 0.3 * dd$x))
  expect_warning(frm(bf(y ~ x) + fam, data = dd),
                 "Starting value 0 for mu is -Inf through its log link")
})
