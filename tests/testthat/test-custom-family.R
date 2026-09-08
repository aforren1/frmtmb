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

# --- the extension API's hardening (frmtmb.eam findings 2, 3, 5) ------
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

test_that("required_aterms takes a character vector or a list of them", {
  bad <- function(req) {
    custom_family("bad_req", dpars = "mu", links = list(mu = "log"),
                  lpdf = function(y, dpars, aterms) y,
                  required_aterms = req)
  }
  expect_error(bad(list(1)), "a character vector for the values it needs")
  expect_error(bad(42), "a character vector for the values it needs")
  expect_error(bad(c("vint1", NA)), "missing or empty")
  # the two legal spellings: a conjunction, and a group of alternatives
  expect_s3_class(bad(c("vint1", "vint2")), "frmtmb_family")
  expect_s3_class(bad(list("vint1", c("dec", "vint1"))), "frmtmb_family")
})

# A family whose per-row datum arrives under either of two spellings.
# That is the shape required_aterms could not express: the refusal had
# to be hand-rolled inside valid_y, where it fires after the frame is
# built and where nothing about the framework can check it.
either_fam <- function(req) {
  custom_family(
    "either", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) RTMB::dpois(y, dpars$mu, log = TRUE),
    required_aterms = req, type = "discrete")
}

either_dat <- function(n = 60) {
  set.seed(11)
  dd <- data.frame(size = 3L, z = 1.5, x = stats::rnorm(n))
  dd$y <- stats::rpois(n, exp(0.4 + 0.3 * dd$x))
  dd
}

test_that("an any-of group is met by either alternative", {
  dd <- either_dat()
  req <- list(c("vint1", "vreal1"))
  expect_s3_class(frm(bf(y | vint(size) ~ x) + either_fam(req), data = dd),
                  "frmtmb_fit")
  expect_s3_class(frm(bf(y | vreal(z) ~ x) + either_fam(req), data = dd),
                  "frmtmb_fit")
})

test_that("an unmet any-of group is refused by naming the alternatives", {
  dd <- either_dat()
  err <- tryCatch(frm(bf(y ~ x) + either_fam(list(c("dec", "vint1"))),
                      data = dd),
                  error = conditionMessage)
  expect_match(err, "the density needs one of `dec` or `vint1`",
               fixed = TRUE)
  # the example formula writes one alternative, because a formula has
  # to pick one to be a formula
  expect_match(err, "y | dec(<column>) ~ ...", fixed = TRUE)
})

test_that("all-of and any-of mix in one declaration", {
  dd <- either_dat()
  req <- list(c("dec", "vint1"), "vreal1")
  # nothing supplied: both requirements are reported, the plain one
  # first so the sentence does not trail a bare name off a choice
  err <- tryCatch(frm(bf(y ~ x) + either_fam(req), data = dd),
                  error = conditionMessage)
  expect_match(err, "needs `vreal1`, one of `dec` or `vint1`", fixed = TRUE)
  expect_match(err, "y | vreal(<column>) + dec(<column>) ~ ...",
               fixed = TRUE)
  # the choice met, the plain requirement still missing
  err2 <- tryCatch(frm(bf(y | vint(size) ~ x) + either_fam(req), data = dd),
                   error = conditionMessage)
  expect_match(err2, "the density needs `vreal1`,", fixed = TRUE)
  expect_false(grepl("one of", err2, fixed = TRUE))
  expect_s3_class(
    frm(bf(y | vint(size) + vreal(z) ~ x) + either_fam(req), data = dd),
    "frmtmb_fit")
})

# --- exclusive_aterms ------------------------------------------------
#
# The allow-list cannot close this one, because both spellings are
# legitimately on it. Measured before the argument existed: wiener()
# fitted `rt | dec(u) + vint(1 - u)` with a log-likelihood bit-identical
# to the dec()-only model, so the user said one thing twice, and
# inconsistently, and nothing complained.

test_that("two spellings of one datum are refused together", {
  dd <- either_dat()
  fam <- custom_family(
    "excl", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) RTMB::dpois(y, dpars$mu, log = TRUE),
    required_aterms = list(c("vint1", "vreal1")),
    exclusive_aterms = list(c("vint1", "vreal1")),
    type = "discrete")
  # either alone still fits: exclusivity is at MOST one, not exactly one
  expect_s3_class(frm(bf(y | vint(size) ~ x) + fam, data = dd),
                  "frmtmb_fit")
  expect_s3_class(frm(bf(y | vreal(z) ~ x) + fam, data = dd),
                  "frmtmb_fit")
  err <- tryCatch(frm(bf(y | vint(size) + vreal(z) ~ x) + fam, data = dd),
                  error = conditionMessage)
  expect_match(err, "excl: `vint1` and `vreal1` are spellings of the same",
               fixed = TRUE)
  # it says which to keep, and the group order is that precedence
  expect_match(err, "The density reads `vint1` and would ignore `vreal1`",
               fixed = TRUE)
  expect_match(err, "Keep vint(<column>) and drop vreal(<column>).",
               fixed = TRUE)
})

test_that("an any-of requirement and the same exclusive set mean exactly one", {
  dd <- either_dat()
  fam <- custom_family(
    "excl1", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) RTMB::dpois(y, dpars$mu, log = TRUE),
    required_aterms = list(c("vint1", "vreal1")),
    exclusive_aterms = list(c("vint1", "vreal1")),
    type = "discrete")
  # neither: the required check speaks
  expect_error(frm(bf(y ~ x) + fam, data = dd),
               "the density needs one of `vint1` or `vreal1`", fixed = TRUE)
  # both: the exclusivity check speaks, and it runs after the required
  # one, which is why the two messages never collide
  expect_error(frm(bf(y | vint(size) + vreal(z) ~ x) + fam, data = dd),
               "are spellings of the same datum", fixed = TRUE)
})

test_that("a term the family does not accept is still refused first", {
  # the allow-list runs LAST, so an unaccepted term reaches it only when
  # nothing more specific fired: exclusivity is more specific
  dd <- either_dat()
  fam <- custom_family(
    "excl2", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) RTMB::dpois(y, dpars$mu, log = TRUE),
    accepts_aterms = c("vint", "vreal"),
    exclusive_aterms = c("vint1", "vreal1"),
    type = "discrete")
  expect_error(frm(bf(y | vint(size) + vreal(z) ~ x) + fam, data = dd),
               "are spellings of the same datum", fixed = TRUE)
  expect_error(frm(bf(y | weights(z) ~ x) + fam, data = dd),
               "is not one this family reads", fixed = TRUE)
})

test_that("exclusive_aterms validates its own shape", {
  bad <- function(ex, req = character(0)) {
    custom_family("bad_excl", dpars = "mu", links = list(mu = "log"),
                  lpdf = function(y, dpars, aterms) y,
                  required_aterms = req, exclusive_aterms = ex)
  }
  expect_error(bad(list("dec")), "distinct value")
  expect_error(bad(list(c("dec", "dec"))), "distinct value")
  expect_error(bad(42), "at most one of each set may be supplied")
  expect_error(bad(list(c("dec", NA))), "missing or empty value")
  # a non-character set is refused, not coerced. The validator used to
  # run AFTER as.character(), so this was accepted and became
  # c("1", "2"): a set matching no term value, and a rule that could
  # never fire. required_aterms refuses the same shape.
  expect_error(bad(list(c(1, 2))), "a character vector of term values")
  expect_error(bad(list(c(TRUE, FALSE))), "a character vector of term values")
  # the refusal says WHICH set and why, not just the argument's class
  expect_error(bad(list(c("dec", "vint1"), c("a", NA))), "Set 2")
  # a bare character vector is ONE set, unlike required_aterms
  f <- bad(c("dec", "vint1"))
  expect_identical(f[["exclusive_aterms"]], list(c("dec", "vint1")))
  expect_identical(bad(list())[["exclusive_aterms"]], list())
  # a conjunction that demands both is a family nobody could fit
  expect_error(bad(list(c("dec", "vint1")), req = c("dec", "vint1")),
               "no model could satisfy the family")
  # the same names as an ANY-of group are fine: together they read
  # "exactly one"
  expect_s3_class(bad(list(c("dec", "vint1")), req = list(c("dec", "vint1"))),
                  "frmtmb_family")
})

test_that("a mixture keeps the exclusive sets every component shares", {
  # Nothing covered mixture composition at all, and the first
  # implementation compared whole sets with setequal(), so a component
  # declaring a SUPERSET dropped the rule and the mixture then accepted
  # both spellings together: the failure was open, in exactly the case
  # this argument exists to close.
  mk <- function(ex) {
    custom_family("cmp", dpars = "mu", links = list(mu = "identity"),
                  lpdf = function(y, dpars, aterms) {
                    RTMB::dnorm(y, dpars$mu, 1, log = TRUE)
                  },
                  exclusive_aterms = ex)
  }
  ex_of <- function(...) {
    lapply(mixture(...)[["exclusive_aterms"]], sort)
  }
  pair <- c("dec", "vint1")
  # identical sets, and the same sets written in the other order
  expect_identical(ex_of(mk(list(pair)), mk(list(pair))), list(pair))
  expect_identical(ex_of(mk(list(pair)), mk(list(rev(pair)))), list(pair))
  # a SUPERSET still treats dec and vint1 as one datum, so the shared
  # pair survives rather than the whole rule vanishing
  expect_identical(
    ex_of(mk(list(pair)), mk(list(c("dec", "vint1", "vint2")))),
    list(pair))
  # one component declaring none, and disjoint sets, both keep nothing
  expect_length(ex_of(mk(list(pair)), mk(list())), 0L)
  expect_length(ex_of(mk(list(pair)), mk(list(c("vreal1", "vreal2")))), 0L)
  # a component that SPLITS a set keeps the half it still shares
  expect_identical(
    ex_of(mk(list(c("dec", "vint1", "vint2"))),
          mk(list(c("dec", "vint1"), c("vint2", "vreal1")))),
    list(pair))
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

# ---------------------------------------------------------------------
# se() is gated on the family's DECLARATION, not on its name. The term
# is the one whose whole effect is inside the density, so the family is
# the only thing that can say whether writing it changes anything.
# ---------------------------------------------------------------------

se_reader <- function(declare = c("accepts", "required", "neither")) {
  declare <- match.arg(declare)
  frmtmb_family(
    "se_reader",
    dpars = "mu",
    links = list(mu = "identity"),
    type = "continuous",
    accepts_aterms = if (declare == "accepts") c("weights", "se"),
    required_aterms = if (declare == "required") "se" else character(0),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], aterms[["se"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y)))
}

se_data <- function(seed = 909, n = 80) {
  set.seed(seed)
  s <- runif(n, 0.1, 0.6)
  data.frame(v = rnorm(n, 1, s), s = s)
}

test_that("a family that declares se() is given it", {
  d <- se_data()
  for (how in c("accepts", "required")) {
    fit <- frm(bf(v | se(s) ~ 1) + se_reader(how), data = d)
    # the known-variance mean is the inverse-variance weighted one
    w <- 1 / d$s^2
    expect_equal(unname(fixef(fit)$mu), sum(w * d$v) / sum(w),
                 tolerance = 1e-8, label = how)
    expect_equal(as.numeric(logLik(fit)),
                 sum(stats::dnorm(d$v, sum(w * d$v) / sum(w), d$s,
                                  log = TRUE)),
                 tolerance = 1e-8)
  }
})

test_that("a family that declares nothing is refused se(), by name", {
  d <- se_data()
  expect_error(frm(bf(v | se(s) ~ 1) + se_reader("neither"), data = d),
               "se_reader.*does not declare that it does")
  # the refusal says how to opt in
  expect_error(frm(bf(v | se(s) ~ 1) + se_reader("neither"), data = d),
               "accepts_aterms")
})

test_that("required_aterms = \"se\" also refuses a model without it", {
  d <- se_data()
  expect_error(frm(bf(v ~ 1) + se_reader("required"), data = d),
               "the density needs `se`")
})

test_that("the built-in families that read se() are unchanged", {
  d <- se_data()
  fit <- frm(bf(v | se(s) ~ 1) + gaussian(), data = d)
  ref <- frm(bf(v | se(s) ~ 1) + se_reader("accepts"), data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-8)
  expect_equal(sigma(fit), 0)
  expect_error(frm(bf(v | se(s) ~ 1) + Gamma(), data = d),
               "does not declare that it does")
})

test_that("a mixture does not inherit se() from its components", {
  d <- se_data()
  # a component reads the term, but the dpar se() maps out is named
  # sigma and a mixture's are sigma1 and sigma2, so the fit would be
  # unidentified rather than wrong: refused instead
  expect_error(frm(bf(v | se(s) ~ 1),
                   family = mixture(gaussian, gaussian), data = d),
               "does not declare that it does")
})

# ---------------------------------------------------------------------
# The other half of se(): the residual scale it replaces has to stop
# being free. The core maps out the dpar the convention names, `sigma`,
# so a declaring family whose scale is called anything else would be
# left with a flat direction and a NaN standard error.
# ---------------------------------------------------------------------

se_scale_fam <- function(scale_name) {
  force(scale_name)
  # by NAME, never by position: `dpars` also carries the `.eta_<dpar>`
  # linear predictors, so dpars[[2L]] is not the second declared dpar
  frmtmb_family(
    paste0("scale_", scale_name),
    accepts_aterms = c("weights", "se"),
    dpars = c("mu", scale_name),
    links = stats::setNames(list("identity", "log"), c("mu", scale_name)),
    lpdf = function(y, dpars, aterms) {
      sd <- if (is.null(aterms[["se"]])) {
        dpars[[scale_name]]
      } else if (isTRUE(aterms[["se_sigma"]])) {
        sqrt(dpars[[scale_name]]^2 + aterms[["se"]]^2)
      } else {
        aterms[["se"]]
      }
      RTMB::dnorm(y, dpars[["mu"]], sd, log = TRUE)
    },
    init_dpars = stats::setNames(
      list(function(y, aterms) mean(y),
           function(y, aterms) stats::sd(y)),
      c("mu", scale_name)))
}

se_scale_data <- function(seed = 5, n = 120, extra = 0) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), sev = runif(n, 0.2, 0.8))
  d$y <- 1 + 0.7 * d$x + rnorm(n, 0, sqrt(extra^2 + d$sev^2))
  d
}

test_that("a declared se() with an unmappable scale dpar is refused", {
  d <- se_scale_data()
  expect_error(frm(bf(y | se(sev) ~ x) + se_scale_fam("tau"), data = d),
               "has no `sigma`")
  # the message names the dpar and all three ways out
  msg <- tryCatch(frm(bf(y | se(sev) ~ x) + se_scale_fam("tau"),
                      data = d),
                  error = conditionMessage)
  expect_match(msg, "`tau`")
  expect_match(msg, "name the scale `sigma`")
  expect_match(msg, "pin it in the formula")
  expect_match(msg, "sigma = TRUE")
})

test_that("each way out of the unmappable scale actually works", {
  d <- se_scale_data()
  # 1. name the scale sigma: the core maps it out
  f1 <- frm(bf(y | se(sev) ~ x) + se_scale_fam("sigma"), data = d)
  expect_true(all(is.finite(sqrt(diag(vcov(f1))))))
  # 2. pin it in the formula
  f2 <- frm(bf(y | se(sev) ~ x, tau = 1) + se_scale_fam("tau"), data = d)
  expect_true(all(is.finite(sqrt(diag(vcov(f2))))))
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
  # 3. sigma = TRUE, where the scale stays estimated alongside the
  # known one. Needs data carrying variance beyond se, or the scale
  # sits on its boundary at zero and is flat for a real reason.
  dq <- se_scale_data(seed = 21, n = 400, extra = 0.9)
  f3 <- frm(bf(y | se(sev, sigma = TRUE) ~ x) + se_scale_fam("tau"),
            data = dq)
  expect_true(all(is.finite(sqrt(diag(vcov(f3))))))
  expect_equal(exp(unname(fixef(f3)$tau)), 0.9, tolerance = 0.05)
  # and it is the same model as the one whose scale is named sigma
  f4 <- frm(bf(y | se(sev, sigma = TRUE) ~ x) + se_scale_fam("sigma"),
            data = dq)
  expect_equal(as.numeric(logLik(f3)), as.numeric(logLik(f4)),
               tolerance = 1e-10)
})

test_that("a family whose whole scale is the known one is not refused", {
  # no dpar beyond the primaries, so there is nothing to map out; this
  # is the shape bcm_gaussian_probit() uses
  d <- se_scale_data()
  fit <- frm(bf(y | se(sev) ~ x) + se_reader("accepts"), data = d)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

test_that("the built-in se() families are untouched by the scale rule", {
  d <- se_scale_data()
  fg <- frm(y | se(sev) ~ x, data = d)
  expect_equal(sigma(fg), 0)
  expect_true(all(is.finite(sqrt(diag(vcov(fg))))))
  # student keeps nu free alongside a known se, as it always has
  fs <- frm(bf(y | se(sev) ~ x) + student(), data = d)
  expect_true("nu" %in% names(fixef(fs)))
  expect_true(all(is.finite(sqrt(diag(vcov(fg))))))
  # sigma = TRUE keeps the estimated scale
  fq <- frm(y | se(sev, sigma = TRUE) ~ x, data = d)
  expect_gt(sigma(fq), 0)
})

test_that("the se() refusal names a declaration a caller can write", {
  d <- se_scale_data()
  msg <- tryCatch(frm(bf(y | se(sev) ~ x) + se_reader("neither"),
                      data = d),
                  error = conditionMessage)
  # it must not send an author to an unexported helper
  expect_false(grepl("resid_sd", msg, fixed = TRUE))
  expect_match(msg, "accepts_aterms")
  expect_match(msg, "required_aterms")
  # and it states the sigma = TRUE convention inline
  expect_match(msg, "se_sigma")
})
