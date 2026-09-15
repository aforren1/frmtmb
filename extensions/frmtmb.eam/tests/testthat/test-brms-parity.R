## The parameterization claim: a model written for brms::wiener() means
## the same thing here. Checked without fitting anything, because the
## claim is about names, coding and the boundary convention, not about
## either package's arithmetic.

test_that("the dpar names are brms's", {
  skip_if_not_installed("brms")
  expect_setequal(brms::wiener()$dpars, wiener(max_ndt = 1)$dpars)
  expect_equal(wiener(max_ndt = 1)$dpars, c("mu", "bs", "ndt", "bias"))
})

test_that("the links match brms, except where ndt's bound lives", {
  skip_if_not_installed("brms")
  bw <- brms::wiener()
  ow <- wiener(max_ndt = 0.3)
  nm <- function(x) if (is.list(x)) x$name else x$name
  expect_equal(nm(ow$links$bs), bw$link_bs)       # log
  expect_equal(nm(ow$links$bias), bw$link_bias)   # logit
  expect_equal(nm(ow$links$mu), bw$link)          # identity

  # ndt is the deliberate difference. With no formula on `ndt`, brms
  # declares a log link and then bounds the parameter in Stan's
  # parameter block as `real<lower=0,upper=min_Y> ndt`. frmtmb has no
  # parameter block to bound, so the same constraint is carried by the
  # link itself. Give `ndt` a formula and brms drops that bound
  # entirely; the test below reads both spellings off brms's own
  # generated program.
  expect_equal(nm(ow$links$ndt), "scaled_logit")
  # over the working range the bound holds strictly
  eta <- c(-40, -5, 0, 5, 30)
  expect_true(all(ow$links$ndt$linkinv(eta) < 0.3))
  expect_true(all(ow$links$ndt$linkinv(eta) > 0))
  # past about 37 the logit saturates in double precision and the value
  # rounds to max_ndt exactly. The density is what keeps the optimizer
  # away from there: it falls off a cliff long before the link does, so
  # this is a known edge rather than a reachable state.
  expect_equal(ow$links$ndt$linkinv(40), 0.3)
  expect_lt(ddm_lpdf_lower(1e-8, 1, 1.4, 0.5), -1e6)
})

test_that("brms bounds ndt only where ndt has no formula", {
  skip_if_not_installed("brms")
  # WHY this is asserted rather than remembered. The bound is the whole
  # difference between the two parameterizations, and the sentence above
  # describes only brms's SCALAR spelling. A hierarchical Wiener writes
  # `ndt ~ 1 + (1 | s)`, and there brms emits a plain log link with no
  # upper bound and no uniform prior: what keeps the non-decision time
  # under the response time is `wiener_diffusion_lpdf` itself rejecting
  # a row whose decision time is not positive, which is a PER-ROW
  # constraint and so, through a subject's own rows, a per-subject one.
  # That is `ndt_group()`'s constraint, not the global one this family
  # used to carry, and it decides what a brms cross-check of item 2.1
  # would even be comparing. See dev/eamhier-findings.md.
  set.seed(3)
  d <- ddm_simulate(120, mu = 0.8, bs = 1.4, ndt = 0.25)
  d$s <- factor(rep(1:4, each = 30))
  scalar <- as.character(brms::make_stancode(
    brms::bf(rt | dec(upper) ~ 1, bs ~ 1, bias = 0.5),
    family = brms::wiener(), data = d))
  formula <- as.character(brms::make_stancode(
    brms::bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1 + (1 | s),
             bias = 0.5),
    family = brms::wiener(), data = d))
  # the scalar spelling: bounded in the parameter block, and given the
  # uniform prior that bound implies
  expect_match(scalar, "real<lower=0,upper=min_Y> ndt")
  expect_match(scalar, "uniform_lpdf(ndt | 0, min_Y)", fixed = TRUE)
  # the formula spelling: neither. Both of the assertions above fail on
  # this program, which is the point of making them separately.
  expect_false(grepl("upper=min_Y> ndt", formula, fixed = TRUE))
  expect_false(grepl("uniform_lpdf(ndt", formula, fixed = TRUE))
  expect_match(formula, "ndt = exp(ndt)", fixed = TRUE)
})

test_that("dec() and vint() carry the same 0/1 coding", {
  skip_if_not_installed("brms")
  skip_if_not_installed("RWiener")
  set.seed(1)
  d <- ddm_simulate(40, mu = 0.8, bs = 1.4, ndt = 0.25)
  d$dec <- d$upper
  d$x <- stats::rnorm(40)
  sd <- brms::make_standata(brms::bf(rt | dec(dec) ~ x),
                            family = brms::wiener(), data = d)
  expect_equal(as.integer(sd$dec), as.integer(d$upper))
})

test_that("dec == 1 is the upper boundary in both packages", {
  skip_if_not_installed("RWiener")
  # brms's generated Stan reads
  #   if (dec == 1) wiener_lpdf(y | alpha, tau, beta, delta)
  #   else          wiener_lpdf(y | alpha, tau, 1 - beta, -delta)
  # and Stan's wiener_lpdf is the density at the UPPER boundary. So
  # dec == 1 is the upper boundary there. It is here too, and RWiener
  # is the arbiter.
  t <- 0.6; v <- 0.9; a <- 1.3; w <- 0.4
  expect_equal(ddm_lpdf_both(t, v, a, w, 1),
               log(RWiener::dwiener(t + 1e-12, a, 1e-12, w, v,
                                    resp = "upper")),
               tolerance = 1e-10)
  expect_equal(ddm_lpdf_both(t, v, a, w, 0),
               log(RWiener::dwiener(t + 1e-12, a, 1e-12, w, v,
                                    resp = "lower")),
               tolerance = 1e-10)
})
