## frmtmb.sample's default on this package's ndt. Without ndt_group()
## the bounded link's image is (0, fastest response), brms's interval,
## and the default is brms's uniform(0, min_Y). With ndt_group() ndt is
## a FRACTION in (0, 1) of each group's fastest response; the default
## is the uniform over that, and through frmtmb.sample's first lane
## build it was uniform(0, min_Y) on the fraction, which pinned the
## posterior below min_Y (0.29 to 0.34 against an ML fraction of 0.818,
## dev/reviews/2026-10-14-sampleprior.md, finding 9).

snd_skip <- function() {
  skip_on_cran()
  skip_if_not_installed("frmtmb.sample")
  skip_if_not_installed("tmbstan")
  skip_if(frmtmb.sample:::tmbstan_build_broken(),
          "this tmbstan samples a standard normal instead of the model")
}

# the default on ndt, its bound, and whether the posterior holds the
# ML fraction inside its central 98 percent
snd_check <- function(form, fam, d, top) {
  fit <- suppressWarnings(frm(form, family = fam, data = d))
  defs <- unclass(frmtmb.sample:::default_priors_for(fit))
  s <- Filter(function(z) identical(z$dpar, "ndt"), defs)[[1L]]
  expect_identical(s$dist$kind, "uniform")
  expect_identical(s$dist$upper, top)
  # the link holds ndt inside the uniform's interval, so no bound is
  # handed to the sampler
  expect_true(is.na(s$ub))
  ml <- frm_linpred(fit, dpar = "ndt", type = "response")[1L]
  ds <- allow_warnings(
    suppressMessages(frmtmb.sample::frm_sample(
      fit, chains = 2, iter = 400, seed = 4, refresh = 0)),
    c("Effective Samples Size", "The largest R-hat is",
      "divergent transitions", "Examine the pairs", "maximum treedepth"))
  v <- as.matrix(ds)[, "ndt"]
  q <- stats::quantile(v, c(0.01, 0.99))
  expect_gt(ml, q[[1L]])
  expect_lt(ml, q[[2L]])
}

test_that("a wiener ndt fraction under ndt_group() takes uniform(0, 1)", {
  snd_skip()
  set.seed(20261120)
  n <- 120
  d <- ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d$subj <- factor(rep(1:4, length.out = n))
  d$x <- stats::rnorm(n)
  snd_check(bf(rt | dec(upper) + ndt_group(subj) ~ x), wiener(), d, 1)
  # and without ndt_group() the uniform is brms's, up to the fastest
  # response, which is also the top of the link's image
  fit <- suppressWarnings(frm(bf(rt | dec(upper) ~ x), family = wiener(),
                              data = d))
  s <- Filter(function(z) identical(z$dpar, "ndt"),
              unclass(frmtmb.sample:::default_priors_for(fit)))[[1L]]
  expect_identical(s$dist$upper, min(d$rt))
  expect_true(is.na(s$ub))
})

# lba() and rdm() sample slowly (about 6 minutes for lba() here), so
# they run in the gated tier
snd_skip_slow <- function() {
  snd_skip()
  skip_if_not(identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true"),
              "slow sampler test: set FRMTMB_BRMS_FIT_TESTS=true")
}

test_that("an lba ndt fraction under ndt_group() takes uniform(0, 1)", {
  snd_skip_slow()
  # the chains of this lba() posterior do not mix in 2 x 400 (R-hat 2.0
  # to 2.9 on three data seeds, dev/sampleprior-ndtgroup.R), so where
  # the ML fraction falls in them is the chains' luck; on one of those
  # seeds it fell outside their central 98 percent. What the default
  # decides, and what is checked, is that they are not pinned below
  # min_Y, where the first lane build held them (max 0.29 against an ML
  # fraction of 0.55 on these data)
  set.seed(20261121)
  d <- lba_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.4, ndt = 0.2)
  d$subj <- factor(rep(1:4, length.out = nrow(d)))
  fit <- suppressWarnings(frm(bf(rt | vint(choice) + ndt_group(subj) ~ 1),
                              family = lba(2), data = d))
  s <- Filter(function(z) identical(z$dpar, "ndt"),
              unclass(frmtmb.sample:::default_priors_for(fit)))[[1L]]
  expect_identical(s$dist$upper, 1)
  expect_true(is.na(s$ub))
  ml <- frm_linpred(fit, dpar = "ndt", type = "response")[1L]
  expect_gt(ml, min(d$rt))
  ds <- allow_warnings(
    suppressMessages(frmtmb.sample::frm_sample(
      fit, chains = 2, iter = 400, seed = 4, refresh = 0)),
    c("Effective Samples Size", "The largest R-hat is",
      "divergent transitions", "Examine the pairs", "maximum treedepth"))
  v <- as.matrix(ds)[, "ndt"]
  expect_gt(stats::quantile(v, 0.5), min(d$rt))
})

test_that("an rdm ndt fraction at its ML edge still samples", {
  snd_skip_slow()
  # on these data rdm()'s ML ndt fraction is 1, the top of its support:
  # the support is closed, so chain 1 starts there, where a NaN at the
  # edge stopped the sampler ("no more scalars to read";
  # dev/sampleprior-rdm.R), and the draws reach up to it
  set.seed(20261121)
  d <- rdm_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.5, ndt = 0.2)
  d$subj <- factor(rep(1:4, length.out = nrow(d)))
  fit <- suppressWarnings(frm(bf(rt | vint(choice) + ndt_group(subj) ~ 1),
                              family = rdm(2), data = d))
  s <- Filter(function(z) identical(z$dpar, "ndt"),
              unclass(frmtmb.sample:::default_priors_for(fit)))[[1L]]
  expect_identical(s$dist$upper, 1)
  ml <- frm_linpred(fit, dpar = "ndt", type = "response")[1L]
  ds <- allow_warnings(
    suppressMessages(frmtmb.sample::frm_sample(
      fit, chains = 2, iter = 300, seed = 4, refresh = 0)),
    c("Effective Samples Size", "The largest R-hat is",
      "divergent transitions", "Examine the pairs", "maximum treedepth"))
  v <- as.matrix(ds)[, "ndt"]
  expect_gt(stats::quantile(v, 0.99), 0.99 * ml)
})
