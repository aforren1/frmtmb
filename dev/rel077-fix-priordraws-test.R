# ---- the 0.77.0 consolidation: the four lanes together -----------------
# priordraws' review asked for these on the merge with reslot, sigmainn
# and mointerop (dev/priordraws-findings.md, section 13; the probe is
# dev/rel077-overlap.R)

pd_ac_data <- function() {
  set.seed(1)
  d <- data.frame(x = stats::rnorm(60), h = factor(rep(1:5, each = 12)),
                  t = rep(1:12, 5))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(60)
  d
}

test_that("default ar, ma and cosy get no prior draws, and \"only\" names them", {
  skip_sampler()
  # brms 2.23.0 declares ar and ma (cov = TRUE) and cosy on their
  # interval with no prior statement, so "yes" draws none and "only" is
  # refused naming the parameter (dev/review-priordraws-acdefault.R).
  # The priordraws lane alone drew them, from the 0.76.0 uniform default
  d <- pd_ac_data()
  pr <- set_prior("normal(0, 1)", class = "b")
  for (cs in list(list(f = y ~ x + ar(time = t, gr = h, cov = TRUE),
                       p = "ar"),
                  list(f = y ~ x + ma(time = t, gr = h, cov = TRUE),
                       p = "ma"),
                  list(f = y ~ x + cosy(time = t, gr = h), p = "cosy"))) {
    ds <- pd_quiet(frm_sample(bf(cs$f), data = d, family = gaussian(),
                              prior = pr, sample_prior = "yes", chains = 1,
                              iter = 300, seed = 1, refresh = 0))
    nm <- names(prior_draws(ds))
    expect_identical(nm, c("Intercept", "b", "sigma"), label = cs$p)
    expect_error(pd_quiet(frm_sample(bf(cs$f), data = d, family = gaussian(),
                                     prior = pr, sample_prior = "only",
                                     chains = 1, iter = 100, refresh = 0)),
                 paste0("Error occurred for parameter '", cs$p, "'"),
                 fixed = TRUE)
  }
})

test_that("\"only\" draws sigma under ar(cov = TRUE) from its own prior", {
  skip_sampler()
  # sigma is brms's innovation sd under cov = TRUE, and its prior is on
  # that parameter, so "only" gives exponential(1) on the sigma column.
  # With sigma the marginal sd (the priordraws lane alone) the mean was
  # 7.95 Monte Carlo standard errors from 1 (priordraws' review, O3)
  d <- pd_ac_data()
  pr <- set_prior("normal(0, 1)", class = "b") +
    set_prior("normal(1, 2)", class = "Intercept") +
    set_prior("exponential(1)", class = "sigma") +
    set_prior("normal(0, 0.5)", class = "ar")
  do <- pd_quiet(frm_sample(bf(y ~ x + ar(time = t, gr = h, cov = TRUE)),
                            data = d, family = gaussian(), prior = pr,
                            sample_prior = "only", chains = 2, iter = 2000,
                            seed = 9, refresh = 0))
  s <- posterior::subset_draws(posterior::as_draws_array(do), "sigma")
  expect_lt(abs(mean(s) - 1), 4 * posterior::mcse_mean(s))
  for (k in c(0.25, 0.5, 0.75)) {
    expect_lt(abs(stats::quantile(s, k) - stats::qexp(k)),
              4 * posterior::mcse_quantile(s, k), label = paste("q", k))
  }
})

test_that("a two-value mo() draws its one prior weight as 1", {
  skip_sampler()
  # brms's simplex[1]: mointerop's two-value mo() under "yes"
  set.seed(3)
  dm <- data.frame(m2 = rep(1:2, 30))
  dm$y <- 0.5 * dm$m2 + stats::rnorm(60)
  ds <- pd_quiet(frm_sample(bf(y ~ mo(m2)), data = dm, family = gaussian(),
                            prior = set_prior("normal(0, 1)", class = "b"),
                            sample_prior = "yes", chains = 1, iter = 300,
                            seed = 7, refresh = 0))
  pd <- prior_draws(ds)
  expect_identical(names(pd), c("Intercept", "simo_mom21[1]", "bsp",
                                "sigma"))
  expect_identical(unique(pd[["simo_mom21[1]"]]), 1)
})

test_that("under prior = \"flat\" nothing is drawn and nothing has a ratio", {
  skip_sampler()
  # reslot's prior = "flat" holds (flat) entries, which are no proper
  # prior: "yes" draws no column for them, as brms draws none for "" in
  # every slot, the evidence ratio is NA as brms's, and "only" is
  # refused with brms's words
  set.seed(2)
  dg <- data.frame(x = stats::rnorm(60), g = factor(rep(1:6, 10)))
  dg$y <- 1 + 0.5 * dg$x + stats::rnorm(6)[dg$g] + stats::rnorm(60)
  ff <- bf(y ~ x + (1 + x | g))
  ds <- pd_quiet(frm_sample(ff, data = dg, family = gaussian(),
                            prior = "flat", sample_prior = "yes",
                            chains = 1, iter = 300, seed = 5, refresh = 0))
  expect_identical(ncol(prior_draws(ds)), 0L)
  h <- pd_quiet(hypothesis(ds, "x = 0"))
  expect_true(is.na(h$hypothesis$Evid.Ratio))
  h2 <- pd_quiet(hypothesis(ds, "sd_g__Intercept = 0.5", class = ""))
  expect_true(is.na(h2$hypothesis$Evid.Ratio))
  expect_error(pd_quiet(frm_sample(ff, data = dg, family = gaussian(),
                                   prior = "flat", sample_prior = "only",
                                   chains = 1, iter = 100, refresh = 0)),
               "Error occurred for parameter 'b'", fixed = TRUE)
  # the guard absent: a proper prior on b gives its column and a ratio
  dp <- pd_quiet(frm_sample(ff, data = dg, family = gaussian(),
                            prior = set_prior("normal(0, 1)", class = "b"),
                            sample_prior = "yes", chains = 1, iter = 300,
                            seed = 5, refresh = 0))
  expect_true("b" %in% names(prior_draws(dp)))
  hp <- pd_quiet(hypothesis(dp, "x = 0"))
  expect_true(is.finite(hp$hypothesis$Evid.Ratio))
})
