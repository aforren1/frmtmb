# brms's sample_prior on frm_sample(): prior draws ("yes"), the prior
# of the whole model ("only"), prior_draws(), and update(sample_prior =).
#
# The expected names are brms 2.23.0's, read from brms fits of the same
# models (dev/priordraws-brms-names.R; dev/priordraws-findings.md,
# section 1). Each draw column is checked against the exact prior it
# should come from, truncation included, by a one-sample KS test on a
# fixed seed. The agreement with brms's own draws, "only" runs and
# evidence ratios is measured in dev/priordraws-evidence.R.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

pd_data <- function() {
  set.seed(1)
  d <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60),
                  g = factor(rep(1:6, 10)), m = sample(1:4, 60, TRUE))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(60)
  d
}

pd_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

pd_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      d <- pd_data()
      pr <- set_prior("normal(0, 1)", class = "b") +
        set_prior("normal(0, 2)", class = "b", coef = "z") +
        set_prior("normal(1, 3)", class = "Intercept") +
        set_prior("exponential(2)", class = "sigma") +
        set_prior("lkj(2)", class = "cor") +
        set_prior("cauchy(0, 1)", class = "sd", group = "g", coef = "x")
      f <- bf(y ~ x + z + (1 + x | g))
      yes <- pd_quiet(frm_sample(f, data = d, family = gaussian(),
                                 prior = pr, sample_prior = "yes",
                                 chains = 2, iter = 1000, seed = 11,
                                 refresh = 0))
      no <- pd_quiet(frm_sample(f, data = d, family = gaussian(),
                                prior = pr, chains = 2, iter = 1000,
                                seed = 11, refresh = 0))
      cache <<- list(d = d, pr = pr, f = f, yes = yes, no = no)
    }
    cache
  }
})

test_that("sample_prior takes brms's values and refuses others", {
  csp <- frmtmb.sample:::check_sample_prior
  expect_identical(csp(NULL), "no")
  expect_identical(csp(TRUE), "yes")
  expect_identical(csp(FALSE), "no")
  expect_identical(csp("only"), "only")
  # match.arg() in brms, so a unique abbreviation is accepted
  expect_identical(csp("y"), "yes")
  expect_error(csp("prior"), "\"no\", \"yes\" or \"only\"")
  expect_error(csp(NA), "\"no\", \"yes\" or \"only\"")
})

test_that("the prior draws carry brms's names, in brms's order", {
  cs <- pd_case()
  # brms 2.23.0 on this model and prior: a coefficient-specific prior
  # on z makes one statement per coefficient of b (b_x, b_z), and the
  # coefficient-specific sd prior one per coefficient of sd_g
  expect_identical(names(prior_draws(cs$yes)),
                   c("Intercept", "b_x", "b_z", "sigma", "sd_g__Intercept",
                     "sd_g__x", "cor_g"))
  expect_identical(nrow(prior_draws(cs$yes)), ndraws(cs$yes))
  # the columns sit in the draws, under prior_, as in a brmsfit
  expect_true(all(paste0("prior_", names(prior_draws(cs$yes))) %in%
                    variables(cs$yes)))
  # and are not parameters: summary() leaves them out
  expect_false(any(grepl("^prior_", rownames(summary(cs$yes)))))
})

test_that("each prior column is a draw from its exact prior", {
  cs <- pd_case()
  pd <- prior_draws(cs$yes)
  ks <- function(v, cdf) stats::ks.test(v, cdf)$p.value
  # half-distributions where brms declares a lower bound of 0
  half_t <- function(q) 2 * stats::pt(q / 2.5, 3) - 1
  expect_gt(ks(pd$Intercept, function(q) stats::pnorm(q, 1, 3)), 0.001)
  expect_gt(ks(pd$b_x, stats::pnorm), 0.001)
  expect_gt(ks(pd$b_z, function(q) stats::pnorm(q, 0, 2)), 0.001)
  expect_gt(ks(pd$sigma, function(q) stats::pexp(q, 2)), 0.001)
  expect_gt(ks(pd$sd_g__Intercept, function(q) pmax(half_t(q), 0)), 0.001)
  expect_gt(ks(pd$sd_g__x, function(q) pmax(2 * stats::pcauchy(q) - 1, 0)),
            0.001)
  # one element of an LKJ(2) correlation matrix of dimension 2:
  # 2 Beta(2, 2) - 1
  expect_gt(ks(pd$cor_g, function(q) stats::pbeta((q + 1) / 2, 2, 2)),
            0.001)
  expect_true(all(pd$sigma > 0) && all(pd$sd_g__x > 0))
  # the construction's control: the half-t check rejects the full t
  set.seed(2)
  expect_lt(ks(2.5 * stats::rt(nrow(pd), 3), function(q) {
    pmax(half_t(q), 0)
  }), 0.001)
})

test_that("a bounded class prior is drawn truncated, as brms rejects", {
  skip_sampler()
  d <- pd_data()
  ds <- pd_quiet(frm_sample(
    bf(y ~ x + z), data = d, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "b", lb = 0, ub = 1),
    sample_prior = "yes", chains = 1, iter = 2000, seed = 3, refresh = 0))
  b <- prior_draws(ds)$b
  expect_true(all(b >= 0 & b <= 1))
  tn <- function(q) {
    (stats::pnorm(pmin(pmax(q, 0), 1)) - 0.5) / (stats::pnorm(1) - 0.5)
  }
  expect_gt(stats::ks.test(b, tn)$p.value, 0.001)
})

test_that("a mo() simplex and nu get brms's draws", {
  skip_sampler()
  d <- pd_data()
  ds <- pd_quiet(frm_sample(
    bf(y ~ mo(m) + x), data = d, family = student(),
    prior = set_prior("normal(0, 1)", class = "b"), sample_prior = "yes",
    chains = 1, iter = 2000, seed = 4, refresh = 0))
  pd <- prior_draws(ds)
  # brms: Intercept, b, simo_mom1[1..3], bsp, sigma, nu
  expect_identical(names(pd),
                   c("Intercept", "b", paste0("simo_mom1[", 1:3, "]"),
                     "bsp", "sigma", "nu"))
  w <- as.matrix(pd[, grep("^simo", names(pd))])
  expect_equal(rowSums(w), rep(1, nrow(w)))
  # each weight of a Dirichlet(1, 1, 1) is Beta(1, 2)
  expect_gt(stats::ks.test(w[, 1], function(q) stats::pbeta(q, 1, 2))$p.value,
            0.001)
  # nu has brms's lower bound of 1: gamma(2, 0.1) truncated there
  g1 <- stats::pgamma(1, 2, 0.1)
  expect_true(all(pd$nu >= 1))
  expect_gt(stats::ks.test(pd$nu, function(q) {
    (stats::pgamma(pmax(q, 1), 2, 0.1) - g1) / (1 - g1)
  })$p.value, 0.001)
})

test_that("prior_draws(variable = ) reads brms's longest match", {
  cs <- pd_case()
  pd <- prior_draws(cs$yes)
  one <- prior_draws(cs$yes, "b_x")
  # the coefficient's own column, in random order, with brms's rows
  expect_identical(dimnames(one),
                   list(as.character(seq_len(ndraws(cs$yes))), "b_x"))
  expect_identical(sort(one$b_x), sort(pd$b_x))
  # sd_g__x reads sd_g__x, not sd_g; a name nothing starts is left out
  two <- prior_draws(cs$yes, c("sd_g__x", "sigma", "nothing"))
  expect_identical(names(two), c("sd_g__x", "sigma"))
  expect_identical(sort(two$sd_g__x), sort(pd$sd_g__x))
  # brms's deprecated alias, with brms's warning
  expect_warning(pp <- prior_draws(cs$yes, pars = "b_x"),
                 "Argument 'pars' is deprecated")
  expect_identical(names(pp), "b_x")
})

test_that("the overall intercept has no prior draws, as in brms", {
  cs <- pd_case()
  # here b is drawn per coefficient, so no column starts b_Intercept and
  # brms leaves it out before its ignore_prior() check, with no warning;
  # test-evidence-ratio.R has the class-wide b, where brms warns
  expect_no_warning(out <- prior_draws(cs$yes, c("b_Intercept", "b_x")))
  expect_identical(names(out), "b_x")
})

test_that("no prior draws without sample_prior = \"yes\"", {
  cs <- pd_case()
  expect_identical(prior_draws(cs$no), data.frame(NULL))
  expect_false(any(grepl("^prior_", variables(cs$no))))
  # the prior draws leave the posterior draws as they were
  expect_identical(cs$yes$draws[, colnames(cs$no$draws)], cs$no$draws)
})

test_that("the seed fixes the prior draws", {
  cs <- pd_case()
  again <- pd_quiet(frm_sample(cs$f, data = cs$d, family = gaussian(),
                               prior = cs$pr, sample_prior = "yes",
                               chains = 2, iter = 1000, seed = 11,
                               refresh = 0))
  expect_identical(prior_draws(again), prior_draws(cs$yes))
})

test_that("update(sample_prior = ) adds the prior draws, as in brms", {
  cs <- pd_case()
  up <- pd_quiet(update(cs$no, sample_prior = "yes"))
  expect_identical(up$draws, cs$yes$draws)
  # what update() cannot change is refused by name
  expect_error(update(cs$no, newdata = cs$d), "new model")
  expect_error(update(cs$no, family = poisson()), "changes the model")
  expect_error(update(cs$no, recompile = TRUE), "compiles no Stan")
})

test_that("sample_prior = \"only\" refuses a flat prior with brms's words", {
  skip_sampler()
  d <- pd_data()
  # brms 2.23.0: "Error occurred for parameter 'b'." under its defaults
  expect_error(frm_sample(bf(y ~ x + (1 | g)), data = d,
                          family = gaussian(), sample_prior = "only",
                          chains = 1, iter = 100, refresh = 0),
               "Error occurred for parameter 'b'", fixed = TRUE)
  # brms samples a coefficient left flat beside one with a prior; this
  # refuses it
  expect_error(
    frm_sample(bf(y ~ x + z), data = d, family = gaussian(),
               prior = set_prior("normal(0, 1)", class = "b", coef = "z"),
               sample_prior = "only", chains = 1, iter = 100, refresh = 0),
    "no proper priors")
})

test_that("sample_prior = \"only\" samples the prior of the whole model", {
  skip_sampler()
  d <- pd_data()
  do <- pd_quiet(frm_sample(
    bf(y ~ x + (1 | g)), data = d, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "b"), sample_prior = "only",
    chains = 2, iter = 2000, seed = 5, refresh = 0))
  # no prior draws of its own, as in brms
  expect_identical(prior_draws(do), data.frame(NULL))
  bx <- as.matrix(do)[, "b_x"]
  # b_x is its normal(0, 1) prior whatever the data say (the data's
  # slope is 0.5): the mean and the sd each within 4 of their Monte
  # Carlo errors of 0 and 1
  bm <- matrix(bx, ncol = 2)
  expect_lt(abs(mean(bx)) / posterior::mcse_mean(bm), 4)
  expect_lt(abs(stats::sd(bx) - 1) / posterior::mcse_sd(bm), 4)
  # the prior predictive distribution has the draws' shape
  yrep <- posterior_predict(do)
  expect_identical(dim(yrep), c(ndraws(do), nrow(d)))
})

test_that("the \"only\" refusals do not fire on a correct model", {
  skip_sampler()
  d <- pd_data()
  # every class this model has, with proper priors: correlated group
  # terms, a predicted sigma, a mo() simplex and nu
  pr <- set_prior("normal(0, 1)", class = "b") +
    set_prior("normal(0, 1)", class = "b", dpar = "sigma")
  expect_no_error(do <- pd_quiet(frm_sample(
    bf(y ~ x + mo(m) + (1 + x | g), sigma ~ z), data = d,
    family = student(), prior = pr, sample_prior = "only", chains = 1,
    iter = 300, seed = 6, refresh = 0)))
  expect_identical(ndraws(do), 150L)
})

test_that("missing values of an mi() term are refused under \"only\"", {
  skip_sampler()
  d <- pd_data()
  d$xm <- d$x
  d$xm[1:5] <- NA
  expect_error(
    frm_sample(bf(y ~ mi(xm)) + bf(xm | mi() ~ 1) + set_rescor(FALSE),
               data = d, family = gaussian(),
               prior = set_prior("normal(0, 1)", class = "b", resp = "y"),
               sample_prior = "only", chains = 1, iter = 100, refresh = 0),
    "5 missing values of xm have no density", fixed = TRUE)
})

test_that("a mi(sdy) response with no missing value samples under \"only\"", {
  skip_sampler()
  # brms 2.23.0 keeps normal_lpdf(Y[Jme] | Yl[Jme], noise[Jme]) outside
  # `if (!prior_only)` (dev/review-priordraws-mi.R), so the latent
  # responses are drawn around the observed ones with sd sdy
  d <- pd_data()
  d$sdy <- 0.3
  do <- pd_quiet(frm_sample(
    bf(y | mi(sdy) ~ x), data = d, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "b"), sample_prior = "only",
    chains = 2, iter = 1000, seed = 9, refresh = 0))
  bx <- matrix(as.matrix(do)[, "b_x"], ncol = 2)
  expect_lt(abs(mean(bx)) / posterior::mcse_mean(bx), 4)
  # the latent value of row 1 is normal(y[1], 0.3) under its measurement
  j <- frmtmb.sample:::draws_par_index(do$fit)[["miss"]][1L]
  yl <- matrix(do$draws[, j], ncol = 2)
  expect_lt(abs(mean(yl) - d$y[1]) / posterior::mcse_mean(yl), 4)
  expect_lt(abs(stats::sd(yl) - 0.3) / posterior::mcse_sd(yl), 4)
})

test_that("an me() model is refused under \"only\", naming brms's classes", {
  skip_sampler()
  d <- pd_data()
  d$sdx <- 0.3
  # brms 2.23.0: "Error occurred for parameter 'meanme_1'" until meanme
  # and sdme have proper priors (dev/priordraws-brms-me.R); frmtmb has
  # no slot for them
  expect_error(
    frm_sample(bf(y ~ me(x, sdx)), data = d, family = gaussian(),
               prior = set_prior("normal(0, 1)", class = "b"),
               sample_prior = "only", chains = 1, iter = 100, refresh = 0),
    "parameter 'meanme_1'. frmtmb has no prior slot yet for brms's classes",
    fixed = TRUE)
})

test_that("ordinal threshold priors of their own are drawn as brms's", {
  skip_sampler()
  set.seed(1)
  d <- data.frame(xo = stats::rnorm(60), sdx = 0.3, z = stats::rnorm(60))
  d$x <- d$xo + stats::rnorm(60, 0, 0.3)
  d$y <- 1 + 0.5 * d$xo + stats::rnorm(60)
  d$yo <- cut(d$y, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
  # brms 2.23.0 on these data (dev/priordraws-brms-me.R): b,
  # Intercept__1 from normal(0, 2), Intercept__2 and __3 from the class
  # default student_t(3, 0, 2.5), in that order
  ds <- pd_quiet(frm_sample(
    bf(yo ~ z), data = d, family = cumulative(),
    prior = set_prior("normal(0, 2)", class = "Intercept", coef = "1") +
      set_prior("normal(0, 1)", class = "b"),
    sample_prior = "yes", chains = 1, iter = 2000, seed = 10, refresh = 0))
  pd <- prior_draws(ds)
  expect_identical(names(pd), c("b", "Intercept__1", "Intercept__2",
                                "Intercept__3"))
  expect_gt(stats::ks.test(pd$Intercept__1, function(q) {
    stats::pnorm(q, 0, 2)
  })$p.value, 0.001)
  expect_gt(stats::ks.test(pd$Intercept__2, function(q) {
    stats::pt(q / 2.5, 3)
  })$p.value, 0.001)
})

test_that("update() names an unnamed argument that is not a formula", {
  cs <- pd_case()
  expect_error(update(cs$no, 5), "the first unnamed argument is the new ",
               fixed = TRUE)
})

test_that("a smooth in two predictors gets one sds column for each", {
  skip_sampler()
  d <- pd_data()
  # brms 2.23.0: Intercept sds_sx sds_sz Intercept_sigma sds_sigma_sx
  # (dev/priordraws-names2.R). s(x) in the location and s(x) in sigma
  # share a label and are two terms
  ds <- pd_quiet(frm_sample(bf(y ~ s(x) + s(z), sigma ~ s(x)), data = d,
                            family = gaussian(), sample_prior = "yes",
                            chains = 1, iter = 200, seed = 7, refresh = 0))
  expect_identical(names(prior_draws(ds)),
                   c("Intercept", "sds_sx", "sds_sz", "Intercept_sigma",
                     "sds_sigma_sx"))
})

test_that("a simplex of one weight is drawn as 1, as brms draws it", {
  # a two-value mo() has simplex[1] in brms, and dirichlet_rng() of
  # length one is 1 in every draw
  st <- frmtmb.sample:::prior_stmt("simo_mox1", "simo_mox1", "dirichlet(1)",
                                   type = "simplex", dim = 1L)
  m <- frmtmb.sample:::prior_draws_make(list(st), 5L, NULL)
  expect_identical(colnames(m), "prior_simo_mox1[1]")
  expect_identical(unname(m[, 1]), rep(1, 5))
})

test_that("a flat sd has no prior draws, and \"only\" refuses it as brms", {
  skip_sampler()
  d <- pd_data()
  pr <- set_prior("normal(0, 1)", class = "b") +
    set_prior("", class = "sd", lb = 0)
  # brms 2.23.0 on (1 | g) with a flat sd: "Error occurred for parameter
  # 'sd_1'" (dev/priordraws-brms-only.R)
  expect_error(frm_sample(bf(y ~ x + (1 | g)), data = d,
                          family = gaussian(), prior = pr,
                          sample_prior = "only", chains = 1, iter = 100,
                          refresh = 0),
               "Error occurred for parameter 'sd_1'", fixed = TRUE)
  ds <- pd_quiet(frm_sample(bf(y ~ x + (1 | g)), data = d,
                            family = gaussian(), prior = pr,
                            sample_prior = "yes", chains = 1, iter = 200,
                            seed = 8, refresh = 0))
  expect_identical(names(prior_draws(ds)), c("Intercept", "b", "sigma"))
})

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
