# The latent residual against brms 2.23.0's own Stan program. Gated, and the
# programs come from the suite's compile cache. test-latres.R holds what needs
# no brms.

lrb_data <- function() {
  set.seed(167)
  d <- expand.grid(t = 1:5, g = factor(1:30))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  e <- as.vector(replicate(30, stats::arima.sim(list(ar = 0.5), 5, sd = 0.5)))
  eta <- 0.8 + 0.4 * d$x + e
  d$count <- stats::rpois(n, exp(eta))
  d$bin <- stats::rbinom(n, 1, stats::plogis(eta - 0.8))
  d$trials <- stats::rbinom(n, 10, stats::plogis(eta - 0.8))
  d$tr <- 10L
  d$pos <- stats::rgamma(n, 4, 4 / exp(eta))
  d$unit <- stats::rbeta(n, 8 * stats::plogis(eta - 0.8),
                         8 * (1 - stats::plogis(eta - 0.8)))
  d$ord <- factor(findInterval(eta + stats::rlogis(n), c(0, 1, 2)) + 1L,
                  ordered = TRUE)
  d
}

lrb_terms <- c("ar(t, g, p = 1, cov = TRUE)", "ar(t, g, p = 1)",
               "ar(t, g, p = 2)", "ma(t, g, q = 1, cov = TRUE)",
               "arma(t, g, p = 1, q = 1, cov = TRUE)", "cosy(t, g)",
               "unstr(t, g)")
lrb_fams <- list(
  list(brms::brmsfamily("poisson"), poisson(), "count"),
  list(brms::brmsfamily("negbinomial"), negbinomial(), "count"),
  list(brms::brmsfamily("bernoulli"), bernoulli(), "bin"),
  list(brms::brmsfamily("binomial"), binomial(), "trials | trials(tr)"),
  list(brms::brmsfamily("Gamma", "log"), Gamma(link = "log"), "pos"),
  list(brms::brmsfamily("lognormal"), lognormal(), "pos"),
  list(brms::brmsfamily("beta"), Beta(), "unit"),
  list(brms::brmsfamily("cumulative"), cumulative(), "ord"))

test_that("the joint at the modes is an identity with brms's log_prob", {
  skip_unless_brms_fit()
  # frm() refuses the post-link forms, and frm_sample() reads this objective
  local_mocked_bindings(latres_refuse_laplace = function(spec, frame) NULL,
                        latres_refuse_margin = function(fit) NULL,
                        .package = "frmtmb")
  d <- lrb_data()
  at_edge <- 0L
  for (tm in lrb_terms) {
    for (fm in lrb_fams) {
      f <- stats::as.formula(paste(fm[[3]], "~ x +", tm))
      fit <- suppressWarnings(frm(bf(f), family = fm[[2]], data = d))
      # a chain at the numeric edge of its range has no Cholesky factor in
      # brms's coordinates, so its program has no point to compare at
      r <- tryCatch(
        brms_lp_check(brms::bf(f), fm[[1]], d, fit, joint = TRUE),
        error = function(e) {
          if (!grepl("leading minor", conditionMessage(e), fixed = TRUE)) stop(e)
          NULL
        })
      if (is.null(r)) {
        at_edge <- at_edge + 1L
        next
      }
      expect_lt(abs(r$measured_const), 1e-9 * max(1, abs(r$ours)))
    }
  }
  expect_lte(at_edge, 2L)
})

test_that("frmtmb's table of built-in densities is brms's", {
  skip_unless_brms()
  d <- lrb_data()
  d$y <- d$count + 1L
  d$b <- d$bin
  # brms's program writes `mu = inv_link(mu);` before the ARMA loop exactly
  # where stan_eta_transform() says so, and that is what the rule stands for
  arms <- list(
    list("poisson", "log", "y"), list("poisson", "sqrt", "y"),
    list("poisson", "identity", "y"),
    list("negbinomial", "log", "y"), list("geometric", "log", "y"),
    list("bernoulli", "logit", "b"), list("bernoulli", "probit", "b"),
    list("binomial", "logit", "trials | trials(tr)"),
    list("binomial", "cloglog", "trials | trials(tr)"),
    list("zero_inflated_poisson", "log", "y"),
    list("zero_inflated_negbinomial", "log", "y"),
    list("hurdle_poisson", "log", "y"),
    list("hurdle_negbinomial", "log", "y"),
    list("Gamma", "log", "y"), list("Gamma", "inverse", "y"),
    list("lognormal", "identity", "y"), list("weibull", "log", "y"),
    list("exponential", "log", "y"), list("beta", "logit", "unit"),
    list("cumulative", "logit", "ord"), list("acat", "logit", "ord"))
  for (a in arms) {
    fb <- brms::brmsfamily(a[[1]], a[[2]])
    code <- brms::stancode(stats::as.formula(paste(a[[3]], "~ x + ar(t, g)")),
                           data = d, family = fb)
    l <- strsplit(code, "\n", fixed = TRUE)[[1]]
    arma <- grep("include ARMA terms", l, fixed = TRUE)
    inv <- grep("^ *mu = ", l)
    brms_says <- length(inv) > 0L && inv[1L] < arma[1L]
    fm <- as_frmtmb_family(eval(parse(text = sprintf(
      "%s(link = \"%s\")", if (a[[1]] == "beta") "Beta" else a[[1]],
      a[[2]]))))
    expect_identical(latres_on_response(fm, list()), brms_says,
                     label = paste(a[[1]], a[[2]]))
  }
  # no family has a built-in cdf on these, so a censored row is post-link
  fm <- as_frmtmb_family(poisson())
  expect_true(latres_on_response(fm, list(cens = 0L)))
  expect_false(latres_on_response(fm, list()))
})

test_that("the methods read brms's err at the same parameter values", {
  skip_unless_brms_fit()
  d <- lrb_data()
  f <- count ~ x + ar(t, g, cov = TRUE)
  fit <- suppressWarnings(frm(bf(f), family = poisson(), data = d))
  fb <- brms::bf(f)
  fam <- brms::brmsfamily("poisson")
  bfit <- suppressWarnings(brms_fixed_fit(fb, fam, d, fit))
  expect_exact_num(fitted(fit)[, 1L], brms::posterior_epred(bfit)[1L, ])
  expect_exact_num(fitted(fit, re_formula = NA)[, 1L],
                   brms::posterior_epred(bfit, re_formula = NA)[1L, ])
  expect_exact_num(frm_linpred(fit), brms::posterior_linpred(bfit)[1L, ])
  expect_exact_num(residuals(fit, type = "response")[, 1L],
                   residuals(bfit, method = "posterior_epred",
                             type = "ordinary", summary = FALSE)[1L, ])
  expect_exact_num(frm_row_loglik(fit), brms::log_lik(bfit)[1L, ])
  # a grid drops the term on both sides
  ea <- as.data.frame(emmeans::emmeans(fit, ~ x))
  eb <- as.data.frame(emmeans::emmeans(bfit, ~ x))
  expect_exact_num(ea$rate %||% ea$emmean %||% ea[[2L]], eb[[2L]])
  # the names of its parameters are brms's
  expect_true(all(c("b_Intercept", "b_x", "ar[1]", "sderr") %in%
                    brms::variables(bfit)))
  expect_true(all(c("b_Intercept", "b_x", "ar[1]", "sderr") %in%
                    variables(fit)))
})

test_that("fitted() on new rows is the mean of brms's posterior_epred()", {
  skip_unless_brms_fit()
  d <- lrb_data()
  nd <- data.frame(t = 1:5, g = factor(rep("new", 5L)), x = 0)
  arms <- list(
    list(count ~ x + ar(t, g, cov = TRUE), poisson(),
         brms::brmsfamily("poisson")),
    list(bin ~ x + ar(t, g, cov = TRUE), bernoulli(),
         brms::brmsfamily("bernoulli")))
  for (a in arms) {
    fit <- suppressWarnings(frm(bf(a[[1L]]), family = a[[2L]], data = d))
    bfit <- suppressWarnings(
      brms_fixed_fit(brms::bf(a[[1L]]), a[[3L]], d, fit, ndraws = 4000L))
    set.seed(167L)
    ep <- brms::posterior_epred(bfit, newdata = nd, allow_new_levels = TRUE)
    # bernoulli at this seed sits at a limit, where a row has no standard error
    mine <- allow_warnings(
      as.numeric(fitted(fit, newdata = nd, allow_new_levels = TRUE)[, 1L]),
      "predictions move along a direction")
    # against the Monte Carlo error of brms's own draws, one latent per draw
    se <- apply(ep, 2L, stats::sd) / sqrt(nrow(ep))
    expect_lt(max(abs(colMeans(ep) - mine) / se), 4.5)
  }
})

test_that("a nonlinear dpar beside a linear mu is brms's model", {
  skip_unless_brms_fit()
  d <- lrb_data()
  f <- count ~ x + ar(t, g, cov = TRUE)
  fit <- allow_boundary(suppressWarnings(
    frm(bf(f) + nlf(shape ~ exp(a)) + lf(a ~ 1), family = negbinomial(),
        data = d)))
  r <- brms_lp_check(brms::bf(f) + brms::nlf(shape ~ exp(a)) + brms::lf(a ~ 1),
                     brms::brmsfamily("negbinomial"), d, fit, joint = TRUE)
  expect_lt(abs(r$measured_const), 1e-9 * max(1, abs(r$ours)))
})

test_that("brms's cov = FALSE methods read a mean its model does not have", {
  skip_unless_brms_fit()
  # brms-45: predictor_ac() adds err and stops, so posterior_epred() leaves out
  # the AR part mu[n] += Err[n, 1:Kar] * ar of the Stan program it summarizes.
  # log_lik() summed over rows then misses log_prob. frmtmb's readers follow the
  # model; this pins both halves, and it fails the day brms changes.
  d <- lrb_data()
  f <- count ~ x + ar(t, g, p = 1)
  fit <- suppressWarnings(frm(bf(f), family = poisson(), data = d))
  fb <- brms::bf(f)
  fam <- brms::brmsfamily("poisson")
  bfit <- suppressWarnings(brms_fixed_fit(fb, fam, d, fit))
  bk <- frame_latres(fit$frame)[[1L]]
  th <- fit$estimates$theta[bk$theta_idx]
  z <- fit$estimates$b[bk$b_idx]
  lp <- suppressWarnings(brms_lp_check(fb, fam, d, fit, joint = TRUE))
  latent_part <- sum(stats::dnorm(z, log = TRUE))
  frm_rows <- sum(frm_row_loglik(fit))
  brms_rows <- sum(brms::log_lik(bfit)[1L, ])
  expect_lt(abs(frm_rows + latent_part - lp$lp), 1e-8 * abs(lp$lp))
  expect_gt(abs(brms_rows + latent_part - lp$lp), 1e-3)
  # brms's err alone
  eta <- as.numeric(cbind(1, d$x) %*% fit$estimates$beta) + exp(th[1L]) * z
  expect_exact_num(exp(eta), brms::posterior_epred(bfit)[1L, ])
})

test_that("brms's prior table is frm()'s", {
  skip_unless_brms()
  d <- lrb_data()
  rows <- function(p) {
    p <- as.data.frame(p)
    p <- p[p$class %in% c("ar", "ma", "cosy", "cortime", "sderr"), ]
    p <- p[order(p$class), c("class", "lb", "ub")]
    p$lb[is.na(p$lb)] <- ""
    p$ub[is.na(p$ub)] <- ""
    rownames(p) <- NULL
    p
  }
  for (tm in lrb_terms) {
    f <- stats::as.formula(paste("count ~ x +", tm))
    a <- rows(default_prior(bf(f), family = poisson(), data = d))
    b <- rows(brms::default_prior(brms::bf(f), data = d,
                                  family = brms::brmsfamily("poisson")))
    expect_equal(a, b, label = tm)
  }
})
