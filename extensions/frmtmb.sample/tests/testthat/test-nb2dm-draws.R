# brmsfamily("negbinomial2") and dirichlet_multinomial() on posterior
# draws: brms 2.23.0's default priors, and every draws method a
# negbinomial() or multinomial() fit answers, with brms's shapes.
# Tolerances are multiples of Monte Carlo standard errors the run
# measures from its own draws.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

nbdm_rows <- function(tab) {
  tab <- as.data.frame(tab)
  stats::setNames(tab$prior, paste(tab$class, tab$coef, tab$group, tab$resp,
                                   tab$dpar, sep = "|"))
}

nbdm_draws_data <- function(seed, n = 200) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                   g = factor(rep(1:8, length.out = n)))
  dd$y <- stats::rnbinom(n, mu = exp(0.8 + 0.5 * dd$x), size = 2.5)
  E <- cbind(0, 0.3 + 0.6 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$N <- sample(4:20, n, TRUE)
  Y <- t(vapply(seq_len(n), function(i) {
    p <- stats::rgamma(3, P[i, ] * 6)
    as.numeric(stats::rmultinom(1, dd$N[i], p / sum(p)))
  }, numeric(3)))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}

nbdm_sample <- function(fit) {
  allow_warnings(
    ds <- suppressMessages(frm_sample(fit, chains = 1, iter = 400,
                                      warmup = 200, refresh = 0, seed = 9)),
    c("Effective Samples Size", "R-hat"))
  ds
}

test_that("the default priors are brms's default_prior()", {
  skip_if_not_installed("brms")
  dd <- nbdm_draws_data(1)
  cases <- list(
    list(frmtmb::bf(y ~ x), brms::bf(y ~ x), "negbinomial2"),
    list(frmtmb::bf(y ~ x, sigma ~ z), brms::bf(y ~ x, sigma ~ z),
         "negbinomial2"),
    list(frmtmb::bf(y ~ x + (1 | g)), brms::bf(y ~ x + (1 | g)),
         "negbinomial2"),
    list(frmtmb::bf(Y | trials(N) ~ x), brms::bf(Y | trials(N) ~ x),
         "dirichlet_multinomial"),
    list(frmtmb::bf(Y | trials(N) ~ x, phi ~ z),
         brms::bf(Y | trials(N) ~ x, phi ~ z), "dirichlet_multinomial"),
    list(frmtmb::bf(Y | trials(N) ~ x + (1 | g)),
         brms::bf(Y | trials(N) ~ x + (1 | g)), "dirichlet_multinomial"))
  n_dens <- 0L
  for (cs in cases) {
    ours <- nbdm_rows(default_prior(cs[[1]], dd,
                                    family = frmtmb::brmsfamily(cs[[3]]),
                                    route = "sample"))
    theirs <- nbdm_rows(brms::default_prior(cs[[2]], dd,
                                            family = brms::brmsfamily(cs[[3]])))
    # the rows where brms states a density; its "" is a flat or an
    # inherited slot
    dens <- theirs[nzchar(theirs)]
    expect_identical(ours[names(dens)], dens)
    n_dens <- n_dens + length(dens)
  }
  # the count the cases hold (dev/nb2dm-out/priors.txt): a case that
  # lost its rows would otherwise pass on nothing
  expect_identical(n_dens, 18L)
  # brms's own hook for negbinomial2: gamma, not negbinomial's inv_gamma
  pr <- nbdm_rows(default_prior(y ~ x, dd,
                                family = frmtmb::brmsfamily("negbinomial2"),
                                route = "sample"))
  expect_identical(pr[["sigma||||"]], "gamma(0.4, 0.3)")
})

test_that("negbinomial2 draws answer negbinomial's methods", {
  skip_sampler()
  dd <- nbdm_draws_data(2)
  fit <- frm(y ~ x, family = frmtmb::brmsfamily("negbinomial2"), data = dd)
  ds <- nbdm_sample(fit)
  expect_true(all(c("b_Intercept", "b_x", "sigma") %in%
                    posterior::variables(posterior::as_draws_df(ds))))
  S <- posterior::ndraws(posterior::as_draws_df(ds))
  n <- nrow(dd)
  expect_identical(dim(log_lik(ds)), c(S, n))
  ep <- posterior_epred(ds)
  pp <- posterior_predict(ds)
  expect_identical(dim(ep), c(S, n))
  expect_identical(dim(pp), c(S, n))
  # the draws' log-likelihood is the family's density at each draw
  dr <- posterior::as_draws_df(ds)
  ll1 <- stats::dnbinom(dd$y[1], size = 1 / dr$sigma,
                        mu = exp(dr$b_Intercept + dr$b_x * dd$x[1]),
                        log = TRUE)
  expect_lt(max(abs(log_lik(ds)[, 1] - ll1) / abs(ll1)), 1e-10)
  # predictive draws scatter about their own draw's mean
  se <- apply(pp - ep, 2, stats::sd) / sqrt(S)
  z <- stats::qnorm(1 - 1e-6 / 2)
  expect_true(all(abs(colMeans(pp) - colMeans(ep)) <= z * se))
  # sigma's draws are 1 / shape's: the predictive variance mu + sigma mu^2
  mu1 <- ep[, 1]
  v1 <- mu1 + dr$sigma * mu1^2
  set.seed(3)
  pv <- posterior_predict(ds)[, 1]
  d2 <- (pv - mu1)^2 - v1
  expect_lt(abs(mean(d2)), z * stats::sd(d2) / sqrt(S))
})

test_that("dirichlet_multinomial draws have brms's shapes and trials", {
  skip_sampler()
  dd <- nbdm_draws_data(3)
  fit <- frm(Y | trials(N) ~ x, family = frmtmb::dirichlet_multinomial(),
             data = dd)
  ds <- nbdm_sample(fit)
  vars <- posterior::variables(posterior::as_draws_df(ds))
  expect_true(all(c("b_mub_Intercept", "b_mub_x", "b_muc_Intercept",
                    "b_muc_x", "phi") %in% vars))
  S <- posterior::ndraws(posterior::as_draws_df(ds))
  n <- nrow(dd)
  # brms's posterior_epred(): S x N x K, named by the response's columns,
  # each draw's expected counts summing to the row's trials
  ep <- posterior_epred(ds)
  expect_identical(dim(ep), c(S, n, 3L))
  expect_identical(dimnames(ep)[[3L]], c("a", "b", "c"))
  tot <- apply(ep, 1:2, sum)
  expect_lt(max(abs(sweep(tot, 2, dd$N, "/") - 1)), 1e-12)
  pp <- posterior_predict(ds)
  expect_identical(dim(pp), c(S, n, 3L))
  expect_true(all(apply(pp, 1:2, sum) == rep(dd$N, each = S)))
  pl <- posterior_linpred(ds)
  expect_identical(dim(pl), c(S, n, 2L))
  expect_identical(dim(log_lik(ds)), c(S, n))
  # the draws' log-likelihood is the family's density at each draw
  dr <- posterior::as_draws_df(ds)
  e2 <- dr$b_mub_Intercept + dr$b_mub_x * dd$x[1]
  e3 <- dr$b_muc_Intercept + dr$b_muc_x * dd$x[1]
  y1 <- dd$Y[1, ]
  ref <- vapply(seq_len(S), function(s) {
    p <- exp(c(0, e2[s], e3[s]))
    a <- dr$phi[s] * p / sum(p)
    lgamma(dr$phi[s]) + lgamma(sum(y1) + 1) - lgamma(sum(y1) + dr$phi[s]) +
      sum(lgamma(y1 + a) - lgamma(a) - lgamma(y1 + 1))
  }, 0)
  expect_lt(max(abs(log_lik(ds)[, 1] - ref) / abs(ref)), 1e-10)
  # newdata: trials() from newdata, in the predictive draws too
  nd <- data.frame(x = c(-1, 1), N = c(3L, 200L))
  pn <- posterior_predict(ds, newdata = nd)
  expect_true(all(apply(pn[, 1, ], 1, sum) == 3L))
  expect_true(all(apply(pn[, 2, ], 1, sum) == 200L))
  en <- posterior_epred(ds, newdata = nd)
  se <- apply(pn[, 2, ] - en[, 2, ], 2, stats::sd) / sqrt(S)
  z <- stats::qnorm(1 - 1e-6 / 2)
  expect_true(all(abs(colMeans(pn[, 2, ]) - colMeans(en[, 2, ])) <= z * se))
  # the polytomous refusals, in the draws methods' words
  expect_error(predictive_error(ds), "'dirichlet_multinomial' family",
               fixed = TRUE)
  expect_error(pp_check(ds, type = "error_binned"),
               paste("the 'dirichlet_multinomial' response is a set of",
                     "counts over categories"), fixed = TRUE)
})

test_that("the count families here are core's multinomial_families", {
  # frmtmb.sample keeps its own copy because core's list is internal; a
  # family core adds there must be added here too
  expect_identical(frmtmb.sample:::draws_count_cats,
                   frmtmb:::multinomial_families)
})
