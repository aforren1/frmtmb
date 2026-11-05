# brms's class "sds" on the sampling route: the default, the table, the
# draws of the prior, and a smooth that class "sd" does not reach
# (dev/refuse-brms-sds.R, dev/refuse-table-compare.R).

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

sds_data <- function(n = 80) {
  set.seed(43)
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  f = factor(rep(c("a", "b"), length.out = n)),
                  g = factor(rep(1:8, length.out = n)))
  d$y <- sin(2 * pi * d$x) + d$z + stats::rnorm(n, 0, 0.3)
  d
}

test_that("a smooth's default is brms's, on class sds", {
  skip_if_not_installed("brms")
  dd <- sds_data()
  cases <- list(
    list(bf(y ~ s(x)), brms::bf(y ~ s(x))),
    list(bf(y ~ t2(x, z)), brms::bf(y ~ t2(x, z))),
    list(bf(y ~ f + s(x, by = f)), brms::bf(y ~ f + s(x, by = f))),
    list(bf(y ~ s(x) + (1 | g)), brms::bf(y ~ s(x) + (1 | g))),
    list(bf(y ~ s(x), sigma ~ s(z)), brms::bf(y ~ s(x), sigma ~ s(z)))
  )
  for (cc in cases) {
    ours <- as.data.frame(default_prior(cc[[1L]], data = dd, route = "sample"))
    theirs <- as.data.frame(brms::default_prior(cc[[2L]], data = dd))
    o <- ours[ours$class == "sds", c("prior", "coef", "dpar")]
    t <- theirs[theirs$class == "sds", c("prior", "coef", "dpar")]
    o$prior[o$prior == "(flat)"] <- ""
    rownames(o) <- rownames(t) <- NULL
    expect_equal(as.data.frame(lapply(o, as.character)),
                 as.data.frame(lapply(t, as.character)))
    # a class row carries the density and each term row inherits it
    expect_true(any(nzchar(o$prior) & !nzchar(o$coef)))
    expect_true(all(!nzchar(o$prior[nzchar(o$coef)])))
    # beside a random effect, sd stays sd
    if ("g" %in% ours$group) expect_true("sd" %in% ours$class)
  }
})

test_that("frm_sample() samples a smooth under sds, and draws its prior", {
  skip_sampler()
  dd <- sds_data()
  s1 <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ s(x)), data = dd, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "sds"),
    chains = 1, iter = 400, seed = 5, refresh = 0, sample_prior = "yes")))
  ps <- as.data.frame(prior_summary(s1))
  sr <- ps[ps$class == "sds", , drop = FALSE]
  expect_true(any(sr$prior == "normal(0, 1)" & !nzchar(sr$coef)))
  # brms's name for the prior draw of a smooth's sd, and its density:
  # half normal(0, 1), whose mean is sqrt(2 / pi)
  pd <- prior_draws(s1)
  expect_true("sds_sx" %in% colnames(pd))
  v <- pd[["sds_sx"]]
  expect_true(all(v > 0))
  se <- stats::sd(v) / sqrt(length(v))
  expect_lt(abs(mean(v) - sqrt(2 / pi)), 5 * se)
})

test_that("class sd is refused on a smooth, as brms refuses it", {
  skip_sampler()
  dd <- sds_data()
  expect_error(suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ s(x)), data = dd, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "sd"),
    chains = 1, iter = 20, seed = 5, refresh = 0))),
    "class \"sds\"", fixed = TRUE)
  # brms refuses a bound with a coef, on this class as on the others
  expect_error(suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ s(x)), data = dd, family = gaussian(),
    prior = set_prior("normal(0, 1)", class = "sds", coef = "s(x)", lb = 0.1),
    chains = 1, iter = 20, seed = 5, refresh = 0))),
    "may not be specified when using boundaries", fixed = TRUE)
})

test_that("mi() the Laplace route refuses is sampled, on its own chart", {
  skip_sampler()
  set.seed(177)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  y <- stats::rgamma(n, 3, 3 / exp(0.5 + 0.4 * d$x))
  d$resp <- d$z + 0.5 * y + stats::rnorm(n, 0, 0.5)
  d$y <- y
  d$y[sample(n, 6)] <- NA
  fit <- suppressWarnings(suppressMessages(frm_sample(
    mvbf(bf(resp ~ mi(y) + z),
         bf(y | mi() ~ x, family = stats::Gamma(link = "log")), rescor = FALSE),
    data = d, chains = 1, iter = 300, seed = 3, refresh = 0)))
  dm <- as.matrix(as_draws_matrix(fit))
  rows <- fit$fit$frame$mi_map$y$rows
  ep <- posterior_epred(fit, resp = "resp")
  # brms's name, on the response's own scale, which the predictor reads
  expect_false(any(startsWith(colnames(dm), "miss_")))
  y1 <- dm[, "Ymi_y[1]"]
  expect_true(all(y1 > 0))
  i <- rows[1L]
  hand <- dm[, "b_resp_Intercept"] + dm[, "b_resp_z"] * d$z[i] +
    dm[, "bsp_resp_miy"] * y1
  expect_lt(max(abs(ep[, i] - hand)) / max(abs(hand)), 1e-8)

  d$w <- d$x + stats::rexp(n) - stats::rexp(n)
  d$w[1:5] <- NA
  al <- suppressWarnings(suppressMessages(frm_sample(
    bf(w | mi() ~ x), family = asym_laplace(), data = d,
    chains = 1, iter = 200, seed = 4, refresh = 0)))
  expect_true(all(c("Ymi[1]", "Ymi[5]") %in%
                    colnames(as.matrix(as_draws_matrix(al)))))
  # frm() has no value to integrate when no density reads it
  ml <- suppressWarnings(frm(bf(w | mi() ~ x), family = asym_laplace(),
                             data = d))
  expect_null(ml$estimates$miss)
  # and frm_sample() of that fit holds the values a sampler reports
  fr <- suppressWarnings(suppressMessages(frm_sample(
    ml, chains = 1, iter = 200, seed = 4, refresh = 0)))
  expect_true(all(c("Ymi[1]", "Ymi[5]") %in%
                    colnames(as.matrix(as_draws_matrix(fr)))))
})
