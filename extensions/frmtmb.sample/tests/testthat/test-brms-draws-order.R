# The draws in brms's order and under brms's names: variables(), the
# as_draws_*() objects and as.matrix() list brms's reorder_pars() order,
# a smooth's basis weights are s_<label>[k], an exact gp()'s z values
# zgp_<term>[k], ar(cov = FALSE)'s states zerr[k], and a one-response
# draws matrix carries dimnames list(NULL, NULL). The expected orders are
# brms 2.23.0's on the same shapes (dev/names-brms.R,
# dev/names-out/brms-*.txt; dev/names-findings.md).

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

ord_data <- function(n = 60, seed = 7) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(sprintf("g%d", 1:6), each = n / 6)),
                  t = rep(seq_len(n / 6), 6))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, exp(0.2 * d$z))
  d$cnt <- stats::rpois(n, exp(0.5 + 0.3 * d$x))
  d
}

ord_sample <- function(f, family, d, seed = 3) {
  skip_sampler()
  suppressWarnings(suppressMessages(frm_sample(
    f, family = family, data = d, chains = 1, iter = 120, warmup = 60,
    refresh = 0, seed = seed)))
}

test_that("variables() and the draws objects list brms's order", {
  d <- ord_data()
  ds <- ord_sample(bf(y ~ x + (1 + x | g), sigma ~ z), gaussian(), d)
  v <- variables(ds)
  # a class b intercept first, then the other coefficients, as brms
  expect_identical(v[1:4], c("b_Intercept", "b_sigma_Intercept", "b_x",
                             "b_sigma_z"))
  # every level of one coefficient before the next: r_1_1, then r_1_2
  r <- v[startsWith(v, "r_")]
  expect_identical(r, c(paste0("r_g[g", 1:6, ",Intercept]"),
                        paste0("r_g[g", 1:6, ",x]")))
  expect_lt(max(which(startsWith(v, "sd_"))), min(which(startsWith(v, "r_"))))
  expect_identical(utils::tail(v, 1L), "lp__")
  # one order everywhere a user reads the variables
  expect_identical(setdiff(names(posterior::as_draws_df(ds)),
                           c(".chain", ".iteration", ".draw")), v)
  expect_identical(colnames(as.matrix(ds)), v)
  expect_identical(dimnames(as.array(ds))[[3L]], v)
  # the same columns as the stored matrix, which keeps the template order
  # its positional readers rely on
  expect_setequal(colnames(ds$draws), v)
  m <- as.matrix(ds)
  for (nm in c("b_x", "r_g[g3,x]", "sd_g__x")) {
    expect_identical(unname(m[, nm]), unname(ds$draws[, nm]))
  }
  # brms's deprecated `pars =` regex selects in variables() order too,
  # as brms's extract_pars(pars, variables(x)) does
  want <- v[grepl("^b_", v)]
  dep <- "Argument 'pars' is deprecated"
  pm <- allow_warnings(as.matrix(ds, pars = "^b_"), dep, require = dep)
  expect_identical(colnames(pm), want)
  pa <- allow_warnings(as.array(ds, pars = "^b_"), dep, require = dep)
  expect_identical(dimnames(pa)[[3L]], want)
  pd <- allow_warnings(as.data.frame(ds, pars = "^b_"), dep, require = dep)
  expect_identical(names(pd), want)
  pi <- allow_warnings(posterior_interval(ds, pars = "^b_"), dep)
  expect_identical(rownames(pi), want)
  skip_if_not_installed("coda")
  expect_identical(colnames(as.mcmc(ds, pars = "^b_")[[1L]]), want)
})

test_that("a smooth's draws are brms's s_, an exact gp()'s zgp_", {
  d <- ord_data()
  d$ys <- sin(d$x) + stats::rnorm(nrow(d), 0, 0.3)
  ds <- ord_sample(bf(ys ~ s(x)), gaussian(), d)
  v <- variables(ds)
  k <- length(ds$fit$frame$par_template$b)
  expect_identical(v[startsWith(v, "s_")], paste0("s_sx_1[", seq_len(k), "]"))
  expect_false(any(grepl("^b[[]", v)))
  # the columns are the basis weights the fit's b holds, in its order
  bi <- draws_index(ds)$b
  expect_identical(colnames(ds$draws)[bi], paste0("s_sx_1[", seq_len(k), "]"))
  expect_identical(v[1:4], c("b_Intercept", "bs_sx_1", "sds_sx_1", "sigma"))

  dg <- ord_sample(bf(ys ~ gp(x)), gaussian(), d)
  vg <- variables(dg)
  n_z <- length(dg$fit$frame$par_template$b)
  expect_identical(vg[startsWith(vg, "zgp_")],
                   paste0("zgp_gpx[", seq_len(n_z), "]"))
  expect_lt(which(vg == "sdgp_gpx"), which(vg == "lscale_gpx"))
  expect_lt(which(vg == "sigma"), min(which(startsWith(vg, "zgp_"))))
})

test_that("ar(cov = FALSE)'s states are brms's zerr", {
  d <- ord_data()
  ds <- ord_sample(bf(cnt ~ x + ar(t, g)), poisson(), d)
  v <- variables(ds)
  expect_identical(v[startsWith(v, "zerr[")],
                   paste0("zerr[", seq_len(nrow(d)), "]"))
  expect_false(any(grepl("^latent_state", v)))
  expect_identical(v[1:4], c("b_Intercept", "b_x", "ar[1]", "sderr"))
})

test_that("one response's draws matrices carry brms's dimnames", {
  d <- ord_data()
  ds <- ord_sample(bf(y ~ x), gaussian(), d)
  expect_identical(dimnames(posterior_predict(ds)), list(NULL, NULL))
  expect_identical(dimnames(predict(ds, summary = FALSE)), list(NULL, NULL))
  expect_identical(dimnames(predictive_error(ds)), list(NULL, NULL))
  expect_identical(dimnames(log_lik(ds)), list(NULL, NULL))
  # brms leaves the expectation's margins unnamed
  expect_null(dimnames(posterior_epred(ds)))
})

test_that("a categorical posterior_epred() leaves the rows unnamed, as brms", {
  d <- ord_data()
  set.seed(4)
  d$yc <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
  ds <- ord_sample(bf(yc ~ x), categorical(), d)
  pe <- posterior_epred(ds)
  # brms 2.23.0: list(NULL, NULL, c("a", "b", "c"))
  expect_identical(dimnames(pe), list(NULL, NULL, c("a", "b", "c")))
  expect_identical(dimnames(posterior_predict(ds)), list(NULL, NULL))
})
