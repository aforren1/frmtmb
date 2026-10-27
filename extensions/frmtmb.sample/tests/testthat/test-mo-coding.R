# mo() of an integer on the sampler, coded from its minimum as brms
# 2.23.0 codes it (get_mo_values()). Up to frmtmb 0.75.0 the integers
# 1..4 had a simplex of 4 weights, one for the category 0 that no row is
# in, and frm_sample() sampled that model: Intercept 0.912 (sd 1.330)
# against brms's 2.098 (0.117), with only ESS warnings
# (dev/mocode-findings.md, seed 11).

mo_code_data <- function() {
  set.seed(5)
  n <- 300
  d <- data.frame(k = sample(1:4, n, replace = TRUE))
  d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
  d
}

test_that("a 1..4 mo() is drawn as brms's 3 simo_ weights", {
  skip_on_cran()
  skip_sampler()
  d <- mo_code_data()
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ mo(k), data = d, chains = 1, iter = 300, seed = 1, refresh = 0)),
    "Effective Samples Size")
  expect_identical(grep("^simo_", variables(ds), value = TRUE),
                   paste0("simo_mok1[", 1:3, "]"))
  # new data are coded from the fitted minimum, and a value outside the
  # fitted range is refused, as brms's validate_newdata() refuses it
  e1 <- posterior_epred(ds, newdata = data.frame(k = 1:4))
  e2 <- posterior_epred(ds, newdata = data.frame(k = 2:4))
  expect_identical(unname(e1[, 2:4]), unname(e2))
  expect_error(posterior_epred(ds, newdata = data.frame(k = 0L)),
               "outside the fitted range 1..4")
  expect_error(posterior_predict(ds, newdata = data.frame(k = 5L)),
               "outside the fitted range 1..4")
})

mo_code_models <- new.env(parent = emptyenv())

# brms's own program on the same data, compiled once into the stan cache
# (FRMTMB_STAN_CACHE), as test-default-priors-fit.R compiles its own
mo_code_brms_draws <- function(d) {
  code <- brms::stancode(y ~ mo(k), data = d)
  sdat <- brms::standata(y ~ mo(k), data = d)
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  dir <- Sys.getenv("FRMTMB_STAN_CACHE", "")
  if (!nzchar(dir)) dir <- file.path(tempdir(), "frmtmb-stan-cache")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(key, ".rds"))
  mod <- mo_code_models[[key]]
  if (is.null(mod) && file.exists(path)) mod <- readRDS(path)
  if (is.null(mod)) {
    mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
    saveRDS(mod, path)
  }
  mo_code_models[[key]] <- mod
  # brms's own chains: rstan's ESS warnings only, by name
  sf <- allow_warnings(rstan::sampling(mod, data = sdat, chains = 2,
                                       iter = 2000, seed = 3, refresh = 0),
                       "Effective Samples Size")
  m <- as.matrix(sf, pars = c("b_Intercept", "bsp", "sigma", "simo_1"))
  colnames(m) <- c("b_Intercept", "bsp_mok", "sigma",
                   paste0("simo_mok1[", 1:3, "]"))
  m
}

test_that("frm_sample() of a 1..4 mo() agrees with brms", {
  skip_unless_brms_fit()
  skip_sampler()
  skip_if_not_installed("posterior")
  d <- mo_code_data()
  mb <- mo_code_brms_draws(d)
  # a divergence or any other warning would fail the test; on 0.75.0's
  # coding the ESS warnings were the only sign of the wrong model
  ds <- allow_warnings(suppressMessages(frm_sample(
    y ~ mo(k), data = d, chains = 2, iter = 2000, seed = 3, refresh = 0)),
    "Effective Samples Size")
  mf <- as.matrix(ds)
  expect_true(all(colnames(mb) %in% colnames(mf)))
  mc <- function(v) posterior::mcse_mean(matrix(v, ncol = 2))
  for (p in colnames(mb)) {
    # an absolute R-hat bound, the convention of
    # test-default-priors-fit.R: it only checks that brms's reference
    # chains mixed, and is not a tolerance on an agreement
    expect_lt(posterior::rhat(matrix(mb[, p], ncol = 2)), 1.05, label = p)
    # the posterior means agree within 4 pooled Monte Carlo standard
    # errors, and the spreads within a quarter of each other
    z <- (mean(mf[, p]) - mean(mb[, p])) / sqrt(mc(mf[, p])^2 +
                                                  mc(mb[, p])^2)
    expect_lt(abs(z), 4, label = paste(p, "mean z"))
    r <- stats::sd(mf[, p]) / stats::sd(mb[, p])
    expect_gt(r, 0.8, label = paste(p, "sd ratio"))
    expect_lt(r, 1.25, label = paste(p, "sd ratio"))
  }
})
