# A missing value's draw under brms's name and on brms's scale: frm_sample()
# against brms 2.23.0's Ymi, on a log chart and on the window of von_mises,
# and from a frm() fit that left the row out of its objective (gated).

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

# through the compile cache, as test-default-priors-fit.R reads it
mi_models <- new.env()
mi_stan_model <- function(code) {
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  dir <- Sys.getenv("FRMTMB_STAN_CACHE", "")
  if (!nzchar(dir)) dir <- file.path(tempdir(), "frmtmb-stan-cache")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(key, ".rds"))
  # a program compiled in this session does not load again from its own RDS
  mod <- mi_models[[key]]
  if (is.null(mod) && file.exists(path)) mod <- readRDS(path)
  if (is.null(mod)) {
    mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
    stan_cache_write(mod, path)
  }
  mi_models[[key]] <- mod
  mod
}

mi_ks <- function(case, seed = 9700, from_fit = FALSE) {
  set.seed(seed)
  n <- 120
  d <- data.frame(x = stats::rnorm(n))
  eta <- 0.4 + 0.5 * d$x
  d$y <- switch(case,
    Gamma = stats::rgamma(n, 2.5, 2.5 / exp(eta)),
    von_mises = ((2 * atan(0.5 * eta) + stats::rnorm(n, 0, 0.5) + pi) %%
                   (2 * pi)) - pi)
  d$y[sample(n, 15)] <- NA
  ff <- switch(case, Gamma = stats::Gamma(link = "log"),
               von_mises = von_mises())
  bfam <- switch(case, Gamma = brms::brmsfamily("gamma", "log"),
                 von_mises = brms::brmsfamily("von_mises"))
  go <- function(x, ...) {
    suppressWarnings(suppressMessages(frm_sample(
      x, chains = 2, iter = 2000, seed = 11, refresh = 0, ...)))
  }
  ours <- if (from_fit) {
    go(suppressWarnings(frm(bf(y | mi() ~ x), family = ff, data = d)))
  } else {
    go(bf(y | mi() ~ x), family = ff, data = d)
  }
  fb <- brms::bf(y | mi() ~ x)
  pr <- brms::default_prior(fb, family = bfam, data = d)
  code <- brms::stancode(fb, family = bfam, data = d, prior = pr)
  sdat <- brms::standata(fb, family = bfam, data = d, prior = pr)
  sf <- suppressWarnings(rstan::sampling(mi_stan_model(code), data = sdat,
                                         chains = 2, iter = 2000, seed = 12,
                                         refresh = 0))
  ym <- as.matrix(rstan::extract(sf, "Ymi")$Ymi)
  dm <- as.matrix(as_draws_matrix(ours))
  cn <- paste0("Ymi[", seq_len(15), "]")
  expect_true(all(cn %in% colnames(dm)))
  expect_false(any(startsWith(colnames(dm), "miss_")))
  ks <- vapply(cn, function(k) {
    j <- match(k, cn)
    suppressWarnings(stats::ks.test(dm[, k], ym[, j])$statistic)
  }, 0)
  # the two-sample bound at the draws' own count, widened for their
  # autocorrelation by the factor the review read (median 0.017 to 0.021 at
  # 2000 each)
  list(ks = ks, bound = 1.36 * sqrt(2 / nrow(ym)), range = range(dm[, cn]))
}

test_that("Ymi matches brms's under a log chart", {
  skip_sampler()
  skip_unless_brms_fit()
  r <- mi_ks("Gamma")
  expect_lt(stats::median(r$ks), r$bound)
  expect_gt(r$range[1L], 0)
})

test_that("Ymi matches brms's in the window of von_mises", {
  skip_sampler()
  skip_unless_brms_fit()
  r <- mi_ks("von_mises")
  expect_lt(stats::median(r$ks), r$bound)
  expect_gte(r$range[1L], -pi)
  expect_lte(r$range[2L], pi)
})

test_that("frm_sample() of a frm() fit holds what the fit left out", {
  skip_sampler()
  skip_unless_brms_fit()
  r <- mi_ks("Gamma", from_fit = TRUE)
  expect_lt(stats::median(r$ks), r$bound)
})
