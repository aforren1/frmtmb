# frm_sample()'s posterior of a dispersion parameter against brms
# 2.23.0's, each under its own default priors, on small data where the
# prior matters: negbinomial shape (inv_gamma(0.4, 0.3)), beta phi and
# dirichlet phi (gamma(0.01, 0.01)), student nu (gamma(2, 0.1)). Through
# frmtmb.sample 0.19.0 all four were flat on the link scale, and the
# posteriors differed by 19 to 30 Monte Carlo standard errors (beta,
# dirichlet), or had no proper upper tail at all (negbinomial: q95 of
# log(shape) 33.4 against brms's 2.77; student: every draw of nu at its
# bound). Data and numbers: dev/sampleprior-agree-common.R,
# dev/sampleprior-agree.R and dev/sampleprior-findings.md, section 1.
# The negbinomial is compared with a quadrature of brms's model rather
# than with brms's chains (dpf_nb_reference()).
#
# Seed policy as test-loo.R's brms agreement: 2 chains of 2000 on seed 3
# on both sides. brms's program is compiled once into the stan cache
# (FRMTMB_STAN_CACHE), as core's log-density harness does.

skip_on_cran()

dpf_models <- new.env(parent = emptyenv())

dpf_brms_draws <- function(form, fam, data, par) {
  bform <- brms::bf(stats::as.formula(form))
  family <- brms::brmsfamily(fam)
  code <- brms::stancode(bform, data = data, family = family)
  sdat <- brms::standata(bform, data = data, family = family)
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  dir <- Sys.getenv("FRMTMB_STAN_CACHE", "")
  if (!nzchar(dir)) dir <- file.path(tempdir(), "frmtmb-stan-cache")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(key, ".rds"))
  # a program compiled in this session does not load again from its own
  # RDS, so the session keeps the object it has
  mod <- dpf_models[[key]]
  if (is.null(mod) && file.exists(path)) mod <- readRDS(path)
  if (is.null(mod)) {
    mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
    saveRDS(mod, path)
  }
  dpf_models[[key]] <- mod
  sf <- suppressWarnings(rstan::sampling(mod, data = sdat, chains = 2,
                                         iter = 2000, seed = 3,
                                         refresh = 0))
  as.matrix(sf, pars = par)[, 1L]
}

# log(parameter) on both sides: its posterior is the closer to normal.
# The means agree within 4 pooled Monte Carlo standard errors, and the
# central 90 percent intervals are within a quarter of each other's
# width; brms's own chains must have mixed for either to mean anything
dpf_expect_agree <- function(form, fam, data, par) {
  lb <- log(dpf_brms_draws(form, fam, data, par))
  ds <- suppressWarnings(suppressMessages(frm_sample(
    bf(stats::as.formula(form)), family = frmtmb::brmsfamily(fam),
    data = data, chains = 2, iter = 2000, seed = 3, refresh = 0)))
  lf <- log(as.matrix(ds)[, par])
  mc <- function(v) posterior::mcse_mean(matrix(v, ncol = 2))
  expect_lt(posterior::rhat(matrix(lb, ncol = 2)), 1.05)
  z <- (mean(lf) - mean(lb)) / sqrt(mc(lf)^2 + mc(lb)^2)
  expect_lt(abs(z), 4, label = paste(fam, par, "mean z"))
  w <- function(v) diff(stats::quantile(v, c(0.05, 0.95)))
  r <- unname(w(lf) / w(lb))
  expect_gt(r, 0.8, label = paste(fam, par, "width ratio"))
  expect_lt(r, 1.25, label = paste(fam, par, "width ratio"))
}

# brms's negbinomial posterior of log(shape) by quadrature, brms's model
# exactly: the centered Intercept under brms's own default (read off
# brms::default_prior()), a flat slope, shape ~ inv_gamma(0.4, 0.3) and
# the Jacobian of log(shape). brms's own chains are no reference here:
# on these data they stick in the Poisson limit on seeds 101 and 202
# (R-hat 2.24 and Inf, bulk ESS 3 and 2; dev/sampleprior-nbref.R)
dpf_nb_reference <- function(d) {
  pri <- with(as.data.frame(suppressMessages(brms::default_prior(
    y ~ x, data = d, family = brms::negbinomial()))),
    prior[class == "Intercept"])
  nums <- as.numeric(strsplit(sub("^student_t[(](.*)[)]$", "\\1", pri),
                              ",")[[1]])
  xc <- d$x - mean(d$x)
  fit <- suppressWarnings(stats::glm(y ~ xc, data = d,
                                     family = stats::poisson()))
  se <- sqrt(diag(stats::vcov(fit)))
  mu <- mean(d$y)
  ug <- seq(-8, 40, by = 0.05)
  # per log(shape), a grid over (Intercept, slope) as wide as that
  # shape's spread: the variance is mu + mu^2 / shape, so the standard
  # errors are the Poisson ones times sqrt(1 + mu / shape)
  lpu <- vapply(ug, function(u) {
    k <- sqrt(1 + mu / exp(u))
    ag <- seq(stats::coef(fit)[1] - 12 * k * se[1],
              stats::coef(fit)[1] + 12 * k * se[1], length.out = 49)
    bg <- seq(stats::coef(fit)[2] - 12 * k * se[2],
              stats::coef(fit)[2] + 12 * k * se[2], length.out = 49)
    ab <- expand.grid(a = ag, b = bg)
    eta <- outer(ab$a, rep(1, nrow(d))) + outer(ab$b, xc)
    ll <- rowSums(matrix(stats::dnbinom(rep(d$y, each = nrow(ab)),
                                        size = exp(u),
                                        mu = exp(as.vector(eta)),
                                        log = TRUE), nrow(ab)))
    lp <- ll + stats::dt((ab$a - nums[2]) / nums[3], nums[1], log = TRUE)
    m <- max(lp)
    m + log(sum(exp(lp - m))) + log(diff(ag)[1]) + log(diff(bg)[1])
  }, 1)
  lpu <- lpu + 0.4 * log(0.3) - lgamma(0.4) - 1.4 * ug - 0.3 / exp(ug) + ug
  w <- exp(lpu - max(lpu))
  w <- w / sum(w)
  cdf <- cumsum(w)
  q <- function(p) stats::approx(cdf, ug, xout = p, ties = "ordered")$y
  list(mean = sum(w * ug), q = c(q(0.05), q(0.95)))
}

test_that("negbinomial shape: brms's posterior under inv_gamma(0.4, 0.3)", {
  skip_unless_brms_fit()
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(20261014)
  d <- data.frame(x = stats::rnorm(15))
  d$y <- stats::rnbinom(15, mu = exp(1 + 0.3 * d$x), size = 0.3)
  ref <- dpf_nb_reference(d)
  ds <- suppressWarnings(suppressMessages(frm_sample(
    bf(y ~ x), family = frmtmb::brmsfamily("negbinomial"), data = d,
    chains = 2, iter = 2000, seed = 3, refresh = 0)))
  lf <- log(as.matrix(ds)[, "shape"])
  # the reference has no Monte Carlo error, so the yardstick is the
  # draws' own; on rellib-r9, flat on log(shape), z was 3.8 and the
  # width ratio 8.9 on this seed (8.3 and 8.2 on seed 101)
  z <- (mean(lf) - ref$mean) / posterior::mcse_mean(matrix(lf, ncol = 2))
  expect_lt(abs(z), 4)
  r <- unname(diff(stats::quantile(lf, c(0.05, 0.95))) / diff(ref$q))
  expect_gt(r, 0.8)
  expect_lt(r, 1.25)
})

test_that("beta phi: brms's posterior under gamma(0.01, 0.01)", {
  skip_unless_brms_fit()
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(20261015)
  d <- data.frame(x = stats::rnorm(10))
  mu <- stats::plogis(0.3 + 0.5 * d$x)
  d$y <- stats::rbeta(10, mu * 300, (1 - mu) * 300)
  dpf_expect_agree("y ~ x", "beta", d, "phi")
})

test_that("student nu: brms's posterior under gamma(2, 0.1)", {
  skip_unless_brms_fit()
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(20261016)
  d <- data.frame(x = stats::rnorm(20))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20)
  dpf_expect_agree("y ~ x", "student", d, "nu")
})

test_that("dirichlet phi: brms's posterior under gamma(0.01, 0.01)", {
  skip_unless_brms_fit()
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(20261017)
  d <- data.frame(x = stats::rnorm(10))
  e <- cbind(0, 0.3 * d$x, -0.2 + 0.2 * d$x)
  p <- exp(e) / rowSums(exp(e))
  g <- matrix(stats::rgamma(30, 300 * p), 10)
  d$Y <- g / rowSums(g)
  colnames(d$Y) <- c("a", "b", "c")
  dpf_expect_agree("Y ~ x", "dirichlet", d, "phi")
})
