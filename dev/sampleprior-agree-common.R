# Lane sampleprior: the four posterior-agreement cases, shared by the
# dev script and read by the gated test's construction. Each case is a
# small data set on which the dispersion prior matters, with the
# dispersion parameter brms reports under the same name.

agree_cases <- function() {
  list(
    # overdispersed counts on 15 rows: the likelihood of shape levels off
    # toward the Poisson as shape grows, so inv_gamma(0.4, 0.3) is what
    # bounds the upper tail
    negbinomial = local({
      set.seed(20261014)
      d <- data.frame(x = stats::rnorm(15))
      d$y <- stats::rnbinom(15, mu = exp(1 + 0.3 * d$x), size = 0.3)
      list(form = "y ~ x", fam = "negbinomial", data = d, par = "shape")
    }),
    # a precise beta response on 10 rows: gamma(0.01, 0.01) pulls a
    # phi near 300 down
    beta = local({
      set.seed(20261015)
      d <- data.frame(x = stats::rnorm(10))
      mu <- stats::plogis(0.3 + 0.5 * d$x)
      d$y <- stats::rbeta(10, mu * 300, (1 - mu) * 300)
      list(form = "y ~ x", fam = "beta", data = d, par = "phi")
    }),
    # gaussian errors on 20 rows: nu is unbounded above in the
    # likelihood, and gamma(2, 0.1) is what holds it
    student = local({
      set.seed(20261016)
      d <- data.frame(x = stats::rnorm(20))
      d$y <- 1 + 0.5 * d$x + stats::rnorm(20)
      list(form = "y ~ x", fam = "student", data = d, par = "nu")
    }),
    # a three-part composition on 10 rows, precision near 300
    dirichlet = local({
      set.seed(20261017)
      d <- data.frame(x = stats::rnorm(10))
      e <- cbind(0, 0.3 * d$x, -0.2 + 0.2 * d$x)
      p <- exp(e) / rowSums(exp(e))
      g <- matrix(stats::rgamma(30, 300 * p), 10)
      d$Y <- g / rowSums(g)
      colnames(d$Y) <- c("a", "b", "c")
      list(form = "Y ~ x", fam = "dirichlet", data = d, par = "phi")
    }))
}

agree_models <- new.env(parent = emptyenv())

# brms 2.23.0's own program for a case, with brms's default priors,
# compiled once into the content-addressed cache and sampled with
# brms's standata on the given seed
agree_brms <- function(cs, seed, chains, iter, cache) {
  bform <- brms::bf(stats::as.formula(cs$form))
  fam <- brms::brmsfamily(cs$fam)
  code <- brms::stancode(bform, data = cs$data, family = fam)
  sdat <- brms::standata(bform, data = cs$data, family = fam)
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", utils::packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache, paste0(key, ".rds"))
  # a program compiled in this session does not load again from its own
  # RDS (core's tests/testthat/helper-brms.R, brms_stan_model()), so the
  # session keeps the object it has
  mod <- agree_models[[key]] %||% if (file.exists(path)) readRDS(path) else {
    m <- rstan::stan_model(model_code = code, save_dso = TRUE)
    saveRDS(m, path)
    m
  }
  agree_models[[key]] <- mod
  sf <- rstan::sampling(mod, data = sdat, chains = chains, iter = iter,
                        seed = seed, refresh = 0)
  as.matrix(sf, pars = cs$par)[, 1L]
}

agree_frm <- function(cs, seed, chains, iter) {
  ds <- suppressWarnings(suppressMessages(frmtmb.sample::frm_sample(
    frmtmb::bf(stats::as.formula(cs$form)),
    family = frmtmb::brmsfamily(cs$fam), data = cs$data, chains = chains,
    iter = iter, seed = seed, refresh = 0)))
  as.matrix(ds)[, cs$par]
}

# the mean's Monte Carlo standard error from the draws of all chains
agree_mcse <- function(v, chains) {
  a <- matrix(v, ncol = chains)
  posterior::mcse_mean(a)
}
