# Shared helpers for lane loglik's brms probes: brms 2.23.0 fits whose
# draws are chosen parameter vectors (Stan's Fixed_param, one chain per
# vector), compiled through the same content-addressed cache the test
# helpers use (tests/testthat/helper-brms.R), so a program compiled here
# is reused by the tests.

.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))

ll_cache <- "C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/stan-cache"
ll_models <- new.env()

ll_stan_model <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  unlink(f)
  if (!is.null(ll_models[[key]])) return(ll_models[[key]])
  path <- file.path(ll_cache, paste0(key, ".rds"))
  if (file.exists(path)) {
    mod <- try(readRDS(path), silent = TRUE)
    if (!inherits(mod, "try-error")) {
      ll_models[[key]] <- mod
      return(mod)
    }
  }
  cat("COMPILE", key, "\n")
  mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(mod, path)
  ll_models[[key]] <- mod
  mod
}

ll_flat_prior <- function(bform, data, family, ...) {
  p <- brms::get_prior(bform, data = data, family = family, ...)
  p$prior <- ""
  p
}

# A brmsfit whose k-th draw is inits[[k]].
ll_fixed_fit <- function(bform, family, data, inits, ...) {
  prior <- ll_flat_prior(bform, data = data, family = family, ...)
  code <- brms::make_stancode(bform, data = data, family = family,
                              prior = prior, ...)
  sdat <- suppressMessages(brms::make_standata(bform, data = data,
                                               family = family,
                                               prior = prior, ...))
  bfit <- suppressMessages(brms::brm(bform, data = data, family = family,
                                     prior = prior, empty = TRUE, ...))
  sf <- suppressMessages(rstan::sampling(
    ll_stan_model(code), data = sdat, chains = length(inits), iter = 1,
    warmup = 0, algorithm = "Fixed_param", init = inits, refresh = 0,
    seed = 1))
  bfit$fit <- sf
  brms::rename_pars(bfit)
}

# The leave-one-out conditional of y_i given the other rows, computed
# directly from the joint covariance S (or scale matrix, for the t),
# without the precision matrix: the reference the probes compare with.
ll_dense_cond <- function(y, m, S, nu = NULL) {
  n <- length(y)
  vapply(seq_len(n), function(i) {
    o <- seq_len(n)[-i]
    A <- solve(S[o, o, drop = FALSE])
    r <- y[o] - m[o]
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% r)
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    if (is.null(nu)) return(dnorm(y[i], cm, sqrt(cv), log = TRUE))
    d2 <- as.numeric(t(r) %*% A %*% r)
    df <- nu + n - 1
    sc <- sqrt((nu + d2) / df * cv)
    dt((y[i] - cm) / sc, df, log = TRUE) - log(sc)
  }, 0)
}
