## brms 2.23.0's own names for every case of dev/brmsnames-1024-cases.R.
##
##   Rscript dev/brmsnames-1024-brms.R <case number> [...]
##
## Writes dev/brmsnames-1024-out/brms-<n>.rds. Nothing here is recalled:
## standata(), stancode() and default_prior() run without a compile, and
## variables(), fixef(), vcov(), coef(), hypothesis() and
## conditional_effects() run on a real brmsfit, built the way brms's own
## documentation of brm(empty = TRUE) builds one (help("brm"), fit8):
## an empty brmsfit, a stanfit of brms's Stan code and standata (one
## chain, 40 iterations, seed 1), then brms::rename_pars(). The compiled
## program is cached in dev/stan-cache under the key the test suite's
## brms_stan_model() uses (md5 of the code plus the rstan version).
## Needs R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ids <- as.integer(commandArgs(trailingOnly = TRUE))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(brms))
q(library(rstan))
source("dev/brmsnames-1024-cases.R")
d <- brmsnames_data()
cases <- brmsnames_cases(brms::bf, brms::mvbf, brms::set_rescor)
dir.create("dev/brmsnames-1024-out", showWarnings = FALSE)
cache <- normalizePath("dev/stan-cache", winslash = "/")

stan_model_cached <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache, paste0(key, ".rds"))
  if (file.exists(path)) {
    m <- try(readRDS(path), silent = TRUE)
    if (!inherits(m, "try-error")) return(m)
  }
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}

try_names <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}

for (i in ids) {
  cs <- cases[[i]]
  cat("case", i, cs$id, "\n")
  fam <- cs$family
  out <- list(id = cs$id, brms = format(packageVersion("brms")))
  sdat <- q(brms::standata(cs$form, data = d, family = fam))
  out$X <- lapply(sdat[grepl("^X(sp|cs|s)?(_|$)", names(sdat))], colnames)
  code <- q(brms::stancode(cs$form, data = d, family = fam))
  pr <- as.data.frame(q(brms::default_prior(cs$form, data = d,
                                            family = fam)))
  out$prior <- pr[, c("class", "coef", "group", "resp", "dpar", "nlpar")]
  # every prior row brms lists with a coefficient, set alone, must pass
  # brms's own validate_prior()
  rows <- which(nzchar(pr$coef))
  out$prior_ok <- vapply(rows, function(r) {
    a <- brmsnames_prior_args(pr[r, ])
    p <- do.call(brms::set_prior, a)
    v <- tryCatch(q(brms::validate_prior(p, cs$form, data = d, family = fam)),
                  error = function(e) e)
    if (inherits(v, "error")) conditionMessage(v) else "ok"
  }, "")
  fit <- q(brms::brm(cs$form, data = d, family = fam, empty = TRUE))
  mod <- stan_model_cached(code)
  sf <- q(rstan::sampling(mod, data = sdat, chains = 1, iter = 40,
                          warmup = 20, refresh = 0, seed = 1))
  fit$fit <- sf
  fit <- brms::rename_pars(fit)
  out$variables <- brms::variables(fit)
  out$fixef <- try_names(rownames(brms::fixef(fit)))
  out$vcov <- try_names(rownames(stats::vcov(fit)))
  out$summary <- try_names(rownames(summary(fit)$fixed))
  out$coef <- try_names(lapply(stats::coef(fit), function(a) dimnames(a)[[3]]))
  out$ranef <- try_names(lapply(brms::ranef(fit), function(a) dimnames(a)[[3]]))
  out$varcorr <- try_names(lapply(brms::VarCorr(fit), function(a) {
    rownames(a$sd)
  }))
  out$ce <- try_names(names(q(brms::conditional_effects(fit))))
  hv <- grep("^(b|bsp|bs|bcs|sds|sd)_", out$variables, value = TRUE)
  out$hyp <- vapply(hv, function(v) {
    cls <- sub("_.*$", "", v)
    rest <- sub("^[a-z]+_", "", v)
    h <- tryCatch(q(brms::hypothesis(fit, paste0(rest, " = 0"), class = cls)),
                  error = function(e) e)
    if (inherits(h, "error")) conditionMessage(h) else "ok"
  }, "")
  saveRDS(out, sprintf("dev/brmsnames-1024-out/brms-%02d.rds", i))
}
