## brms 2.23.0's names for dev/review-brmsnames-cases.R, from standata(),
## default_prior() and a real brmsfit (empty fit + stanfit of brms's own
## code and data, 1 chain, 40 iterations, seed 1, rename_pars()).
##   Rscript dev/review-brmsnames-brms.R <case> ...
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ids <- as.integer(commandArgs(trailingOnly = TRUE))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(brms))
q(library(rstan))
source("dev/review-brmsnames-cases.R")
d <- rv_data()
cases <- rv_cases(brms::bf, brms::mvbf)
dir.create("dev/review-brmsnames-out", showWarnings = FALSE)
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
tn <- function(e) tryCatch(e, error = function(err) paste("ERROR:",
                                                         conditionMessage(err)))
for (i in ids) {
  cs <- cases[[i]]
  fam <- cs$family
  out <- list(id = cs$id)
  sdat <- tn(q(brms::standata(cs$form, data = d, family = fam)))
  if (is.character(sdat)) {
    out$error <- sdat
  } else {
    out$X <- lapply(sdat[grepl("^(X|Z)", names(sdat))], function(m) {
      if (is.matrix(m)) colnames(m)
    })
    out$X <- Filter(Negate(is.null), out$X)
    pr <- as.data.frame(q(brms::default_prior(cs$form, data = d,
                                              family = fam)))
    out$prior <- pr[, c("prior", "class", "coef", "group", "resp", "dpar",
                        "nlpar")]
    code <- q(brms::stancode(cs$form, data = d, family = fam))
    fit <- q(brms::brm(cs$form, data = d, family = fam, empty = TRUE))
    mod <- stan_model_cached(code)
    sf <- q(rstan::sampling(mod, data = sdat, chains = 1, iter = 40,
                            warmup = 20, refresh = 0, seed = 1))
    fit$fit <- sf
    fit <- tn(brms::rename_pars(fit))
    if (is.character(fit)) {
      out$error <- fit
    } else {
      out$variables <- brms::variables(fit)
      out$fixef <- tn(rownames(brms::fixef(fit)))
      out$vcov <- tn(rownames(stats::vcov(fit)))
    }
  }
  saveRDS(out, sprintf("dev/review-brmsnames-out/brms-%02d.rds", i))
  cat("== case", i, cs$id, "\n")
  str(out[setdiff(names(out), "prior")])
  if (!is.null(out$prior)) print(out$prior[nzchar(out$prior$coef), -1])
}
