# Review of lane readers: fitted() at a new level against brms 2.23.0.
# Data as dev/readers-item1.R: seed 4101 (poisson (1 | g)), 4103 (bernoulli).
# brms: default priors, 4 chains of 1000 post-warmup draws, sampling seed 3,
# compiled through the Stan cache in dev/stan-cache. posterior_epred() with
# allow_new_levels = TRUE under sample_new_levels = "uncertainty" (brms's
# default) and "gaussian", epred seed 21.
source("dev/review-readers-setup.R")
suppressPackageStartupMessages(library(brms))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
cat("brms", format(packageVersion("brms")), "\n")
cache <- "C:/Users/adf44/source/r/frmtmb-wt-readers/dev/stan-cache"
stan_cached <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache, paste0(key, ".rds"))
  cat("cache key", key, if (file.exists(path)) "hit" else "MISS", "\n")
  if (file.exists(path)) return(readRDS(path))
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}
brms_fit <- function(form, fam, data) {
  b <- brm(form, data = data, family = fam, empty = TRUE)
  mod <- stan_cached(stancode(b))
  b$fit <- rstan::sampling(mod, data = standata(b), chains = 4, iter = 2000,
                           seed = 3, refresh = 0)
  brms:::rename_pars(b)
}

set.seed(4101)
J <- 40; n <- 400
g <- factor(sample(seq_len(J), n, TRUE))
u <- rnorm(J, 0, 1)
x <- runif(n, -1, 1)
dd <- data.frame(x = x, g = g)
dd$y <- rpois(n, exp(0.5 + 0.3 * x + u[g]))
set.seed(4103)
dd$yb <- rbinom(n, 1, plogis(0.2 + 0.5 * x + 1.5 * u[g]))
nd <- data.frame(x = c(-0.5, 0, 0.8), g = "new")

for (cs in list(list("y", poisson()), list("yb", bernoulli()))) {
  fo <- as.formula(paste(cs[[1]], "~ x + (1 | g)"))
  f <- frm(frmtmb::bf(fo), family = cs[[2]], data = dd)
  ff <- fitted(f, newdata = nd, allow_new_levels = TRUE)
  f0 <- fitted(f, newdata = nd, re_formula = NA)
  b <- brms_fit(brms::bf(fo), cs[[2]]$family, dd)
  sdb <- as.matrix(b, variable = "sd_g__Intercept")[, 1]
  cat("\n==", cs[[2]]$family, " frm sd(g)", VarCorr(f)$g$sd[1, 1],
      " brms sd(g) mean", mean(sdb), "median", median(sdb), "\n")
  set.seed(21)
  eu <- posterior_epred(b, newdata = nd, allow_new_levels = TRUE)
  set.seed(21)
  eg <- posterior_epred(b, newdata = nd, allow_new_levels = TRUE,
                        sample_new_levels = "gaussian")
  e0 <- posterior_epred(b, newdata = nd, re_formula = NA)
  out <- data.frame(
    x = nd$x,
    frm_est = ff[, "Estimate"], frm_err = ff[, "Est.Error"],
    frm_q025 = ff[, "Q2.5"], frm_q975 = ff[, "Q97.5"],
    frm_reNA = f0[, "Estimate"],
    unc_mean = colMeans(eu), unc_sd = apply(eu, 2, sd),
    gau_mean = colMeans(eg), gau_sd = apply(eg, 2, sd),
    gau_q025 = apply(eg, 2, quantile, 0.025),
    gau_q975 = apply(eg, 2, quantile, 0.975),
    reNA_mean = colMeans(e0))
  print(out, digits = 5)
  cat("frm_est / gau_mean:", format(out$frm_est / out$gau_mean, digits = 4),
      "\n")
  cat("frm_est / unc_mean:", format(out$frm_est / out$unc_mean, digits = 4),
      "\n")
  cat("frm_err / gau_sd:  ", format(out$frm_err / out$gau_sd, digits = 4),
      "\n")
  cat("frm_err / unc_sd:  ", format(out$frm_err / out$unc_sd, digits = 4),
      "\n")
  # "uncertainty": the mean over the seen levels' effects, which is not the
  # normal integral when the levels' empirical law is not normal
  re <- ranef(b, summary = FALSE)$g[, , "Intercept"]
  fe <- fixef(b, summary = FALSE)
  emp <- sapply(nd$x, function(xx) {
    eta <- fe[, "Intercept"] + fe[, "x"] * xx
    mean(cs[[2]]$linkinv(eta + re))
  })
  cat("brms mean over seen levels' draws:", format(emp, digits = 5), "\n")
  # frmtmb at its own estimate, mean over its own seen levels' modes
  rf <- ranef(f)$g
  rf <- if (length(dim(rf)) == 3L) rf[, "Estimate", 1L] else if (is.matrix(rf)) rf[, 1L] else unlist(rf)
  fef <- fixef(f)[, "Estimate"]
  empf <- sapply(nd$x, function(xx) {
    mean(cs[[2]]$linkinv(fef[1] + fef[2] * xx + rf))
  })
  cat("frm mean over seen levels' modes:", format(empf, digits = 5), "\n")
  # how much of brms's sd the sd's own uncertainty carries: brms gaussian
  # with the sd fixed at each draw's posterior median
}
