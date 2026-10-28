# Reviewer probe 2b: posterior under brms's "" everywhere vs frm_sample's
# prior = "flat", on a three-dimensional random-effect block, where flat
# on the Cholesky factor differs from lkj(1). Arm "lkj1" is frmtmb's
# flat stack with lkj(1) on cor, for contrast.
# Usage: Rscript dev/review-reslot-p2b.R brms|flat|lkj1
# Data seed 41; sampler seed 1; 4 chains of 4000 (1000 warmup).
arm <- commandArgs(trailingOnly = TRUE)[1L]
if (identical(arm, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r14",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  suppressPackageStartupMessages(library(brms))
  cat("arm: brms", format(packageVersion("brms")), "\n")
} else {
  source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R",
         local = TRUE)
  suppressPackageStartupMessages(library(frmtmb.sample))
}
set.seed(41)
G <- 12L; n <- 8L
g <- rep(seq_len(G), each = n)
x <- rnorm(G * n); z <- rnorm(G * n)
y <- 1 + 0.5 * x + rnorm(G, 0, 0.6)[g] + rnorm(G, 0, 0.4)[g] * x +
  rnorm(G, 0, 0.4)[g] * z + rnorm(G * n)
d <- data.frame(y, x, z, g = factor(g))
f <- y ~ x + z + (1 + x + z | g)
vars <- c("cor_g__Intercept__x", "cor_g__Intercept__z", "cor_g__x__z",
          "sd_g__Intercept", "sd_g__x", "sd_g__z", "sigma")
summ <- function(v) {
  sprintf("mean %.4f sd %.4f E[v^2] %.4f q05 %.4f q50 %.4f q95 %.4f",
          mean(v), sd(v), mean(v^2), quantile(v, 0.05), quantile(v, 0.5),
          quantile(v, 0.95))
}
if (identical(arm, "brms")) {
  pr <- c(set_prior("", class = "Intercept"), set_prior("", class = "sigma"),
          set_prior("", class = "sd"), set_prior("", class = "cor"))
  m <- brm(bf(f), data = d, prior = pr, chains = 4, iter = 4000,
           warmup = 1000, seed = 1, refresh = 0, silent = 2)
  cat(grep("lprior|L_1|lkj", strsplit(stancode(m), "\n")[[1]], value = TRUE),
      sep = "\n")
  dr <- as.matrix(m)
  cat("divergent:", sum(rstan::get_divergent_iterations(m$fit)), "\n")
  ess <- posterior::summarise_draws(posterior::as_draws_df(m)[, vars],
                                    "ess_bulk")
} else {
  fit <- suppressWarnings(suppressMessages(frm(f, data = d)))
  prior <- "flat"
  if (identical(arm, "lkj1")) {
    prior <- frmtmb.sample:::default_priors_for(fit, flat = TRUE)
    for (k in seq_along(prior)) {
      if (identical(prior[[k]]$class, "cor")) prior[[k]]$dist$flat <- NULL
    }
  }
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, prior = prior, chains = 4, iter = 4000, warmup = 1000,
               seed = 1, refresh = 0)))
  print(prior_summary(ds))
  dr <- as.matrix(ds)
  cat("divergent:", sum(rstan::get_divergent_iterations(ds$stanfit)), "\n")
  ess <- posterior::summarise_draws(posterior::as_draws_df(dr[, vars]),
                                    "ess_bulk")
}
for (v in vars) cat(sprintf("  %-22s %s\n", v, summ(dr[, v])))
print(as.data.frame(ess))
