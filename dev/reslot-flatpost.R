# Lane reslot, item 3: frm_sample(prior = "flat") against brms 2.23.0
# with "" in every slot (brms's (flat) on its own parameters). Three
# small models where the internal and the natural parameter differ:
#   re:  y ~ 1 + (1 | g), 6 groups of 5 (sd)
#   cor: y ~ x + (1 + x | g), 8 groups of 6 (a 2 x 2 correlation)
#   ar:  y ~ 1 + ar(t, g, cov = TRUE), 6 series of 5 (rho, sigma)
# Usage: Rscript dev/reslot-flatpost.R brms|before|after
# Seeds: data 11, 12, 13; sampler seed 1; 4 chains of 4000 (1000 warmup).
arm <- commandArgs(trailingOnly = TRUE)[1L]
if (identical(arm, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r14",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  suppressPackageStartupMessages(library(brms))
  cat("arm: brms", format(packageVersion("brms")), "\n")
} else {
  source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
  suppressPackageStartupMessages(library(frmtmb.sample))
}
mk <- function(seed, G, n, slope = FALSE, ar = FALSE) {
  set.seed(seed)
  g <- rep(seq_len(G), each = n)
  t <- rep(seq_len(n), G)
  x <- rnorm(G * n)
  e <- if (ar) {
    unlist(lapply(seq_len(G), function(j) {
      as.numeric(stats::arima.sim(list(ar = 0.5), n))
    }))
  } else rnorm(G * n)
  b0 <- rnorm(G, 0, 0.5)[g]
  b1 <- if (slope) rnorm(G, 0, 0.5)[g] else 0
  data.frame(y = 1 + (0.5 + b1) * x * slope + b0 * (!ar) + e, x = x,
             g = factor(g), t = t)
}
cases <- list(
  re = list(d = mk(11, 6, 5), f = y ~ 1 + (1 | g),
            vars = c("sd_g__Intercept", "sigma", "b_Intercept")),
  cor = list(d = mk(12, 8, 6, slope = TRUE), f = y ~ x + (1 + x | g),
             vars = c("cor_g__Intercept__x", "sd_g__Intercept", "sd_g__x")),
  ar = list(d = mk(13, 6, 5, ar = TRUE), f = y ~ 1 + ar(t, g, cov = TRUE),
            vars = c("ar[1]", "sigma")))
summ <- function(v) {
  sprintf("mean %.4f sd %.4f q05 %.4f q50 %.4f q95 %.4f", mean(v), sd(v),
          quantile(v, 0.05), quantile(v, 0.5), quantile(v, 0.95))
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n==", nm, "\n")
  if (identical(arm, "brms")) {
    pr <- c(set_prior("", class = "Intercept"),
            set_prior("", class = "sigma"))
    if (nm != "ar") pr <- c(pr, set_prior("", class = "sd"))
    if (nm == "cor") pr <- c(pr, set_prior("", class = "cor"))
    if (nm == "ar") pr <- c(pr, set_prior("", class = "ar"))
    m <- brm(bf(cs$f), data = cs$d, prior = pr, chains = 4, iter = 4000,
             warmup = 1000, seed = 1, refresh = 0, silent = 2)
    dr <- as.matrix(m)
    cat("divergent:", sum(rstan::get_divergent_iterations(m$fit)), "\n")
    for (v in cs$vars) cat(sprintf("  %-22s %s\n", v, summ(dr[, v])))
    if (nm == "ar") {
      cat(sprintf("  %-22s %s\n", "marginal sd",
                  summ(dr[, "sigma"] / sqrt(1 - dr[, "ar[1]"]^2))))
    }
  } else {
    fit <- suppressWarnings(suppressMessages(frm(cs$f, data = cs$d)))
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, prior = "flat", chains = 4, iter = 4000,
                 warmup = 1000, seed = 1, refresh = 0)))
    dr <- as.matrix(ds)
    cat("divergent:", sum(rstan::get_divergent_iterations(ds$stanfit)),
        "\n")
    if (nm == "ar") {
      rho <- dr[, grep("^ar", colnames(dr))[1L]]
      cat(sprintf("  %-22s %s\n", "ar[1]", summ(rho)))
      # frmtmb 0.76.0 reports the marginal sd as sigma
      cat(sprintf("  %-22s %s\n", "sigma (innovation)",
                  summ(dr[, "sigma"] * sqrt(1 - rho^2))))
      cat(sprintf("  %-22s %s\n", "marginal sd", summ(dr[, "sigma"])))
    } else {
      for (v in cs$vars) cat(sprintf("  %-22s %s\n", v, summ(dr[, v])))
    }
  }
}
