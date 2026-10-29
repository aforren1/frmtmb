# Lane meme: fit agreement of frm_sample() with brms 2.23.0 on me()
# models, brms's default priors on both sides. One arm per process:
#   Rscript dev/meme-agree.R brms   (brms fit)
#   Rscript dev/meme-agree.R lane   (frm_sample, lane library)
#   Rscript dev/meme-agree.R base   (frm_sample, rellib-r15)
# Each writes dev/meme-out/agree-<arm>.rds (per-variable mean, MCSE, sd
# and quantiles); dev/meme-agree-cmp.R prints the z-scores.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
if (identical(arm, "brms")) libs <- libs[2]
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(posterior))

agree_data <- function(seed = 2026, ng = 10, k = 6) {
  set.seed(seed)
  n <- ng * k
  g <- factor(rep(seq_len(ng), each = k))
  tx <- stats::rnorm(ng, 0.5, 1)
  tz <- 0.6 * tx + stats::rnorm(ng, 0, 0.8)
  d <- data.frame(g = g, xg = (tx + stats::rnorm(ng, 0, 0.4))[g],
                  zg = (tz + stats::rnorm(ng, 0, 0.4))[g], sg = 0.4)
  sx <- runif(n, 0.2, 0.5)
  ux <- stats::rnorm(n, 1, 0.8)
  d$x <- ux + stats::rnorm(n, 0, sx)
  d$sx <- sx
  d$y <- 1 + 0.6 * tx[g] - 0.5 * tz[g] + 0.4 * ux + stats::rnorm(n, 0, 0.6)
  d
}
d <- agree_data()
summ <- function(M) {
  out <- lapply(colnames(M), function(v) {
    x <- matrix(M[, v], ncol = 4)
    c(mean = mean(x), mcse = mcse_mean(x), sd = stats::sd(x),
      mcse_sd = mcse_sd(x),
      q05 = unname(stats::quantile(x, 0.05)),
      q50 = unname(stats::quantile(x, 0.5)),
      q95 = unname(stats::quantile(x, 0.95)),
      mcse_q05 = mcse_quantile(x, 0.05), mcse_q95 = mcse_quantile(x, 0.95))
  })
  structure(do.call(rbind, out), dimnames = list(colnames(M), NULL)) |>
    (\(m) { colnames(m) <- c("mean", "mcse", "sd", "mcse_sd", "q05", "q50",
                             "q95", "mcse_q05", "mcse_q95"); m })()
}
keep <- c("b_Intercept", "bsp_mexgsggrEQg", "bsp_mezgsggrEQg", "bsp_mexsx",
          "sigma", "meanme_mexg", "meanme_mezg", "meanme_mex", "sdme_mexg",
          "sdme_mezg", "sdme_mex", "corme_g__mexg__mezg", "Xme_mexg[1]",
          "Xme_mex[1]")
form_txt <- "y ~ me(xg, sg, gr = g) + me(zg, sg, gr = g) + me(x, sx)"
t0 <- proc.time()[["elapsed"]]
if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  fit <- brm(stats::as.formula(form_txt), data = d, chains = 4,
             iter = 4000, warmup = 1000, seed = 1, refresh = 0, cores = 4,
             save_pars = save_pars(latent = TRUE))
  M <- as_draws_matrix(fit)
  cat("brms divergent:", sum(rstan::get_divergent_iterations(fit$fit)),
      "\n")
} else {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  ds <- frm_sample(stats::as.formula(form_txt), data = d,
                   family = gaussian(), chains = 4, iter = 4000,
                   warmup = 1000, seed = 1, refresh = 0, cores = 1)
  M <- as.matrix(ds$draws)
  cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "from",
      find.package("frmtmb.sample"), "\n")
  print(as.data.frame(prior_summary(ds)))
  cat("divergent:", sum(rstan::get_divergent_iterations(ds$stanfit)), "\n")
}
cat("elapsed", proc.time()[["elapsed"]] - t0, "s\n")
miss <- setdiff(keep, colnames(M))
if (length(miss)) cat("missing columns:", miss, "\n")
s <- summ(M[, intersect(keep, colnames(M)), drop = FALSE])
print(round(s, 4))
saveRDS(s, file.path("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/meme-out",
                     paste0("agree-", arm, ".rds")))
