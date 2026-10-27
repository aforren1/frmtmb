# Reviewer: a second seed for the ar comparison of
# dev/review-acprior-direct.R, to see whether the sd gap (z -2.45 on
# seed 3) is noise. Same data; seed 7; 4 chains of 6000.
# Usage: Rscript dev/review-acprior-rerun.R frm|brms|summary
# Log: dev/review-acprior-log/rerun-<arm>.txt
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
lg <- file.path(wt, "dev", "review-acprior-log")
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "frm")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13", ulib))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
set.seed(20261027)
ng <- 8; k <- 6
da <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                 x = rnorm(ng * k))
e <- numeric(nrow(da))
for (i in seq_len(nrow(da))) {
  e[i] <- if (da$t[i] == 1) rnorm(1) else 0.9 * e[i - 1] +
    rnorm(1, 0, sqrt(1 - 0.81))
}
da$y <- 1 + 0.5 * da$x + e
fo <- y ~ x + ar(time = t, gr = g, cov = TRUE)
CH <- 4; IT <- 6000; SEED <- 7
if (identical(arm, "frm")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
  fit <- frm(fo, data = da)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = CH, iter = IT, seed = SEED, refresh = 0)))
  x <- as.matrix(ds)[, "ar[1]"]
  saveRDS(x, file.path(lg, "rerun-frm.rds"))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  code <- stancode(fo, data = da)
  old <- "lprior += student_t_lpdf(sigma | 3, 0, 2.5)"
  new <- paste0("lprior += student_t_lpdf(sigma / sqrt(1 - square(ar[1])) ",
                "| 3, 0, 2.5) - 0.5 * log1m(square(ar[1]))")
  code2 <- sub(old, new, code, fixed = TRUE)
  mod <- rstan::stan_model(model_code = code2)
  sf <- rstan::sampling(mod, data = standata(fo, data = da), chains = CH,
                        iter = IT, seed = SEED, refresh = 0)
  saveRDS(as.matrix(sf, pars = "ar[1]")[, 1], file.path(lg, "rerun-brms.rds"))
} else {
  .libPaths(ulib)
  A <- readRDS(file.path(lg, "rerun-frm.rds"))
  D <- readRDS(file.path(lg, "rerun-brms.rds"))
  st <- function(x) {
    m <- matrix(x, ncol = CH)
    c(mean = mean(x), mcse_mean = unname(posterior::mcse_mean(m)),
      sd = sd(x), mcse_sd = unname(posterior::mcse_sd(m)),
      q05 = unname(quantile(x, 0.05)),
      mcse_q05 = unname(posterior::mcse_quantile(m, 0.05)),
      q95 = unname(quantile(x, 0.95)),
      mcse_q95 = unname(posterior::mcse_quantile(m, 0.95)))
  }
  sa <- st(A); sb <- st(D)
  for (s in c("mean", "sd", "q05", "q95")) {
    z <- (sa[[s]] - sb[[s]]) / sqrt(sa[[paste0("mcse_", s)]]^2 +
                                      sb[[paste0("mcse_", s)]]^2)
    cat(sprintf("seed 7 ar %s: frm %.4f brms(sigma_m prior) %.4f z %.2f\n",
                s, sa[[s]], sb[[s]], z))
  }
}
