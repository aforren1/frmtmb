# Reviewer, fix round 1, item 1: rerun the lane's arma comparison
# against brms's OWN draws (no weights) on a NEW seed. The arma case
# reads every term of c. Data as dev/acprior-innov.R's arma (ARMA(1, 1)
# 0.7 / 0.3, 12 groups of 6, seed 20261027); 4 chains of 3000, seed 9.
# frmtmb's sigma draws carried to brms's innovation sd with c computed
# here from the draws, independently of innov_logc().
# Usage: Rscript dev/review-acprior-r1-arma.R frm|brms|summary
# Log: dev/review-acprior-log/r1-arma-<arm>.txt
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
lg <- file.path(wt, "dev", "review-acprior-log")
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "frm")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13", ulib))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
set.seed(20261027)
ng <- 12; k <- 6
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ng), function(i) {
  as.numeric(stats::arima.sim(list(ar = 0.7, ma = 0.3), n = k))
}))
fo <- y ~ x + arma(time = t, gr = g, cov = TRUE)
CH <- 4; IT <- 3000; SEED <- 9
if (identical(arm, "frm")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
  fit <- frm(fo, data = d)
  w <- character(0)
  ds <- withCallingHandlers(
    suppressMessages(frm_sample(fit, chains = CH, iter = IT, seed = SEED,
                                refresh = 0)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  cat("warnings:", if (length(w)) w else "none", "\n")
  cat("divergent:", sum(unlist(lapply(rstan::get_sampler_params(
    ds$stanfit, inc_warmup = FALSE), function(m) m[, "divergent__"]))), "\n")
  mm <- as.matrix(ds)
  ar <- mm[, "ar[1]"]; ma <- mm[, "ma[1]"]
  cc <- sqrt((1 - ar^2) / (1 + ma^2 + 2 * ar * ma))
  saveRDS(cbind(`ar[1]` = ar, `ma[1]` = ma, sigma = mm[, "sigma"] * cc,
                b_Intercept = mm[, "b_Intercept"], b_x = mm[, "b_x"]),
          file.path(lg, "r1-arma-frm.rds"))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  mod <- rstan::stan_model(model_code = stancode(fo, data = d))
  sf <- rstan::sampling(mod, data = standata(fo, data = d), chains = CH,
                        iter = IT, seed = SEED, refresh = 0)
  cat("divergent:", rstan::get_num_divergent(sf), "\n")
  mm <- as.matrix(sf, pars = c("ar[1]", "ma[1]", "sigma", "b_Intercept",
                               "b[1]"))
  colnames(mm) <- c("ar[1]", "ma[1]", "sigma", "b_Intercept", "b_x")
  saveRDS(mm, file.path(lg, "r1-arma-brms.rds"))
} else {
  .libPaths(ulib)
  F <- readRDS(file.path(lg, "r1-arma-frm.rds"))
  B <- readRDS(file.path(lg, "r1-arma-brms.rds"))
  for (pn in colnames(B)) {
    f <- F[, pn]; b <- B[, pn]
    mf <- posterior::mcse_mean(matrix(f, ncol = CH))
    mb <- posterior::mcse_mean(matrix(b, ncol = CH))
    sf <- posterior::mcse_sd(matrix(f, ncol = CH))
    sb <- posterior::mcse_sd(matrix(b, ncol = CH))
    cat(sprintf("%-12s frm %.4f (sd %.4f) brms %.4f (sd %.4f) | mean z %6.2f | sd z %6.2f\n",
                pn, mean(f), sd(f), mean(b), sd(b),
                (mean(f) - mean(b)) / sqrt(mf^2 + mb^2),
                (sd(f) - sd(b)) / sqrt(sf^2 + sb^2)))
  }
}
