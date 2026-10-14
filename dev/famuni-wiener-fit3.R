# Fit agreement under brms's own ndt parameterization: ndt ~ 1, which
# brms puts on a log link with a flat prior on its intercept, against
# frmtmb.eam's wiener(link_ndt = "log"). Same data as
# dev/famuni-wiener-fit2.R (seed 20261013); brms seed 20261013, rstan,
# 4 chains x 1000 post-warmup draws.
# Output: dev/famuni-out/wiener-fit3.txt; brmsfit cached in
# dev/famuni-out/wiener-brmsfit3.rds.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
suppressMessages(library(brms))
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out"
seed <- 20261013L

## data: identical construction to dev/famuni-wiener-fit2.R
set.seed(seed)
n <- 400
x <- rnorm(n)
truth <- c(b_Intercept = 0.6, b_x = 0.8, b_bs_Intercept = 0.25,
           b_bs_x = 0.2, ndt = 0.3, bias = 0.5)
mu <- truth[["b_Intercept"]] + truth[["b_x"]] * x
bs <- exp(truth[["b_bs_Intercept"]] + truth[["b_bs_x"]] * x)
sims <- lapply(seq_len(n), function(i) {
  RWiener::rwiener(1, alpha = bs[i], tau = truth[["ndt"]],
                   beta = truth[["bias"]], delta = mu[i])
})
dat <- data.frame(x = x, rt = vapply(sims, function(s) s$q, 0),
                  d = vapply(sims, function(s) as.integer(s$resp == "upper"),
                             0L))
cat("n", n, "min rt", min(dat$rt), "share upper", mean(dat$d), "\n")

fit <- frm(frmtmb::bf(rt | dec(d) ~ x, bs ~ x, ndt ~ 1),
           family = frmtmb.eam::wiener(link_ndt = "log"), data = dat)
cat("frmtmb.eam ML: convergence", fit$opt$convergence, "logLik",
    format(as.numeric(logLik(fit)), digits = 10), "\n")
print(summary(fit))
fit0 <- frm(frmtmb::bf(rt | dec(d) ~ x, bs ~ x),
            family = frmtmb.eam::wiener(), data = dat)
cat(sprintf("default-link fit on the same data: logLik %.10f (difference %.3e)\n",
            logLik(fit0), as.numeric(logLik(fit)) - as.numeric(logLik(fit0))))

bpath <- file.path(out_dir, "wiener-brmsfit3.rds")
if (file.exists(bpath)) {
  bfit <- readRDS(bpath)
  cat("\nbrms fit read from cache\n")
} else {
  t_b <- system.time(bfit <- brm(brms::bf(rt | dec(d) ~ x, bs ~ x, ndt ~ 1),
                                 family = brms::wiener(), data = dat,
                                 chains = 4, iter = 2000, warmup = 1000,
                                 seed = seed, backend = "rstan",
                                 refresh = 0, silent = 2))
  cat("\nbrms fit in", t_b[["elapsed"]], "s\n")
  saveRDS(bfit, bpath)
}
print(prior_summary(bfit))
print(summary(bfit))
dr <- as_draws_df(bfit)
cat("max Rhat", max(rhat(bfit), na.rm = TRUE), " min bulk ESS",
    min(posterior::summarise_draws(dr)$ess_bulk, na.rm = TRUE), "\n")
cat("divergences", sum(nuts_params(bfit)$Value[
  nuts_params(bfit)$Parameter == "divergent__"]), "\n")

fe <- frmtmb::fixef(fit)
print(fe)
ix <- c(Intercept = "b_Intercept", x = "b_x", bs_Intercept = "b_bs_Intercept",
        bs_x = "b_bs_x", ndt_Intercept = "b_ndt_Intercept")
stopifnot(all(names(ix) %in% rownames(fe)))
sp <- summary(fit)$spec_pars
tab <- data.frame(
  par = c(names(ix), "bias"),
  ml = c(fe[names(ix), "Estimate"], sp["bias", "Estimate"]),
  se = c(fe[names(ix), "Est.Error"], sp["bias", "Est.Error"]),
  post_median = c(vapply(ix, function(p) stats::median(dr[[p]]), 0),
                  stats::median(dr$bias)),
  post_sd = c(vapply(ix, function(p) stats::sd(dr[[p]]), 0),
              stats::sd(dr$bias)))
tab$z_shift <- (tab$ml - tab$post_median) / tab$post_sd
tab$se_over_sd <- tab$se / tab$post_sd
cat("\n-- ML against posterior median, in posterior SD\n")
print(tab, digits = 4, row.names = FALSE)
cat("max |z_shift|", max(abs(tab$z_shift)), "\n")
