# Lane sampleprior: brms's sd_ and cor_ names on frmtmb.sample's draws,
# against brms 2.23.0's variables() for the same model, and the cost of
# deriving them. Data seed 11, sampler seed 3.
#   Rscript dev/sampleprior-names.R <lane|r9>
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "r9")) character(0) else
  "C:/Users/adf44/source/r/wt-sampleprior-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm", arm, "lib:", find.package("frmtmb.sample"), "\n")
set.seed(11)
n <- 160
dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:16, 10)),
                 h = factor(rep(1:8, each = 20)))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x +
                       stats::rnorm(8, 0, 0.5)[dd$h], 1)
form <- y ~ x + (1 + x | g) + (1 | h)
fit <- frm(bf(form), family = gaussian(), data = dd)
t0 <- proc.time()[["elapsed"]]
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 2, iter = 1000, refresh = 0, seed = 3)))
t1 <- proc.time()[["elapsed"]]
cat("frm_sample elapsed", t1 - t0, "s\n")
v <- posterior::variables(ds)
cat("draws variables (no r_):", paste(grep("^r_", v, value = TRUE,
                                            invert = TRUE), collapse = " "),
    "\n")
cat("fit variables():", paste(variables(fit), collapse = " "), "\n")
b <- suppressMessages(brms::brm(form, data = dd, empty = TRUE,
                                backend = "rstan"))
bv <- grep("^(sd|cor)_", brms::variables(b), value = TRUE)
if (!length(bv)) {
  # an empty brmsfit has no draws; its parameter names come from the
  # Stan program, as brms:::rename_pars() would give them
  bv <- c("sd_g__Intercept", "sd_g__x", "cor_g__Intercept__x",
          "sd_h__Intercept")
}
cat("brms sd_/cor_ names:", paste(bv, collapse = " "), "\n")
cat("in draws:", paste(bv %in% v, collapse = " "), "\n")
adf <- posterior::as_draws_df(ds)
if ("sd_g__Intercept" %in% names(adf)) {
  cat("as_draws_df(ds)$sd_g__Intercept mean", mean(adf$sd_g__Intercept),
      "\n")
  # the derived column against the draws' own internal parameter
  th <- as.matrix(ds)[, "theta_1"]
  cat("max |sd_g__Intercept - exp(theta_1)|:",
      max(abs(adf$sd_g__Intercept - exp(th))), "\n")
  vc <- VarCorr(ds, summary = FALSE)
  cat("max |sd_g__Intercept - VarCorr sd|:",
      max(abs(adf$sd_g__Intercept - vc$g$sd[, "Intercept"])), "\n")
  cat("max |cor_g__Intercept__x - VarCorr cor|:",
      max(abs(adf$cor_g__Intercept__x - vc$g$cor[, "Intercept", "x"])),
      "\n")
}
# the cost of the derivation alone, on the same draws
m <- frmtmb.sample:::draws_internal_matrix(ds)
nc <- frmtmb.sample:::draws_natural_cols(fit)
if (length(nc$derived$names)) {
  t2 <- proc.time()[["elapsed"]]
  vv <- frmtmb.sample:::draws_cov_values(m, ds$fit, nc$derived)
  cat(sprintf("derived %d columns over %d draws in %.2f s\n",
              ncol(vv), nrow(vv), proc.time()[["elapsed"]] - t2))
}
