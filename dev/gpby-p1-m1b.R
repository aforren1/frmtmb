# Punch round 1, m1: the fit route's chain on the review's construction,
# each arm, with the sampler's own diagnostics: does every draw sit at
# one point, and from which start.
# Usage: Rscript dev/gpby-p1-m1b.R <lane|base> [seed]
arm <- commandArgs(TRUE)[1]
sd_seed <- as.integer(commandArgs(TRUE)[2] %||% 4L)
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"),
    "| seed", sd_seed, "\n")
options(mc.cores = 1)
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d$y <- 0.5 + sin(d$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
for (route in c("fit", "formula")) {
  ds <- withCallingHandlers(
    if (route == "fit") {
      frm_sample(fit, chains = 1, iter = 600, refresh = 0, seed = sd_seed)
    } else {
      frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 1,
                 iter = 600, refresh = 0, seed = sd_seed)
    },
    warning = function(w) {
      cat("  warning:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  m <- as.matrix(ds)
  s <- apply(m, 2, stats::sd)
  cat(sprintf("%s route: %d columns, %d with sd 0; sd of the first four %s\n",
              route, ncol(m), sum(s == 0),
              paste(sprintf("%.3g", s[1:4]), collapse = " ")))
  sp <- tryCatch(rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)[[1]],
                 error = function(e) NULL)
  if (!is.null(sp)) {
    cat(sprintf("  accept_stat mean %.3f, stepsize %.3g, divergent %d, treedepth max %d\n",
                mean(sp[, "accept_stat__"]), sp[1, "stepsize__"],
                sum(sp[, "divergent__"]), max(sp[, "treedepth__"])))
  }
}
