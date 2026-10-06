# Punch round 1, m1: which priors the fit route and the formula route of
# frm_sample() put on a gp() term, and whether the fit route's chain
# sticks on the review's construction (y ~ gp(x), 60 points, seed 5,
# chains = 1, iter = 600, seed = 4).
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d$y <- 0.5 + sin(d$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
run <- function(lab, ...) {
  ds <- suppressWarnings(suppressMessages(frm_sample(...)))
  cat("==", lab, "\n")
  ps <- as.data.frame(prior_summary(ds))
  print(ps[, intersect(c("prior", "class", "coef", "source"), names(ps))])
  m <- as.matrix(ds)
  cat(sprintf("sd of draws: sdgp_gpx %.4g, lscale_gpx %.4g, sigma %.4g\n",
              stats::sd(m[, "sdgp_gpx"]), stats::sd(m[, "lscale_gpx"]),
              stats::sd(m[, "sigma"])))
  set.seed(1)
  e <- posterior_epred(ds, newdata = data.frame(x = c(6.7, 7.4)))
  cat(sprintf("across-draw sd of epred at 6.7, 7.4: %.4g %.4g\n",
              stats::sd(e[, 1]), stats::sd(e[, 2])))
}
run("fit route", fit, chains = 1, iter = 600, refresh = 0, seed = 4)
run("formula route", bf(y ~ gp(x)), data = d, family = gaussian(),
    chains = 1, iter = 600, refresh = 0, seed = 4)
