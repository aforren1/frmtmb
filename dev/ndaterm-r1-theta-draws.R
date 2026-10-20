# Lane ndaterm, fix round 1: the gaussian mixture of dev/ndaterm-r1-theta.R
# through frmtmb.sample's posterior_predict(newdata). Seed 7.
# Usage: Rscript dev/ndaterm-r1-theta-draws.R <arm>
arm <- commandArgs(TRUE)[1]
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(7)
n <- 300
d <- data.frame(x = rnorm(n))
k <- rbinom(n, 1, 0.2)
d$y <- ifelse(k == 1, 6, 0) + d$x + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 400, warmup = 200, refresh = 0,
             seed = 7)))
nd <- data.frame(x = c(-1, 0, 1))
set.seed(26)
p <- posterior_predict(ds, newdata = nd)
e <- posterior_epred(ds, newdata = nd)
cat("posterior_predict means", colMeans(p), "\nposterior_epred means",
    colMeans(e), "\n")
sh <- frmtmb.sample:::draws_fit_at(ds, 1, frmtmb.sample:::draws_index(ds))
cat("draw 1 theta1 response", frm_linpred(sh, newdata = nd[1, , drop = FALSE],
    dpar = "theta1", type = "response"), "\n")
