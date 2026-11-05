# M5: a nonlinear dpar beside a linear mu with a latent residual, the identity with
# brms's log_prob. Seed 167, the test file's data. Usage: Rscript dev/latres-m5.R
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({library(testthat); library(frmtmb); library(rstan)})
env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-brms.R", "helper-brms-methods.R", "helper-warnings.R")) {
  sys.source(file.path("tests/testthat", f), envir = env)
}
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30)); d <- d[order(d$g, d$t), ]
n <- nrow(d); d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
d$count <- rpois(n, exp(0.8 + 0.4 * d$x + e))
fz <- count ~ x + ar(t, g, cov = TRUE)
arms <- list(
  list(bf(fz) + nlf(shape ~ exp(a)) + lf(a ~ 1), negbinomial(),
       brms::bf(fz) + brms::nlf(shape ~ exp(a)) + brms::lf(a ~ 1), brms::brmsfamily("negbinomial")),
  list(bf(fz) + nlf(zi ~ exp(a) / (1 + exp(a))) + lf(a ~ 1), zero_inflated_poisson(),
       brms::bf(fz) + brms::nlf(zi ~ exp(a) / (1 + exp(a))) + brms::lf(a ~ 1), brms::brmsfamily("zero_inflated_poisson")))
for (a in arms) {
  fit <- suppressMessages(suppressWarnings(frm(a[[1]], family = a[[2]], data = d)))
  r <- env$brms_lp_check(a[[3]], a[[4]], d, fit, joint = TRUE)
  cat(a[[4]]$family, "ours", r$ours, "brms", r$lp, "constant", r$measured_const, "\n")
}
