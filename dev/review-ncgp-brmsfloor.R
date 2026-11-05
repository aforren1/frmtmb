# Item 2: brms's own behavior on noise-free data and at y * 1000, and the closed form
# at brms's absolute line. Args: label ysc noise iter treedepth
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
suppressPackageStartupMessages(library(brms))
lab <- args[1]; ysc <- as.numeric(args[2]); noise <- as.numeric(args[3]); it <- as.integer(args[4]); td <- as.integer(args[5])
set.seed(4130)
n <- 30L
x <- sort(runif(n, 0, 10)); e <- rnorm(n)
y <- ysc * (sin(2 * x) + noise * e)
d <- data.frame(y = y, x = x)
t0 <- Sys.time()
b <- brm(y ~ gp(x), d, chains = 2, iter = it, seed = 7, refresh = 0, silent = 2, backend = "rstan",
         control = list(adapt_delta = 0.95, max_treedepth = td))
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
sm <- posterior::summarise_draws(posterior::subset_draws(posterior::as_draws_df(b), variable = c("sigma", "sdgp_gpx", "lscale_gpx", "b_Intercept")))
print(as.data.frame(sm)[, c("variable", "q5", "median", "q95", "rhat", "ess_bulk")], digits = 4)
np <- nuts_params(b)
cat(sprintf("RB %s ysc %g noise %g wall %.0f s, at depth %d of %d, divergent %d\n", lab, ysc, noise, el,
            sum(np$Value[np$Parameter == "treedepth__"] >= td), sum(np$Parameter == "treedepth__"), sum(np$Value[np$Parameter == "divergent__"])))
r <- rev_fit(y ~ gp(x), d)
g <- rev_gp_read(r$fit)
cat(sprintf("RB frm sigma %.4g sd %.4g rho(scaled by brms's?) %.4g ll %.4f conv %d w %d\n", g$sigma, g$sd, g$rho, g$ll, g$conv, length(r$warn)))
ep <- posterior_epred(b)
cat(sprintf("RB epred median abs error against y: %.3e; frm's fitted %.3e\n", median(abs(apply(ep, 2, median) - y)), median(abs(fitted(r$fit)[, 1] - y))))
cat("DONE\n")
