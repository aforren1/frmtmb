# Item 7: bernoulli and poisson against brms on a design the lane did not use, brms's default
# priors both sides. Args: family. Seed 7300.
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(REV_WT, "dev/review-ncgp-stan-cache"))
suppressPackageStartupMessages(library(brms))
fam <- args[1]
set.seed(7300)
n <- 80L
x <- round(sort(runif(n, 0, 8)), 3)
eta <- 1.2 * sin(1.3 * x) + 0.2 * x - 0.5
d <- data.frame(x = x, y = if (fam == "bernoulli") rbinom(n, 1, plogis(eta)) else rpois(n, exp(eta + 1)))
fm <- if (fam == "bernoulli") brms::bernoulli() else stats::poisson()
frmfam <- if (fam == "bernoulli") frmtmb::bernoulli() else stats::poisson()
b <- brm(y ~ gp(x), d, family = fm, chains = 4, cores = 4, iter = 1500, seed = 73, refresh = 0, silent = 2, backend = "rstan", control = list(adapt_delta = 0.95, max_treedepth = 12))
sm <- posterior::summarise_draws(posterior::subset_draws(posterior::as_draws_df(b), variable = c("b_Intercept", "sdgp_gpx", "lscale_gpx")))
print(as.data.frame(sm)[, c("variable", "q5", "median", "q95", "rhat", "ess_bulk")], digits = 4)
ppe <- brms::posterior_linpred(b)
nw <- brms::posterior_linpred(b, newdata = data.frame(x = c(0.4, 3.31, 6.77)))
pr <- brms::prior_summary(b)
print(pr[pr$class %in% c("sdgp", "lscale", "Intercept"), c("prior", "class", "coef")])
pl <- get_prior(y ~ gp(x), d, family = fm)
map <- suppressWarnings(frm(y ~ gp(x), d, family = frmfam, prior = c(set_prior(pl$prior[pl$class == "Intercept"][1], class = "Intercept"),
        set_prior(pl$prior[pl$class == "sdgp" & nzchar(pl$prior)][1], class = "sdgp"), set_prior(pl$prior[pl$class == "lscale" & nzchar(pl$prior)][1], class = "lscale", coef = "gpx"))))
lp <- frm_linpred(map, se.fit = TRUE)
z <- (lp$fit - colMeans(ppe)) / apply(ppe, 2, sd)
sr <- lp$se.fit / apply(ppe, 2, sd)
lpn <- frm_linpred(map, newdata = data.frame(x = c(0.4, 3.31, 6.77)), se.fit = TRUE)
cat(sprintf("RNB %s frm MAP conv %d intercept %.3f sd %.3f rho(scaled) %.3f | eta |z| max %.3f median %.3f | se ratio %.3f to %.3f | new x z %s se ratio %s\n", fam, map$opt$convergence, map$estimates$beta[1], exp(map$estimates$theta[1]), exp(map$estimates$theta[2]), max(abs(z)), median(abs(z)), min(sr), max(sr),
            paste(sprintf("%.2f", (lpn$fit - colMeans(nw)) / apply(nw, 2, sd)), collapse = " "), paste(sprintf("%.2f", lpn$se.fit / apply(nw, 2, sd)), collapse = " ")))
cat("DONE\n")
