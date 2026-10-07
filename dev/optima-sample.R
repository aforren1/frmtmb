# Lane optima, item 1: the sampler on a mo() fit after the change of
# simplex coordinates. frmtmb.sample puts no density on a simplex (the
# vignette port's documented gap, dev/brms-vignettes/brms_monotonic.R),
# so the posterior in the coordinates is improper on both builds. This
# records what one chain of brms_monotonic's fit1 (fit route, as the
# port runs it) gives on each arm: divergences, the simplex's draws.
#   Rscript dev/optima-sample.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
dat <- data.frame(income, ls)
fit1 <- frm(bf(ls ~ mo(income)), data = dat, family = gaussian())
ms <- tryCatch(utils::getFromNamespace("mo_simplex", "frmtmb"),
               error = function(e) function(z) {
                 x <- exp(c(0, z))
                 x / sum(x)
               })
cat("ML simplex:", format(ms(fit1$estimates$zeta1), digits = 4), "\n")
s <- suppressWarnings(frm_sample(fit1, chains = 1, iter = 1000,
                                 warmup = 500, seed = 1, cores = 1,
                                 refresh = 0))
dr <- as.matrix(s)
z <- dr[, grep("zeta", colnames(dr)), drop = FALSE]
W <- t(apply(z, 1, ms))
sp <- rstan::get_sampler_params(s$stanfit %||% s$fit, inc_warmup = FALSE)
div <- tryCatch(sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
                error = function(e) NA)
cat("divergent transitions:", div, "\n")
cat("simplex draws: mean", format(colMeans(W), digits = 3), "; sd",
    format(apply(W, 2, stats::sd), digits = 3), "\n")
cat("zeta draws: range", format(range(z), digits = 4), "\n")
