# test-reparam.R's rp_cor_case() on one build: the centered and the
# non-centered chain of y ~ x + (x | g), saved for comparison across
# builds, with the fit the formula route starts from.
#   Rscript dev/moparity-rp-cmp.R base|lane
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
set.seed(31)
ng <- 12L
dd <- data.frame(x = stats::rnorm(120),
                 g = factor(rep(seq_len(ng), length.out = 120)))
u <- matrix(stats::rnorm(2 * ng), 2)
dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x + 0.7 * u[1, dd$g] +
                       0.4 * u[2, dd$g] * dd$x, 1)
mk <- function(ncp) {
  suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (x | g)), family = gaussian(), data = dd,
               chains = 1, iter = 1000, refresh = 0, seed = 33,
               reparameterize = ncp)))
}
cen <- mk(FALSE)
ncp <- mk(TRUE)
saveRDS(list(cen = cen$draws, ncp = ncp$draws, est = cen$fit$estimates,
             init = cen$stanfit@inits, step = rstan::get_sampler_params(cen$stanfit, inc_warmup = TRUE)[[1]][1:5, ],
             fit_par = cen$fit$opt$par, fit_obj = cen$fit$opt$objective),
        file.path("dev/moparity-log", paste0("rp-cmp-", arm, ".rds")))
cat("saved\n")
