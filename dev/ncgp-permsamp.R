# The factor's order on the sampling route, on this lane's tape: the design gpgap's
# review ran (dev/review-gpgap-perm.R in that tree), seed 6601, true sigma 1e-3, one
# chain per process since a socket worker would not carry the patch. With a chain
# seed per process; dev/ncgp-permsamp-sum.R reads the saved draws across chains.
#   Rscript dev/ncgp-permsamp.R <sorted|perm> <chain seed> <iter> <treedepth> [n]
a <- commandArgs(TRUE)
arm <- a[1]
sseed <- as.integer(a[2])
iter <- as.integer(a[3])
td <- as.integer(a[4])
n <- if (length(a) >= 5) as.integer(a[5]) else 15L
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb")), "arm", arm, "\n")
if (arm == "perm") source("dev/ncgp-perm-patch.R")
set.seed(6601)
x <- round(stats::runif(n, 0, 6), 3)
e0 <- stats::rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + 1e-3 * e0)
t0 <- proc.time()[["elapsed"]]
fit <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 1L, iter = iter,
             warmup = iter %/% 2L, seed = sseed, refresh = 0, cores = 1L,
             control = list(adapt_delta = 0.95, max_treedepth = td))))
wall <- proc.time()[["elapsed"]] - t0
dr <- as.data.frame(posterior::as_draws_df(fit))[, c("sigma", "sdgp_gpx", "lscale_gpx")]
sp <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)[[1L]]
dir.create("dev/ncgp-out/permsamp", showWarnings = FALSE)
saveRDS(dr, sprintf("dev/ncgp-out/permsamp/%s-%d-%d.rds", arm, n, sseed))
cat(sprintf("NCGPPS %s n %d seed %d sigma q50 %.3e sdgp %.3f lscale %.3f stepsize %.3g at-depth %d of %d divergent %d mean leapfrog %.0f wall %.0fs\n",
            arm, n, sseed, stats::median(dr$sigma), stats::median(dr$sdgp_gpx),
            stats::median(dr$lscale_gpx), sp[1L, "stepsize__"],
            sum(sp[, "treedepth__"] >= td), nrow(sp), sum(sp[, "divergent__"]),
            mean(sp[, "n_leapfrog__"]), wall))
