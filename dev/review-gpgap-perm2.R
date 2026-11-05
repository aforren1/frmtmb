# Short check of fix round 1, item 1: is a cheap fix there? The sampler's factor over
# a fixed farthest-first order of the positions, patched in this process, one chain
# per process (a socket worker would not carry the patch). b = P L z, so the kernel
# and the model are the same; only the rounding and the basis of z move. Data as in
# review-gpgap-draw.R: seed 6601, n 15. Prints per-chain lines; the sum script reads
# the saved draws.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-perm.R <log> <arm: perm|same> <chain seed> <iter> <td> [n]
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
sseed <- as.integer(a[2])
iter <- as.integer(a[3])
td <- as.integer(a[4])
n <- if (length(a) >= 5) as.integer(a[5]) else 15L
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
if (arm == "perm") {
  maximin_D2 <- function(D2) {
    m <- nrow(D2)
    p <- which.min(rowSums(D2))
    while (length(p) < m) {
      d <- apply(D2[, p, drop = FALSE], 1L, min)
      p <- c(p, which.max(d))
    }
    p
  }
  reg <- frmtmb:::frmtmb_covstruct_registry
  e <- reg$reg$gp
  e$chol_L <- function(theta, blk) {
    K <- frmtmb:::gp_corr(theta, blk)
    D2 <- Reduce(`+`, blk[["aux_D2"]])
    p <- maximin_D2(D2)
    ip <- order(p)
    Lp <- t(chol(K[p, p]))
    Lp[ip, , drop = FALSE] * exp(theta[1])
  }
  reg$reg$gp <- e
}
set.seed(6601)
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + as.numeric(Sys.getenv("REV_S", "1e-3")) * e0)
t0 <- proc.time()[["elapsed"]]
told <- character()
fit <- withCallingHandlers(
  frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 1L, iter = iter,
             warmup = iter %/% 2, seed = sseed, refresh = 0, cores = 1,
             control = list(adapt_delta = 0.95, max_treedepth = td)),
  warning = function(w) {
    told <<- c(told, gsub("[\r\n]+", " ", substr(conditionMessage(w), 1, 80)))
    invokeRestart("muffleWarning")
  })
wall <- proc.time()[["elapsed"]] - t0
dd <- posterior::as_draws_df(fit)
sp <- rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE)[[1]]
saveRDS(as.data.frame(dd[, c("sigma", "sdgp_gpx", "lscale_gpx")]),
        sprintf("dev/review-gpgap-out/permdraw2-%s-%d-%d.rds", arm, n, sseed))
cat(sprintf("PERM %s n %d seed %d sigma q50 %.3e sdgp %.3f lscale %.3f stepsize %.2g at-depth %d of %d divergent %d mean depth %.2f wall %.0f\n",
            arm, n, sseed, median(dd$sigma), median(dd$sdgp_gpx), median(dd$lscale_gpx),
            mean(sp[, "stepsize__"]), sum(sp[, "treedepth__"] >= td), nrow(sp),
            sum(sp[, "divergent__"]), mean(sp[, "treedepth__"]), wall))
