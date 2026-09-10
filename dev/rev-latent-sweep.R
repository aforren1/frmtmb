# Reviewer, item 2.3 claim 1 and item 7: re-run the lane's jitter sweep
# against a package built by THIS review, on the lane's own seeds.
#
# WHY. The lane disclosed that dev/latent-install.ps1 gated on
# PowerShell's `$?` and therefore skipped its install silently, so two
# probes ran against a stale library. It does not name which two. The
# cheapest way to settle that is not to reconstruct the timeline, it is
# to take the numbers again from a build this review made.
#
#   Rscript dev/rev-latent-sweep.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
cat("frmtmb.latent from:", find.package("frmtmb.latent")[[1L]], "\n")

d <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-2p3-repro81.rds")
dat <- d$dat
form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
fit <- frm(form, family = fam, data = dat)
ll0 <- as.numeric(logLik(fit))
GLOBAL <- -1087.99646521
cat("cold logLik:", format(ll0, digits = 12),
    " lane recorded -1096.09575602; relative gap:",
    format(abs(ll0 - d$ll_cold) / abs(ll0), digits = 3), "\n\n")

rows <- list()
for (j in c(0.5, 1, 2, 4, 8)) {
  for (r in 1:5) {
    sd <- 4100L + r
    ms <- suppressWarnings(hmm_starts(fit, n = 8, jitter = j, seed = sd))
    b <- as.numeric(logLik(ms$best))
    rows[[length(rows) + 1L]] <- data.frame(
      jitter = j, rep = r, seed = sd, best = b,
      recovered = as.integer(abs(b - GLOBAL) < 1e-4 * abs(GLOBAL)),
      converged = ms$n_converged, not_conv = ms$n_not_converged,
      err = ms$n_error, spread = ms$spread, modes = nrow(ms$modes))
  }
}
tb <- do.call(rbind, rows)
print(tb, row.names = FALSE, digits = 12)
cat("\nrecovered, by jitter (lane: 1/5 at 0.5, 5/5 at 1, 2, 4, 8):\n")
print(tapply(tb$recovered, tb$jitter, sum))
cat("\ndistinct optima, by jitter (lane: 1-2, 2, 2, 2, 2-4):\n")
print(tapply(tb$modes, tb$jitter, function(z) paste(range(z),
                                                    collapse = "-")))
cat("\nnot converged / errored anywhere:", sum(tb$not_conv), "/",
    sum(tb$err), "  (lane: 0 / 0)\n")
