# Lane `latent`, item 2.3: does `hmm_starts()` escape the mode the
# probe found?
#
# Probe D4's cold start converges 8.099 log-likelihood units below the
# global optimum with every diagnostic clean
# (dev/latent-2p3-repro81.R reproduces it to twelve digits). This script
# runs `hmm_starts()` on that exact fit at several jitter sizes, several
# seeds each, and reports how often it recovers. A remedy that has not
# been shown recovering from the failure it was built for is not a
# remedy.
#
#   Rscript dev/latent-2p3-starts-probe.R [n] [reps]
#
# Seeds: the probe data is seed 2026 (probe D4's own); hmm_starts()
# takes seeds 4100 + r for repeat r.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})

args <- commandArgs(trailingOnly = TRUE)
NSTART <- if (length(args) >= 1L) as.integer(args[[1L]]) else 8L
REPS <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L
JIT <- c(0.5, 1, 2, 4, 8)
GLOBAL <- -1087.99646521

d <- readRDS("dev/latent-2p3-repro81.rds")
dat <- d$dat
form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")
fit <- frm(form, family = fam, data = dat)
ll0 <- as.numeric(logLik(fit))
cat("cold start logLik:", format(ll0, digits = 12),
    " gap to the global optimum:", format(GLOBAL - ll0, digits = 6),
    "\n\n")

ms1 <- hmm_starts(fit, n = NSTART, jitter = 2, seed = 4101)
print(ms1)
cat("\n")

out <- "dev/latent-2p3-starts-probe.tsv"
cat("n\tjitter\trep\tseed\tbest\tgap_to_global\trecovered\tconverged",
    "\tnot_converged\terror\tspread\tspread_all\tmodes\tseconds\n",
    sep = "", file = out)
for (j in JIT) {
  for (r in seq_len(REPS)) {
    sd <- 4100L + r
    t0 <- Sys.time()
    ms <- hmm_starts(fit, n = NSTART, jitter = j, seed = sd)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    b <- as.numeric(logLik(ms$best))
    rec <- abs(b - GLOBAL) < 1e-4 * abs(GLOBAL)
    cat(paste(c(NSTART, j, r, sd,
                formatC(b, digits = 10, format = "f"),
                formatC(GLOBAL - b, digits = 6, format = "g"),
                as.integer(rec), ms$n_converged, ms$n_not_converged,
                ms$n_error,
                formatC(c(ms$spread, ms$spread_all), digits = 6,
                        format = "g"),
                nrow(ms$modes),
                formatC(secs, digits = 4, format = "g")),
              collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
    cat(sprintf("jitter %4.1f rep %d: best %.6f gap %+.4f rec %d",
                j, r, b, GLOBAL - b, as.integer(rec)),
        sprintf(" conv %d/%d nc %d err %d modes %d %.1fs\n",
                ms$n_converged, NSTART, ms$n_not_converged, ms$n_error,
                nrow(ms$modes), secs))
    flush(stdout())
  }
}
