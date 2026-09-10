# Reviewer: guard3's case F and the d4 probe, both re-run against the
# F1 build, so the proposed fix is measured rather than argued.
#
#   Rscript dev/rev-latent-fixcheck.R <library>

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
lib <- commandArgs(trailingOnly = TRUE)[[1L]]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
cat("frmtmb.latent from:", find.package("frmtmb.latent")[[1L]], "\n")

small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}

cat("\n-- the six unimodal fits (guard3 case F) --\n")
res <- do.call(rbind, lapply(1:6, function(k) {
  dd <- small_data(seed = 4500L + k)
  f <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
  ms <- suppressWarnings(hmm_starts(f, n = 8, jitter = 2, seed = 8800L + k))
  out <- utils::capture.output(print(ms))
  data.frame(data_seed = 4500L + k, modes = nrow(ms$modes),
             best_is_original = isTRUE(all.equal(
               as.numeric(logLik(ms$best)), ms$original_logLik)),
             says_local_optimum = any(grepl("found a local optimum", out)),
             says_original_best =
               any(grepl("original fit was the best found", out)))
}))
print(res, row.names = FALSE)
cat("false alarms with the fix:", sum(res$says_local_optimum), "of 6",
    " (without it: 6 of 6)\n")
cat("best left as the incumbent:", sum(res$best_is_original), "of 6\n")

cat("\n-- the d4 probe, where the message must still fire --\n")
d <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/latent-2p3-repro81.rds")
fit <- frm(bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf)),
           family = hmm(K = 2, gaussian(), time = t, group = ID,
                        init = "stationary"),
           data = d$dat)
ms <- suppressWarnings(hmm_starts(fit, n = 6, jitter = 2, seed = 4101))
out <- utils::capture.output(print(ms))
cat("cold logLik :", format(ms$original_logLik, digits = 12), "\n")
cat("best logLik :", format(as.numeric(logLik(ms$best)), digits = 12), "\n")
cat("gap         :",
    format(as.numeric(logLik(ms$best)) - ms$original_logLik, digits = 8),
    "\n")
cat("modes       :", nrow(ms$modes), "\n")
cat("still says 'found a local optimum':",
    any(grepl("found a local optimum", out)), "\n")
cat("margin over the new threshold, as a ratio:",
    format((as.numeric(logLik(ms$best)) - ms$original_logLik) /
             (1e-6 * abs(ms$original_logLik)), digits = 5), "\n")
