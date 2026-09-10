# Lane `latent`, punch round 2: the shrink weight on seeds it was NOT
# tuned on, and where the good region actually is.
#
# Round one chose w = 0.5 as "the smallest weight that loses none" on
# the 200 seeds the sweep ran on. That is a constant picked as the
# unique zero of a loss curve on its own tuning set, and review round
# two showed it does not survive a fresh draw: 199 of 200 out of
# sample, with the miss 256 units deep and every diagnostic clean.
#
# This script re-derives that independently, on this lane's own run,
# and locates BOTH edges of the plateau so the weight can be put in the
# interior rather than on a boundary.
#
#   Rscript dev/latent-2p4-oos.R <first seed> <n seeds> <out.tsv>
#
# Default seeds are 20270401 to 20270600, the block review round two
# used, which is disjoint from this lane's 20260910 to 20261109. The
# poLCA start seeds are 930001 upward, disjoint from both the lane's
# 700001 to 700200 and the review's 910001 to 920xxx.
#
# THE REFERENCE IS ADJUDICATED, not taken on trust. `poLCA(nrep = 10)`
# can miss too, so the reference for a seed is the LARGEST
# log-likelihood anyone reached on it: poLCA's best of ten, or any of
# the starts below. A start "loses" a seed when it is more than 1e-6
# relative below that. The raw poLCA figure is written out too, so the
# adjudication can be redone from the file.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
source("dev/latent-2p4-startfix-fns.R")
suppressMessages(loadNamespace("poLCA"))

args <- commandArgs(trailingOnly = TRUE)
S0 <- if (length(args) >= 1L) as.integer(args[[1L]]) else 20270401L
NS <- if (length(args) >= 2L) as.integer(args[[2L]]) else 200L
OUT <- if (length(args) >= 3L) args[[3L]] else "dev/latent-2p4-oos.tsv"
WS <- c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
wnm <- paste0("ll_w", sub("[.]", "", format(WS)))
cat(paste(c("seed", "ll_polca", wnm, "s_polca", "s_frm"),
          collapse = "\t"), "\n", sep = "", file = OUT)

for (i in seq_len(NS)) {
  seed <- S0 + i - 1L
  s <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(s$Y))
  dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
  for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "),
    ") ~ x1 + x2"))
  set.seed(930000L + i)
  t0 <- Sys.time()
  pl <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                     maxiter = 20000, tol = 1e-12)
  s_pl <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  sl <- slice_kmeans(s$Y, nc, K)
  t0 <- Sys.time()
  lls <- vapply(WS, function(w) {
    f <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight, start = profiles_from(s$Y, nc, K, sl, w = w)))),
      silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
  }, numeric(1))
  s_frm <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cat(paste(c(seed, formatC(c(pl$llik, lls), digits = 6, format = "f"),
              formatC(c(s_pl, s_frm), digits = 4, format = "g")),
            collapse = "\t"), "\n", sep = "", file = OUT, append = TRUE)
  if (i %% 20L == 0L) { cat("  ", i, " seeds\n", sep = ""); flush(stdout()) }
}
cat("wrote", OUT, "\n")
