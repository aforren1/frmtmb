# Lane `latent`, punch round 2: does the recipe `?lca` prints for
# reproducing a 0.2.2 fit actually reproduce it?
#
# `lca()`'s starting values changed in 0.3.0 and the old rule is gone
# from the package, so a user with a published 0.2.2 result cannot refit
# it. The help page now carries the old rule as ten lines of R. A recipe
# on a help page is a claim, so this checks it against the rule as the
# lane's own sweep implemented it (`slice_score()` plus
# `profiles_from(w = 0)` in `dev/latent-2p4-startfix-fns.R`), and
# against the log-likelihood the 200-replicate sweep recorded for the
# old start.
#
#   Rscript dev/latent-2p4-oldstart.R [seed ...]

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")
source("dev/latent-2p4-startfix-fns.R")

args <- commandArgs(trailingOnly = TRUE)
SEEDS <- if (length(args)) as.integer(args) else
  c(20260910L, 20260970L, 20261013L)
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
ref <- utils::read.delim("dev/latent-2p4-lca.tsv", check.names = FALSE)

# VERBATIM from the block in ?lca. Anything changed here has to be
# changed there.
old_start <- function(Y, K) {
  ncat <- apply(Y, 2, max, na.rm = TRUE)
  sc <- rowMeans(sweep(Y - 1, 2, pmax(ncat - 1, 1), "/"), na.rm = TRUE)
  sl <- pmin(1 + ((rank(sc, ties.method = "first") - 1) * K) %/%
               nrow(Y), K)
  stats::setNames(lapply(seq_len(ncol(Y)), function(j) {
    unlist(lapply(seq_len(K), function(k) {
      p <- tabulate(Y[sl == k, j], nbins = ncat[j]) + 1
      p <- p / sum(p)
      log(p[-1]) - log(p[1])
    }))
  }), paste0("pi", seq_len(ncol(Y))))
}

for (seed in SEEDS) {
  s <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(s$Y))
  a <- old_start(s$Y, K)
  b <- profiles_from(s$Y, nc, K, slice_score(s$Y, nc, K), w = 0)
  d <- max(abs(unlist(a) - unlist(b)))
  f <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
        control = tight, start = a)))
  ll <- as.numeric(logLik(f))
  rec <- ref$ll_frm[ref$seed == seed]
  cat(sprintf("seed %d\n", seed))
  cat(sprintf("  help-page start vs the sweep's own, max |d| : %.3g\n", d))
  cat(sprintf("  logLik from it        : %.9f\n", ll))
  cat(sprintf("  the 0.2.2 sweep's own : %.9f\n", rec))
  cat(sprintf("  relative difference   : %.3g\n",
              abs(ll - rec) / abs(ll)))
}
