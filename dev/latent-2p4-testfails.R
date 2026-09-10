# Lane `latent`, item 2.4: see the new identity test FAIL.
#
# test-lca.R's realistic-scale block asserts that lca() and poLCA reach
# the same optimum, with every tolerance a ratio to poLCA's own standard
# errors. A green assertion is only evidence if it can go red, so this
# runs the same block's body on seed 20260970, one of the eight of 200
# where the deterministic start reaches a local optimum, and reports
# each expectation's verdict.
#
#   Rscript dev/latent-2p4-testfails.R [seed]

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.latent)
})
suppressMessages(loadNamespace("poLCA"))
source("dev/latent-lca-sim.R")

args <- commandArgs(trailingOnly = TRUE)
SEED <- if (length(args)) as.integer(args[[1L]]) else 20260970L
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
s <- lca_sim(seed = SEED)
J <- ncol(s$Y)
dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
for (j in seq_len(J)) dp[[paste0("I", j)]] <- s$Y[, j]

fit <- suppressWarnings(suppressMessages(
  frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
      control = tight)))
set.seed(700001L)
pl <- poLCA::poLCA(
  stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(J)), collapse = ", "), ") ~ x1 + x2")),
  dp, nclass = K, nrep = 10L, verbose = FALSE, maxiter = 20000,
  tol = 1e-12)

ll <- as.numeric(logLik(fit))
Pf <- vapply(lca_profiles(fit), function(m) m[, 2L], numeric(K))
Pp <- vapply(pl$probs, function(m) m[, 2L], numeric(K))
perms <- as.matrix(expand.grid(rep(list(seq_len(K)), K)))
perms <- perms[apply(perms, 1L, function(p) length(unique(p)) == K), ]
dst <- apply(perms, 1L, function(p) max(abs(Pf[p, ] - Pp)))
se_pr <- max(vapply(pl$probs.se, function(m) m[, 2L], numeric(K)))

report <- function(label, value, bound) {
  cat(sprintf("  %-34s %12.4g  bound %.0e  %s\n", label, value, bound,
              if (value < bound) "PASSES" else "FAILS"))
}
cat("seed", SEED, "\n")
cat("  frm logLik   ", format(ll, digits = 10), "\n")
cat("  poLCA logLik ", format(pl$llik, digits = 10), "\n")
cat("  gap          ", format(ll - pl$llik, digits = 6), "\n\n")
report("|dll| / |ll|", abs(ll - pl$llik) / abs(ll), 1e-10)
report("profile gap / poLCA profile SE", min(dst) / se_pr, 1e-3)
