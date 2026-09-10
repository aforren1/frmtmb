# Reviewer, punch round 1, attack 1: the new `lca()` start was TUNED and
# EVALUATED on the same 200 seeds. Is the shrink weight an overfit
# constant?
#
# The lane swept `w` over 0, 0.1, 0.25, 0.5, 0.75, 0.9 on seeds
# 20260910 + 0..199, found losses of 6, 5, 2, 0, 1, 0, took 0.5 as the
# smallest zero-loss weight, and then reported "200 of 200" on those
# same 200 seeds. A weight picked as the unique zero-loss point on the
# data it was picked from is the shape of an overfit constant, and the
# lane's own note that the plateau is not flat (0.75 loses one where
# 0.5 and 0.9 lose none) sharpens the worry rather than settling it.
#
# THE SEEDS HERE ARE NEW. The lane's data seeds are 20260910 to
# 20261109 and its poLCA seeds are 700001 to 700200, 800001 to 800010,
# 600xxx and 500xxx. This draws 20270401 to 20270440 with poLCA seeds
# 910001 to 910040, at the identical design, and measures:
#
#   * the SHIPPED start (lca_init_extras(), what the package now does);
#   * the OLD start (score cut, no shrink), so the improvement is
#     measured out of sample and not only the level;
#   * the whole weight grid, so "0.5 is the smallest weight that loses
#     none" can be checked on seeds it was not chosen from.
#
# The poLCA arm IS run here: an out-of-sample check cannot read its
# reference from the lane's in-sample file.
#
#   Rscript dev/rev-latent-oos.R [nseed]

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
source(paste0(D, "latent-2p4-startfix-fns.R"))
suppressMessages(loadNamespace("poLCA"))
cat("frmtmb.latent from:", find.package("frmtmb.latent")[[1L]], "\n")

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args)) as.integer(args[[1L]]) else 40L
K <- 4L
WS <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9)
SEEDS <- 20270401L + seq_len(NSEED) - 1L
PSEEDS <- 910001L + seq_len(NSEED) - 1L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

out <- paste0(D, "rev-latent-oos.tsv")
cat(paste(c("seed", "seed_polca", "ll_polca", "ll_ship", "ll_old",
            paste0("ll_w", sub("[.]", "", WS)), "s_ship", "s_polca",
            "pdhess", "maxgrad_rel",
            paste0("cov_", as.vector(t(outer(2:4, 1:3, function(i, j)
              paste0("c", i, "_", j)))))),
          collapse = "\t"), "\n", sep = "", file = out)

gate_cover <- function(fit, K) {
  # the nine gating coefficients, re-referenced to class 1 the way the
  # lane's own harness does, with the Wald interval built on the
  # contrast the re-referencing implies
  b <- fixef(fit)
  cn <- names(b[[1L]])
  V <- vcov(fit)
  nms <- rownames(V)
  truth <- rbind(c(-0.4, 0.8, -0.5), c(0.2, -0.6, 0.9),
                 c(-0.1, 0.3, 0.4))
  covs <- numeric(0)
  for (k in 2:K) {
    for (j in seq_along(cn)) {
      # theta_{k-1} minus theta_1, in frmtmb's last-class reference
      est <- unname(b[[paste0("theta", k - 1L)]][j]) -
        unname(b[["theta1"]][j])
      L <- rep(0, length(nms))
      i1 <- which(nms == paste0("theta", k - 1L, "_", cn[j]))
      i2 <- which(nms == paste0("theta1_", cn[j]))
      L[i1] <- L[i1] + 1
      L[i2] <- L[i2] - 1
      se <- sqrt(drop(crossprod(L, V %*% L)))
      covs <- c(covs, as.integer(
        abs(est - truth[k - 1L, j]) <= stats::qnorm(0.975) * se))
    }
  }
  covs
}

for (i in seq_len(NSEED)) {
  seed <- SEEDS[i]
  sm <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(sm$Y))
  dp <- data.frame(x1 = sm$dd$x1, x2 = sm$dd$x2)
  for (j in seq_len(ncol(sm$Y))) dp[[paste0("I", j)]] <- sm$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(sm$Y))), collapse = ", "),
    ") ~ x1 + x2"))

  t0 <- Sys.time()
  fs <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = sm$dd,
        control = tight)))
  s_ship <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ll_ship <- as.numeric(logLik(fs))
  dg <- frmtmb::diagnose(fs, quiet = TRUE)

  t1 <- Sys.time()
  set.seed(PSEEDS[i])
  pl <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                     maxiter = 20000, tol = 1e-12)
  s_polca <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  ll_old <- {
    f <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = sm$dd,
          control = tight,
          start = profiles_from(sm$Y, nc, K,
                                slice_score(sm$Y, nc, K), w = 0)))),
      silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
  }

  sl <- slice_kmeans(sm$Y, nc, K)
  ll_w <- vapply(WS, function(w) {
    f <- try(suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = sm$dd,
          control = tight,
          start = profiles_from(sm$Y, nc, K, sl, w = w)))), silent = TRUE)
    if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
  }, numeric(1))

  cat(paste(c(seed, PSEEDS[i],
              formatC(c(pl$llik, ll_ship, ll_old, ll_w), digits = 8,
                      format = "f"),
              formatC(c(s_ship, s_polca), digits = 4, format = "g"),
              isTRUE(dg$pdHess),
              formatC(dg$max_grad / abs(ll_ship), digits = 4,
                      format = "e"),
              gate_cover(fs, K)),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  cat(sprintf("%3d seed %d: polca %.4f  ship %+.4f  old %+.4f\n",
              i, seed, pl$llik, ll_ship - pl$llik, ll_old - pl$llik))
  flush(stdout())
}

## ------------------------------------------------------------ summary
d <- utils::read.delim(out)
tol <- 1e-6 * abs(d$ll_polca)
lost <- function(v) sum(v < d$ll_polca - tol, na.rm = TRUE)
cat("\n==== OUT OF SAMPLE, seeds", min(SEEDS), "to", max(SEEDS),
    "====\n")
cat("replicates:", nrow(d), "\n")
cat(sprintf("  %-28s %6s %8s\n", "starting rule", "lost", "reaches"))
cat(sprintf("  %-28s %6d %8s\n", "OLD score cut", lost(d$ll_old),
            paste0(nrow(d) - lost(d$ll_old), " of ", nrow(d))))
cat(sprintf("  %-28s %6d %8s\n", "SHIPPED lca_init_extras()",
            lost(d$ll_ship),
            paste0(nrow(d) - lost(d$ll_ship), " of ", nrow(d))))
cat("\n  the weight grid, out of sample:\n")
for (k in seq_along(WS)) {
  v <- d[[paste0("ll_w", sub("[.]", "", WS[k]))]]
  cat(sprintf("  %-6s lost %2d   worst gap %10s   errored %d\n",
              format(WS[k]), lost(v),
              format(min(v - d$ll_polca, na.rm = TRUE), digits = 5),
              sum(is.na(v))))
}
cat("\n  does the shipped default match the w = 0.5 helper?\n")
cat("    max |ll_ship - ll_w05| :",
    format(max(abs(d$ll_ship - d$ll_w05), na.rm = TRUE), digits = 4),
    "\n")
cat("\n  where the two rules differ (out of sample):\n")
diff_rows <- which(abs(d$ll_ship - d$ll_old) > tol)
if (length(diff_rows)) {
  print(d[diff_rows, c("seed", "ll_polca", "ll_ship", "ll_old")],
        row.names = FALSE, digits = 10)
} else {
  cat("    none\n")
}
cat("\n  pdHess:", sum(d$pdhess), "of", nrow(d),
    "  maxgrad_rel range:", format(range(d$maxgrad_rel), digits = 3),
    "\n")
cc <- as.matrix(d[, grep("^cov_", names(d))])
p <- mean(cc)
se <- sqrt(p * (1 - p) / length(cc))
cat("  Wald coverage:", sprintf("%.4f (%d/%d) MC %.4f to %.4f",
                                p, sum(cc), length(cc),
                                p - 1.96 * se, p + 1.96 * se), "\n")
cat("  median fit seconds, shipped:",
    format(stats::median(d$s_ship), digits = 4), "\n")
