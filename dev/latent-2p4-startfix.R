# Lane `latent`, punch round 1: can `lca()`'s DEFAULT START be made to
# carry the 8 of 200, instead of a help-page recipe?
#
# The coordinator's option 1. Before choosing it, two things have to be
# measured, and the first is a diagnosis rather than a fix.
#
# WHY THE SHIPPED START FAILS ON THIS DESIGN. `lca_init_extras()` scores
# each subject by the MEAN of its item codes and cuts the scores into K
# equal-count slices. On the plan's design every class endorses exactly
# three of the ten items at 0.85 and the other seven at 0.2, so every
# class has the SAME expected mean score. The scoring statistic is blind
# to the classes here by construction, and the slices it produces are
# four samples of the same mixture. Section 1 below measures that.
#
# THE CANDIDATE. Slice on the response PATTERN instead of on a scalar
# score: one-hot the item codes, seed K centers by farthest-point (the
# row furthest from the column means, then repeatedly the row furthest
# from every center already chosen), run hard Lloyd iterations to a
# fixed point, and take each cluster's smoothed empirical proportions as
# that class's starting profile, exactly as today. It is deterministic,
# it costs O(n K J) and no fitting, and the clusters are ORDERED BY THE
# SAME MEAN SCORE afterwards so that the documented "class 1 is the
# low-score end" labeling survives wherever that score means anything.
#
# THE EVALUATION IS FREE OF NEW poLCA RUNS. `dev/latent-2p4-lca.tsv`
# already holds `ll_polca` for all 200 seeds from the 200-replicate
# sweep, so the reference is READ rather than recomputed; only the
# frmtmb fit is re-run, at 0.38 s median.
#
#   Rscript dev/latent-2p4-startfix.R [nseed]
#
# Seeds: the sweep's own, 20260910 + 0..199.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-lca-sim.R")

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args)) as.integer(args[[1L]]) else 200L
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

## ---- the two slicing rules, side by side ----------------------------

# today's: the mean item code, rank-cut into K equal-count slices
slice_score <- function(y, ncat, K) {
  n <- nrow(y)
  span <- pmax(ncat - 1L, 1L)
  sc <- rowMeans(sweep(y - 1, 2L, span, "/"), na.rm = TRUE)
  fin <- is.finite(sc)
  sc[!fin] <- if (any(fin)) stats::median(sc[fin]) else 0
  pmin(1L + ((rank(sc, ties.method = "first") - 1L) * K) %/% n, K)
}

# the candidate: hard k-means on the one-hot response pattern, seeded by
# farthest point, relabelled by the same mean score
slice_kmeans <- function(y, ncat, K, iter = 25L) {
  n <- nrow(y)
  J <- ncol(y)
  X <- matrix(0, n, sum(ncat))
  off <- 0L
  for (j in seq_len(J)) {
    v <- y[, j]
    md <- which.max(tabulate(v[!is.na(v)], nbins = ncat[j]))
    v[is.na(v)] <- md
    X[cbind(seq_len(n), off + v)] <- 1
    off <- off + ncat[j]
  }
  mu0 <- colMeans(X)
  d2 <- rowSums(sweep(X, 2L, mu0, "-")^2)
  idx <- integer(K)
  idx[1L] <- which.max(d2)
  if (K > 1L) {
    dmin <- rowSums(sweep(X, 2L, X[idx[1L], ], "-")^2)
    for (k in 2:K) {
      idx[k] <- which.max(dmin)
      dmin <- pmin(dmin, rowSums(sweep(X, 2L, X[idx[k], ], "-")^2))
    }
  }
  C <- X[idx, , drop = FALSE]
  cl <- rep(1L, n)
  for (it in seq_len(iter)) {
    dd <- vapply(seq_len(K),
                 function(k) rowSums(sweep(X, 2L, C[k, ], "-")^2),
                 numeric(n))
    new <- max.col(-dd, ties.method = "first")
    if (identical(new, cl) && it > 1L) break
    cl <- new
    for (k in seq_len(K)) {
      if (any(cl == k)) C[k, ] <- colMeans(X[cl == k, , drop = FALSE])
    }
  }
  # a cluster that emptied keeps its seed row, so every class starts
  # from a real response pattern rather than from a degenerate profile
  for (k in seq_len(K)) if (!any(cl == k)) cl[idx[k]] <- k
  # relabel by the same mean score the shipped rule uses, so the
  # documented low-to-high class ordering survives where it means
  # anything and is at least deterministic where it does not
  span <- pmax(ncat - 1L, 1L)
  sc <- rowMeans(sweep(y - 1, 2L, span, "/"), na.rm = TRUE)
  sc[!is.finite(sc)] <- 0
  ord <- order(vapply(seq_len(K), function(k) mean(sc[cl == k]), 0))
  match(cl, ord)
}

# the profile logits a slicing produces, the shipped arithmetic exactly,
# with an optional shrink toward the POOLED profile. `w = 0` is the
# shipped rule; a positive `w` keeps the classes separated but stops the
# start from being as confident as a hard partition makes it.
profiles_from <- function(y, ncat, K, slice, w = 0) {
  J <- ncol(y)
  out <- vector("list", J)
  for (j in seq_len(J)) {
    aj <- y[, j]
    aj <- aj[!is.na(aj)]
    pool <- (tabulate(aj, nbins = ncat[j]) + 1) /
      sum(tabulate(aj, nbins = ncat[j]) + 1)
    v <- numeric(0)
    for (k in seq_len(K)) {
      cj <- y[slice == k, j]
      cj <- cj[!is.na(cj)]
      cnt <- tabulate(cj, nbins = ncat[j]) + 1
      p <- cnt / sum(cnt)
      p <- (1 - w) * p + w * pool
      v <- c(v, log(p[-1L]) - log(p[1L]))
    }
    out[[j]] <- v
  }
  stats::setNames(out, paste0("pi", seq_len(J)))
}

## ---- 1. the diagnosis ------------------------------------------------

s <- lca_sim(seed = 20260970L)
Y <- s$Y
ncat <- rep(2L, ncol(Y))
cat("== 1. is the shipped scoring statistic informative here? ==\n")
span <- pmax(ncat - 1L, 1L)
sc <- rowMeans(sweep(Y - 1, 2L, span, "/"))
cat("mean item-code score, by TRUE class:\n")
print(round(tapply(sc, s$cl, mean), 4))
cat("its standard deviation within class:",
    formatC(round(tapply(sc, s$cl, stats::sd), 4), format = "g"), "\n")
ss <- stats::anova(stats::lm(sc ~ factor(s$cl)))
cat("between-class share of the score's variance (R squared):",
    formatC(ss[1, 2] / sum(ss[, 2]), digits = 4, format = "g"), "\n")
cat("\nshipped slice against the true class:\n")
print(table(slice = slice_score(Y, ncat, K), true = s$cl))
cat("\ncandidate slice against the true class:\n")
print(table(slice = slice_kmeans(Y, ncat, K), true = s$cl))

## ---- 2. all 200 seeds, reference read from the sweep -----------------

ref <- utils::read.delim("dev/latent-2p4-lca.tsv", check.names = FALSE)
out <- "dev/latent-2p4-startfix.tsv"
cat("seed\tll_polca\tll_shipped\tll_km\tll_km05\tgap_shipped",
    "\tgap_km\tgap_km05\tseconds\n", sep = "", file = out)
cat("\n== 2. the candidates on every seed of the sweep ==\n")
fit_from <- function(dd, st) {
  f <- try(suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = dd,
        control = tight, start = st))), silent = TRUE)
  if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
}
for (i in seq_len(min(NSEED, nrow(ref)))) {
  seed <- ref$seed[i]
  sm <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(sm$Y))
  sl <- slice_kmeans(sm$Y, nc, K)
  t0 <- Sys.time()
  ll_km <- fit_from(sm$dd, profiles_from(sm$Y, nc, K, sl, w = 0))
  ll_k5 <- fit_from(sm$dd, profiles_from(sm$Y, nc, K, sl, w = 0.5))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(paste(c(seed,
              formatC(c(ref$ll_polca[i], ref$ll_frm[i], ll_km, ll_k5),
                      digits = 6, format = "f"),
              formatC(c(ref$ll_frm[i] - ref$ll_polca[i],
                        ll_km - ref$ll_polca[i],
                        ll_k5 - ref$ll_polca[i], secs),
                      digits = 6, format = "g")),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  if (i %% 25L == 0L) {
    cat("  ", i, " seeds done\n", sep = "")
    flush(stdout())
  }
}

d <- utils::read.delim(out)
tol <- 1e-6 * abs(d$ll_polca)
lost <- function(g) sum(g < -tol, na.rm = TRUE)
cat("\nseeds:", nrow(d), "\n")
cat("  shipped start below poLCA           :", lost(d$gap_shipped), "\n")
cat("  kmeans start below poLCA            :", lost(d$gap_km), "\n")
cat("  kmeans shrunk 0.5 below poLCA       :", lost(d$gap_km05), "\n")
best <- pmax(d$gap_km, d$gap_km05, na.rm = TRUE)
cat("  better of the two kmeans starts     :", lost(best), "\n")
best3 <- pmax(d$gap_shipped, d$gap_km, d$gap_km05, na.rm = TRUE)
cat("  better of all three starts          :", lost(best3), "\n")
cat("  any start ABOVE poLCA by 1e-6 rel   :",
    sum(best3 > tol, na.rm = TRUE), "\n")
cat("  errors: km", sum(is.na(d$ll_km)), " km05", sum(is.na(d$ll_km05)),
    "\n")
bad_old <- d$seed[d$gap_shipped < -tol]
cat("\nthe eight the shipped start loses:\n")
print(d[d$seed %in% bad_old,
        c("seed", "gap_shipped", "gap_km", "gap_km05")],
      row.names = FALSE)
for (nm in c("gap_km", "gap_km05")) {
  bn <- d$seed[!is.na(d[[nm]]) & d[[nm]] < -tol & !(d$seed %in% bad_old)]
  cat("\n", nm, " loses ", length(bn),
      " seeds the shipped start does NOT:\n", sep = "")
  if (length(bn)) {
    print(d[d$seed %in% bn,
            c("seed", "gap_shipped", "gap_km", "gap_km05")],
          row.names = FALSE)
  }
}
cat("\ntwo frm fits with the candidate starts, seconds: min",
    formatC(min(d$seconds), digits = 3, format = "g"), " median",
    formatC(stats::median(d$seconds), digits = 3, format = "g"),
    " max", formatC(max(d$seconds), digits = 3, format = "g"), "\n")
