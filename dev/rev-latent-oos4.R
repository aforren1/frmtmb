# Reviewer: before recommending a different shrink weight, probe the
# right-hand end of the grid, which the lane's sweep stops short of.
#
# WHY. Pooled over 400 seeds, `w = 0.9` loses 0 where the shipped
# `w = 0.5` loses 1. But `w = 1` makes every class profile the POOLED
# profile, which is the symmetry axis `?hmm` warns about for mixtures:
# a start with every component equal is a start the optimizer cannot
# leave. So 0.9 sits next to a known degeneracy and a recommendation to
# move there is worth nothing until the degeneracy is located.
#
# Weights 0.9, 0.95, 0.99 and 1.0, on:
#   * the 17 seeds where SOMETHING is known to go wrong (the 8 the old
#     rule lost, the 6 `w = 0` lost, the 3 out-of-sample failures);
#   * the 40 fresh seeds that carry a poLCA(nrep = 10) reference.
#
#   Rscript dev/rev-latent-oos4.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
source(paste0(D, "latent-2p4-startfix-fns.R"))
K <- 4L
WS <- c(0.5, 0.9, 0.95, 0.99, 1)
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

hard <- c(20260970L, 20260990L, 20260996L, 20260999L, 20261010L,
          20261013L, 20261042L, 20261043L,          # the old rule's 8
          20260911L, 20260953L, 20261009L, 20261071L, 20261100L,
          20261109L,                                # w = 0's 6
          20270476L, 20270506L, 20270591L)          # out of sample
ref_in <- utils::read.delim(paste0(D, "latent-2p4-lca.tsv"))
ref_oos <- utils::read.delim(paste0(D, "rev-latent-oos.tsv"))
polca_of <- function(seed) {
  i <- match(seed, ref_in$seed)
  if (!is.na(i)) return(ref_in$ll_polca[i])
  j <- match(seed, ref_oos$seed)
  if (!is.na(j)) return(ref_oos$ll_polca[j])
  NA_real_
}
extra_polca <- c(`20270476` = -11507.6035576, `20270506` = -11273.6911314,
                 `20270591` = -11437.8232163)

run <- function(seeds, label) {
  cat("\n======== ", label, ": ", length(seeds), " seeds ========\n",
      sep = "")
  M <- matrix(NA_real_, length(seeds), length(WS),
              dimnames = list(as.character(seeds), format(WS)))
  pv <- numeric(length(seeds))
  for (i in seq_along(seeds)) {
    seed <- seeds[i]
    s <- lca_sim(seed = seed)
    nc <- rep(2L, ncol(s$Y))
    sl <- slice_kmeans(s$Y, nc, K)
    pv[i] <- polca_of(seed)
    if (is.na(pv[i]) && as.character(seed) %in% names(extra_polca)) {
      pv[i] <- extra_polca[[as.character(seed)]]
    }
    for (k in seq_along(WS)) {
      f <- try(suppressWarnings(suppressMessages(
        frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
            control = tight,
            start = profiles_from(s$Y, nc, K, sl, w = WS[k])))),
        silent = TRUE)
      M[i, k] <- if (inherits(f, "try-error")) NA_real_ else
        as.numeric(logLik(f))
    }
  }
  ref <- ifelse(is.na(pv), apply(M, 1L, max, na.rm = TRUE), pv)
  tol <- 1e-6 * abs(ref)
  G <- sweep(M, 1L, ref, "-")
  cat("gaps to the reference (poLCA where available):\n")
  print(round(G, 4))
  cat("\nlost, by weight:\n")
  for (k in seq_along(WS)) {
    cat(sprintf("  w = %-5s lost %2d of %d   worst %10s   errored %d\n",
                format(WS[k]), sum(G[, k] < -tol, na.rm = TRUE),
                nrow(G), format(min(G[, k], na.rm = TRUE), digits = 5),
                sum(is.na(M[, k]))))
  }
  invisible(G)
}

run(hard, "the 17 seeds where something is known to go wrong")
run(ref_oos$seed, "the 40 fresh seeds with a poLCA reference")

cat("\n======== is w = 1 the symmetry axis? ========\n")
s <- lca_sim(seed = 20260910L)
nc <- rep(2L, ncol(s$Y))
sl <- slice_kmeans(s$Y, nc, K)
for (w in c(0.9, 0.99, 1)) {
  st <- profiles_from(s$Y, nc, K, sl, w = w)
  m <- do.call(rbind, lapply(st, function(v) v))
  cat(sprintf("w = %-5s : largest spread across the K classes in any",
              format(w)),
      "item logit:", format(max(apply(matrix(unlist(st), ncol = K,
                                             byrow = TRUE), 1L,
                                      function(z) diff(range(z)))),
                            digits = 4), "\n")
}
