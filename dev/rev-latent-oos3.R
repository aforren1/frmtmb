# Reviewer, punch round 1, attack 1 concluded: adjudicate the three
# out-of-sample seeds where a start lost, against poLCA, and pool the
# weight grid over all 400 seeds.
#
# `dev/rev-latent-oos2.R` screened 160 fresh seeds against the best of
# eight starts rather than against poLCA, which is enough to say a
# start was BEATEN but not enough to say by how much against the
# reference the help page quotes. Three seeds came out of that screen:
#
#   20270476  the SHIPPED start lost by 256 units
#   20270506  the OLD start lost by 256 units
#   20270591  the OLD start lost by 236 units
#
# Each gets a poLCA(nrep = 10) arm here, and the shipped fit's own
# diagnostics, because "silently" is half the claim.
#
#   Rscript dev/rev-latent-oos3.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
source(paste0(D, "latent-2p4-startfix-fns.R"))
suppressMessages(loadNamespace("poLCA"))
K <- 4L
WS <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9)
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))

for (seed in c(20270476L, 20270506L, 20270591L)) {
  cat("\n======== seed", seed, "========\n")
  s <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(s$Y))
  dp <- data.frame(x1 = s$dd$x1, x2 = s$dd$x2)
  for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
  ff <- stats::as.formula(paste0("cbind(",
    paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "),
    ") ~ x1 + x2"))
  set.seed(920000L + (seed %% 1000L))
  pl <- poLCA::poLCA(ff, dp, nclass = K, nrep = 10L, verbose = FALSE,
                     maxiter = 20000, tol = 1e-12)
  fs <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
        control = tight)))
  ll <- as.numeric(logLik(fs))
  dg <- frmtmb::diagnose(fs, quiet = TRUE)
  msg <- utils::capture.output(frmtmb::diagnose(fs))
  cat("poLCA(nrep = 10)     :", format(pl$llik, digits = 12), "\n")
  cat("SHIPPED start        :", format(ll, digits = 12),
      sprintf("  gap %+.3f\n", ll - pl$llik))
  ll_old <- {
    f <- suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight,
          start = profiles_from(s$Y, nc, K, slice_score(s$Y, nc, K),
                                w = 0))))
    as.numeric(logLik(f))
  }
  cat("OLD score cut        :", format(ll_old, digits = 12),
      sprintf("  gap %+.3f\n", ll_old - pl$llik))
  sl <- slice_kmeans(s$Y, nc, K)
  for (w in WS) {
    f <- suppressWarnings(suppressMessages(
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
          control = tight, start = profiles_from(s$Y, nc, K, sl, w = w))))
    cat(sprintf("  w = %-5s          : %s  gap %+.3f\n", format(w),
                format(as.numeric(logLik(f)), digits = 12),
                as.numeric(logLik(f)) - pl$llik))
  }
  cat("the SHIPPED fit's own diagnostics:\n")
  cat("  nlminb", dg$convergence, " pdHess", isTRUE(dg$pdHess),
      " bad_se", length(dg$bad_se), " flat", length(dg$flat),
      " max_grad/|ll|", format(dg$max_grad / abs(ll), digits = 3), "\n")
  cat(paste0("  | ", msg, collapse = "\n"), "\n")
  cat("  prints 'No convergence problems detected':",
      any(grepl("No convergence problems detected", msg)), "\n")
}

## --------------------------------------- the grid, pooled over 400
cat("\n\n======== the weight grid over ALL 400 seeds ========\n")
cols <- paste0("ll_w", sub("[.]", "", WS))
ins <- utils::read.delim(paste0(D, "latent-2p4-shrink.tsv"))
gcols <- paste0("gap_w", sub("[.]", "", WS))
oosa <- utils::read.delim(paste0(D, "rev-latent-oos.tsv"))
oosb <- utils::read.delim(paste0(D, "rev-latent-oos2.tsv"))
oos <- rbind(oosa[, c("seed", "ll_ship", "ll_old", cols)],
             oosb[, c("seed", "ll_ship", "ll_old", cols)])
Mo <- as.matrix(oos[, c("ll_ship", "ll_old", cols)])
besto <- apply(Mo, 1L, max, na.rm = TRUE)
tolo <- 1e-6 * abs(besto)
toli <- 1e-6 * abs(ins$ll_polca)
cat(sprintf("  %-8s %14s %16s %14s\n", "w", "in-sample lost",
            "out-of-sample lost", "pooled lost"))
for (k in seq_along(WS)) {
  li <- sum(ins[[gcols[k]]] < -toli, na.rm = TRUE)
  lo <- sum((Mo[, cols[k]] - besto) < -tolo, na.rm = TRUE)
  cat(sprintf("  %-8s %14s %16s %14s\n", format(WS[k]),
              paste0(li, " / 200"), paste0(lo, " / 200"),
              paste0(li + lo, " / 400")))
}
li_old <- 8L
lo_old <- sum((Mo[, "ll_old"] - besto) < -tolo, na.rm = TRUE)
cat(sprintf("  %-8s %14s %16s %14s\n", "OLD",
            paste0(li_old, " / 200"), paste0(lo_old, " / 200"),
            paste0(li_old + lo_old, " / 400")))
li_sh <- sum(ins[[gcols[4L]]] < -toli, na.rm = TRUE)
lo_sh <- sum((Mo[, "ll_ship"] - besto) < -tolo, na.rm = TRUE)
cat(sprintf("  %-8s %14s %16s %14s   <- what ships\n", "SHIPPED",
            paste0(li_sh, " / 200"), paste0(lo_sh, " / 200"),
            paste0(li_sh + lo_sh, " / 400")))
