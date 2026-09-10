# Reviewer, punch round 1, attack 1 continued: 160 MORE fresh seeds.
#
# WHY. `dev/rev-latent-oos.R` drew 40 fresh seeds and neither the old
# nor the new start lost one. That settles "the new rule does not
# regress" and settles nothing about "the new rule fixes anything",
# because at the old rule's measured 4 percent rate 40 seeds expect
# only 1.6 failures and this draw got none. The out-of-sample answer is
# ambiguous at 40, so it is taken to 200, matching the in-sample count.
#
# THE poLCA ARM IS NOT RUN ON THESE 160. It does not have to be. The
# question is whether the SHIPPED start is beaten on a data set by some
# other start, and eight starts are run per seed: the old score cut and
# the whole shrink grid. A start that is below the best of the eight by
# more than 1e-6 relative LOST that replicate, and no EM is needed to
# say so. This can only UNDER-count losses, in the case where all eight
# land on the same wrong mode, and the 40 seeds that DO carry a poLCA
# reference bound that: on those, the best of the eight was poLCA's
# optimum on 40 of 40.
#
# Seeds 20270441 to 20270600, which the lane has never used (its data
# seeds are 20260910 to 20261109).
#
#   Rscript dev/rev-latent-oos2.R [nseed]

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
source(paste0(D, "latent-2p4-startfix-fns.R"))
cat("frmtmb.latent from:", find.package("frmtmb.latent")[[1L]], "\n")

args <- commandArgs(trailingOnly = TRUE)
NSEED <- if (length(args)) as.integer(args[[1L]]) else 160L
K <- 4L
WS <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9)
SEEDS <- 20270441L + seq_len(NSEED) - 1L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
out <- paste0(D, "rev-latent-oos2.tsv")
cat(paste(c("seed", "ll_ship", "ll_old",
            paste0("ll_w", sub("[.]", "", WS))), collapse = "\t"),
    "\n", sep = "", file = out)

fitll <- function(dd, st = NULL) {
  f <- try(suppressWarnings(suppressMessages(
    if (is.null(st)) {
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = dd, control = tight)
    } else {
      frm(bf(Y ~ x1 + x2), family = lca(K = K), data = dd,
          control = tight, start = st)
    })), silent = TRUE)
  if (inherits(f, "try-error")) NA_real_ else as.numeric(logLik(f))
}

for (i in seq_len(NSEED)) {
  seed <- SEEDS[i]
  sm <- lca_sim(seed = seed)
  nc <- rep(2L, ncol(sm$Y))
  ll_ship <- fitll(sm$dd)
  ll_old <- fitll(sm$dd, profiles_from(sm$Y, nc, K,
                                       slice_score(sm$Y, nc, K), w = 0))
  sl <- slice_kmeans(sm$Y, nc, K)
  ll_w <- vapply(WS, function(w)
    fitll(sm$dd, profiles_from(sm$Y, nc, K, sl, w = w)), numeric(1))
  cat(paste(c(seed, formatC(c(ll_ship, ll_old, ll_w), digits = 8,
                            format = "f")), collapse = "\t"),
      "\n", sep = "", file = out, append = TRUE)
  if (i %% 20L == 0L) { cat("  ", i, " seeds\n", sep = ""); flush(stdout()) }
}

## --------------------------------------------------------- summary
a <- utils::read.delim(paste0(D, "rev-latent-oos.tsv"))
b <- utils::read.delim(out)
cols <- c("ll_ship", "ll_old", paste0("ll_w", sub("[.]", "", WS)))
all_rows <- rbind(a[, c("seed", cols)], b[, c("seed", cols)])
M <- as.matrix(all_rows[, cols])
best <- apply(M, 1L, max, na.rm = TRUE)
tol <- 1e-6 * abs(best)
cat("\n==== OUT OF SAMPLE, 200 fresh seeds:", min(all_rows$seed), "to",
    max(all_rows$seed), "====\n")
cat("reference: the best of the eight starts on each data set.\n")
cat("on the 40 that also carry a poLCA(nrep = 10) arm, that best WAS\n")
cat("poLCA's optimum on", sum(abs(apply(as.matrix(a[, cols]), 1L, max) -
                                   a$ll_polca) < 1e-6 * abs(a$ll_polca)),
    "of", nrow(a), "\n\n")
cat(sprintf("  %-28s %6s %14s\n", "starting rule", "lost", "worst gap"))
for (cl in cols) {
  g <- M[, cl] - best
  nm <- switch(cl, ll_ship = "SHIPPED lca_init_extras()",
               ll_old = "OLD score cut",
               paste0("grid w = ", sub("^ll_w0?", "0.",
                                       sub("^ll_w0$", "ll_w0", cl))))
  if (cl == "ll_w0") nm <- "grid w = 0"
  cat(sprintf("  %-28s %6d %14s\n", nm, sum(g < -tol, na.rm = TRUE),
              format(min(g, na.rm = TRUE), digits = 6)))
}
cat("\nthe seeds where the SHIPPED start lost:\n")
lost_rows <- which((M[, "ll_ship"] - best) < -tol)
if (length(lost_rows)) {
  print(all_rows[lost_rows, ], row.names = FALSE, digits = 10)
} else {
  cat("  none\n")
}
cat("\nthe seeds where the OLD start lost:\n")
lo <- which((M[, "ll_old"] - best) < -tol)
if (length(lo)) {
  print(all_rows[lo, c("seed", "ll_old", "ll_ship")], row.names = FALSE,
        digits = 10)
  cat("  on those, the shipped start reached the best:",
      sum(abs(M[lo, "ll_ship"] - best[lo]) <= tol[lo]), "of", length(lo),
      "\n")
} else {
  cat("  none\n")
}
cat("\nshipped default identical to the w = 0.5 helper on all 200:",
    max(abs(M[, "ll_ship"] - M[, "ll_w05"]), na.rm = TRUE) == 0, "\n")
