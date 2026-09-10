# Reviewer, punch round 1, attacks 1b and 4: recompute the re-measured
# recovery table from raw columns, and check the diagnosis of the rule
# that was replaced.
#
#  1b  the old score cut is claimed to explain 0.000814 of its own
#      variance between the true classes on this design, because every
#      class endorses three of ten items at 0.85 and seven at 0.2. If
#      that is right the old rule was blind by construction. Computed
#      here from the simulator, on the lane's seed and on others, with
#      no fitting.
#
#   4  the re-measured table: 200 of 200, 1.8e-14 relative, pdHess 200
#      of 200, coverage 94.3 percent at 1698 of 1800. Recomputed from
#      dev/latent-2p4-recheck.tsv the way the first table was, and the
#      denominator change accounted for.
#
#   Rscript dev/rev-latent-arith2.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env2.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
rd <- function(f) utils::read.delim(paste0(D, f), stringsAsFactors = FALSE)
hr <- function(s) cat("\n==== ", s, " ====\n", sep = "")
mci <- function(k, n) {
  p <- k / n
  se <- sqrt(p * (1 - p) / n)
  sprintf("%.4f (%d/%d)  MC 95%%: %.4f to %.4f", p, k, n,
          p - 1.96 * se, p + 1.96 * se)
}

## ------------------------------------------------- 1b: the diagnosis
hr("1b. is the old score blind on this design, by construction?")
base <- lca_truth_base(4L, 10L)
cat("the truth's item table, P(endorse) by class and item:\n")
print(round(base, 2))
cat("\nrow sums (the expected mean item code, minus 1, per class):\n")
print(round(rowMeans(base), 6))
cat("every class endorses exactly", sum(base[1L, ] > 0.5), "of",
    ncol(base), "items at 0.85 and the rest at 0.2, so the EXPECTED\n")
cat("score is identical across classes:",
    length(unique(round(rowMeans(base), 12))) == 1L, "\n")

for (seed in c(20260970L, 20260910L, 20270401L)) {
  s <- lca_sim(seed = seed)
  span <- rep(1, ncol(s$Y))
  sc <- rowMeans(sweep(s$Y - 1, 2L, span, "/"))
  ss_tot <- sum((sc - mean(sc))^2)
  ss_bet <- sum(vapply(sort(unique(s$cl)), function(k) {
    sum(s$cl == k) * (mean(sc[s$cl == k]) - mean(sc))^2
  }, numeric(1)))
  cat(sprintf("\nseed %d: R2 of the score between the TRUE classes = %.6f\n",
              seed, ss_bet / ss_tot))
  cat("  class means:",
      paste(format(vapply(sort(unique(s$cl)),
                          function(k) mean(sc[s$cl == k]), numeric(1)),
                   digits = 4), collapse = " "), "\n")
  cat("  within-class sd:",
      format(sqrt(mean(vapply(sort(unique(s$cl)), function(k)
        stats::var(sc[s$cl == k]), numeric(1)))), digits = 4),
      "  (lane: 0.127 on seed 20260970)\n")
  # what the slices the old rule produced actually contained
  n <- length(sc)
  sl <- pmin(1L + ((rank(sc, ties.method = "first") - 1L) * 4L) %/% n, 4L)
  tb <- table(slice = sl, truth = s$cl)
  cat("  the four slices' composition by true class:\n")
  print(round(prop.table(tb, 1L), 3))
}

## ----------------------------------------------- 4: the new table
hr("4. the re-measured recovery table")
r <- rd("latent-2p4-recheck.tsv")
cat("columns:", paste(names(r), collapse = " "), "\n")
cat("rows:", nrow(r), " distinct seeds:", length(unique(r$seed)),
    " contiguous 20260910..20261109:",
    identical(sort(r$seed), 20260910:20261109), "\n")
num <- function(v) suppressWarnings(as.numeric(trimws(as.character(v))))
for (nm in names(r)) if (!is.numeric(r[[nm]]) && nm != "pdhess") {
  z <- num(r[[nm]]); if (!all(is.na(z))) r[[nm]] <- z
}
cat("\nreaches poLCA's optimum (|d|/|ll| < 1e-6):",
    sum(abs(r$ll_rel_gap) < 1e-6, na.rm = TRUE), "of", nrow(r),
    "  claim 200 of 200\n")
cat("worst |d|/|ll|:", format(max(abs(r$ll_rel_gap), na.rm = TRUE),
                              digits = 3), " claim 1.8e-14\n")
pd <- r$pdhess
cat("pdHess TRUE:", sum(pd %in% c(TRUE, "TRUE")), "of", nrow(r),
    " claim 200 of 200\n")
cat("maxgrad_rel range:", format(range(num(r$maxgrad_rel), na.rm = TRUE),
                                 digits = 3),
    " claim 7.5e-10 to 6.5e-09\n")

cov_cols <- grep("^cov_", names(r), value = TRUE)
cat("\ncoverage columns:", length(cov_cols), ":",
    paste(cov_cols, collapse = " "), "\n")
m <- as.matrix(r[, cov_cols])
m <- matrix(num(m), nrow = nrow(r))
cat("NA cells:", sum(is.na(m)), "  non-0/1 cells:",
    sum(!is.na(m) & !(m %in% c(0, 1))), "\n")
cat("overall:", mci(sum(m, na.rm = TRUE), sum(!is.na(m))),
    " claim 94.3 percent, 1698 of 1800, MC 93.3 to 95.4\n")
cat("\nDENOMINATOR: 200 replicates x 9 coefficients = 1800, against the\n")
cat("first table's 192 x 9 = 1728. The 8 that were excluded as a\n")
cat("different optimum are now IN, because the start reaches the same\n")
cat("optimum on all 200. 1800 - 1728 =", 1800 - 1728, "= 8 x 9.\n")

if ("s_frm" %in% names(r)) {
  cat("\nmedian fit seconds:", format(stats::median(num(r$s_frm),
                                                    na.rm = TRUE),
                                      digits = 4),
      " claim 0.341 (old start 0.376)\n")
}

## the per-coefficient table
truth <- c(c2_1 = -0.4, c2_2 = 0.8, c2_3 = -0.5,
           c3_1 = 0.2, c3_2 = -0.6, c3_3 = 0.9,
           c4_1 = -0.1, c4_2 = 0.3, c4_3 = 0.4)
if (all(paste0("err_", names(truth)) %in% names(r))) {
  tab <- do.call(rbind, lapply(names(truth), function(nm) {
    e <- num(r[[paste0("err_", nm)]]); s <- num(r[[paste0("se_", nm)]])
    data.frame(coef = nm, truth = truth[[nm]], bias = mean(e),
               rmse = sqrt(mean(e^2)), empSD = stats::sd(e),
               meanSE = mean(s),
               cover95 = 100 * mean(abs(e) <= stats::qnorm(0.975) * s),
               nNA = sum(is.na(e) | is.na(s)))
  }))
  cat("\nper-coefficient, recomputed:\n")
  print(tab, row.names = FALSE, digits = 4)
  cat("sum of per-coef covered:",
      sum(round(tab$cover95 / 100 * nrow(r))), " vs matrix total ",
      sum(m, na.rm = TRUE), "\n")
}

## -------------------------------------- the shrink sweep, recomputed
hr("the in-sample weight sweep, recomputed")
sh <- rd("latent-2p4-shrink.tsv")
cat("rows:", nrow(sh), " seeds contiguous:",
    identical(sort(sh$seed), 20260910:20261109), "\n")
tol <- 1e-6 * abs(sh$ll_polca)
ws <- c("0", "01", "025", "05", "075", "09")
lab <- c("0", "0.1", "0.25", "0.5", "0.75", "0.9")
for (k in seq_along(ws)) {
  g <- num(sh[[paste0("gap_w", ws[k])]])
  cat(sprintf("  w %-5s lost %2d  errored %d  worst gap %10s\n",
              lab[k], sum(g < -tol, na.rm = TRUE), sum(is.na(g)),
              format(min(g, na.rm = TRUE), digits = 5)))
}
cat("  lane's claim: 6, 5, 2, 0, 1, 0\n")
