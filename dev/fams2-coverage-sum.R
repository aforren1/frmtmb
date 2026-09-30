# Summarize dev/fams2-cov/*.rds (dev/fams2-coverage.R): per family and
# parameter, the 95% Wald interval's coverage with its Wilson interval,
# the mean estimate against the truth, and the error and warning counts.
# Counted from the files on disk.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-cov"
files <- list.files(dir, pattern = "[.]rds$", full.names = TRUE)
res <- lapply(files, readRDS)
fam <- vapply(res, `[[`, "", "family")
wilson <- function(k, n, z = stats::qnorm(0.975)) {
  p <- k / n
  c <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(c - h, c + h)
}
for (f in unique(fam)) {
  rr <- res[fam == f]
  seeds <- vapply(rr, `[[`, 0L, "seed")
  ok <- vapply(rr, function(r) is.null(r$error), NA)
  nwarn <- sum(vapply(rr[ok], function(r) length(r$warn) > 0L, NA))
  cat(sprintf("\n%s: %d result files, seeds %d..%d (%d distinct), %d errors, %d with a warning\n",
              f, length(rr), min(seeds), max(seeds), length(unique(seeds)),
              sum(!ok), nwarn))
  if (any(!ok)) print(table(vapply(rr[!ok], `[[`, "", "error")))
  if (nwarn) print(table(unlist(lapply(rr[ok], `[[`, "warn"))))
  cov <- do.call(rbind, lapply(rr[ok], `[[`, "cover"))
  est <- do.call(rbind, lapply(rr[ok], `[[`, "est"))
  truth <- rr[[which(ok)[1]]]$truth
  n <- nrow(cov)
  for (p in colnames(cov)) {
    k <- sum(cov[, p])
    w <- wilson(k, n)
    # the mean estimate's distance from the truth in Monte Carlo
    # standard errors of that mean
    bz <- (mean(est[, p]) - truth[[p]]) / (stats::sd(est[, p]) / sqrt(n))
    cat(sprintf(paste0("  %-18s cover %3d/%d = %.3f  Wilson [%.3f, %.3f]",
                       "  mean est %8.4f  truth %8.4f  sd %.4f  bias/mcse %5.2f\n"),
                p, k, n, k / n, w[1], w[2], mean(est[, p]), truth[[p]],
                stats::sd(est[, p]), bz))
  }
}
