# Short check of lane squant fix round 1, item 3: elf_quantile() round
# trip pr -> q -> pr at extreme pr and small shapes, against pbeta() on
# the log scale (each tail from its own side) and against the package's
# own elf_cdf_log(); the review's three cases against its uniroot()
# roots.
# Run (worktree root):
#   Rscript dev/reviews/squant-short-quant.R > dev/reviews/squant-short-quant.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
eq <- ns$elf_quantile
# log F and log S at t from pbeta, each tail from its own side, x and
# 1 - x formed without cancellation
pb_pair <- function(t, p, lam) {
  a <- lam * (1 - p)
  b <- lam * p
  x <- stats::plogis(t)
  xc <- stats::plogis(-t)
  lF <- ifelse(t <= 0, stats::pbeta(x, a, b, log.p = TRUE),
               stats::pbeta(xc, b, a, lower.tail = FALSE, log.p = TRUE))
  lS <- ifelse(t <= 0, stats::pbeta(x, a, b, lower.tail = FALSE, log.p = TRUE),
               stats::pbeta(xc, b, a, log.p = TRUE))
  list(lF = lF, lS = lS)
}
cat("== review cases (mu 0, sigma 1) against the review's uniroot roots ==\n")
cases <- list(c(0.01, 0.25, 0.9, 229.34003), c(0.01, 0.05, 0.5, 68.313613),
              c(0.01, 0.001, 0.5, 68.309686),
              c(0.99, 0.25, 0.1, -229.34003), c(0.99, 0.05, 0.5, -68.313613),
              c(0.99, 0.001, 0.5, -68.309686))
for (cc in cases) {
  w <- character()
  q <- withCallingHandlers(eq(cc[3], 0, 1, cc[1], cc[2]),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  cat(sprintf("p %.2f lam %g pr %.2f: q %.8f (review %.8f, diff %.2g) warnings %d\n",
              cc[1], cc[2], cc[3], q, cc[4], q - cc[4], length(w)))
}
cat("\n== round trip; err = |log pr' - log pr| for pr <= 0.5, else on log(1 - pr) ==\n")
prs <- c(1e-300, 1e-100, 1e-20, 1e-8, 1e-3, 0.1, 0.5, 0.9, 1 - 1e-3,
         1 - 1e-8, 1 - 1e-15)
worst <- 0
worst_pkg <- 0
nbad <- 0
for (p in c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)) {
  for (lam in c(0.001, 0.01, 0.05, 0.25, 1, 4)) {
    for (mu_s in list(c(0, 1), c(3, 2))) {
      w <- character()
      q <- withCallingHandlers(eq(prs, mu_s[1], mu_s[2], p, lam),
                               warning = function(x) {
                                 w <<- c(w, conditionMessage(x))
                                 invokeRestart("muffleWarning")
                               })
      t <- (q - mu_s[1]) / (lam * mu_s[2])
      pb <- suppressWarnings(pb_pair(t, p, lam))
      ck <- ns$elf_cdf_log(q, mu_s[1], mu_s[2], p, 1 - p, lam)
      up <- prs > 0.5
      goal <- ifelse(up, log1p(-prs), log(prs))
      got_pb <- ifelse(up, pb$lS, pb$lF)
      got_pk <- ifelse(up, ck$lS, ck$lF)
      e_pb <- abs(got_pb - goal)
      e_pk <- abs(got_pk - goal)
      worst <- max(worst, e_pb[is.finite(e_pb)])
      worst_pkg <- max(worst_pkg, e_pk, na.rm = FALSE)
      bad <- !is.finite(q) | !is.finite(e_pk) | e_pk > 1e-10
      nbad <- nbad + sum(bad)
      if (mu_s[1] == 0 || any(bad)) {
        cat(sprintf(paste0("p %.3f lam %.3f mu %g sigma %g: max err pkg cdf %.2g, ",
                           "pbeta %.2g | t range %.4g..%.4g | nonfinite q %d ",
                           "| warnings %d%s\n"),
                    p, lam, mu_s[1], mu_s[2], max(e_pk), max(e_pb),
                    min(t), max(t), sum(!is.finite(q)), length(w),
                    if (any(bad)) paste0(" | BAD at pr ",
                                         paste(format(prs[bad], digits = 3),
                                               collapse = ",")) else ""))
      }
    }
  }
}
cat(sprintf("SUMMARY: worst round-trip err on the package cdf %.3g, on pbeta %.3g; points above 1e-10 or nonfinite: %d of %d\n",
            worst_pkg, worst, nbad, length(prs) * 7 * 6 * 2))
cat("DONE\n")
