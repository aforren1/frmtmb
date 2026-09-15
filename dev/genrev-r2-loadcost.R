# genrev round 2: settle whether +0.0167 s (round-1 review, SWAP build)
# was the WHOLE change or only the Imports.
#
# Round 1 took the minimum over 30 replicates per arm.  The lane's
# three runs of the same statistic spread 7 ms, which is the size of
# the disagreement, so a minimum over 30 cannot settle it.  This run:
#  * 80 interleaved ROUNDS, each round one fresh process per arm, the
#    arm order shuffled within the round;
#  * a hi-res clock (Sys.time) inside the child around library() only;
#  * the statistic is the PAIRED within-round difference, summarized by
#    its median with a bootstrap interval, beside the difference of
#    minima and of 10th percentiles, so the three can be compared;
#  * two controls that must read zero: an empty child on each library.
av <- commandArgs(trailingOnly = TRUE)
rounds <- if (length(av)) as.integer(av[1]) else 80L
RS <- file.path(R.home("bin"), "Rscript")
COMMON <- c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6")
LIBS <- c(BASE = "C:/Users/adf44/source/r/rellib-r3",
          SWAP = "C:/Users/adf44/source/r/genrev-lib",
          ACTIVE = "C:/Users/adf44/source/r/genrev2-lib")
arms <- list(
  BASE = list(LIBS[["BASE"]], "library(frmtmb)"),
  IMPORTS = list(LIBS[["BASE"]],
    "loadNamespace('nlme'); loadNamespace('generics'); library(frmtmb)"),
  SWAP = list(LIBS[["SWAP"]], "library(frmtmb)"),
  ACTIVE = list(LIBS[["ACTIVE"]], "library(frmtmb)"),
  NULL_BASE = list(LIBS[["BASE"]], "invisible(1)"),
  NULL_ACTIVE = list(LIBS[["ACTIVE"]], "invisible(1)"))
one <- function(lib, expr) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(sprintf(".libPaths(c(%s))",
                       paste(deparse(c(lib, COMMON)), collapse = "")),
               "a <- Sys.time()", "suppressMessages({", expr, "})",
               "b <- Sys.time()",
               "cat('HR', format(as.numeric(b - a, units = 'secs'), digits = 9), '\\n')"),
             f)
  o <- suppressWarnings(system2(RS, c("--vanilla", shQuote(f)),
                                stdout = TRUE, stderr = TRUE))
  v <- grep("^HR", o, value = TRUE)
  if (!length(v)) return(NA_real_)
  as.numeric(sub("^HR ", "", trimws(v[1])))
}
set.seed(20260915)
res <- matrix(NA_real_, rounds, length(arms), dimnames = list(NULL, names(arms)))
for (i in seq_len(rounds)) {
  for (j in sample(seq_along(arms))) res[i, j] <- one(arms[[j]][[1]], arms[[j]][[2]])
}
saveRDS(res, "C:/Users/adf44/source/r/frmtmb-wt-generics/dev/genrev-out/r2-loadcost.rds")
cat(sprintf("%d rounds, arms shuffled within each round, hi-res clock\n\n", rounds))
cat(sprintf("%-12s %8s %8s %8s\n", "arm", "min", "q10", "median"))
for (nm in names(arms)) {
  v <- res[, nm]
  cat(sprintf("%-12s %8.4f %8.4f %8.4f\n", nm, min(v, na.rm = TRUE),
              quantile(v, 0.1, na.rm = TRUE), median(v, na.rm = TRUE)))
}
cmp <- function(a, b) {
  d <- res[, a] - res[, b]
  d <- d[is.finite(d)]
  bs <- replicate(4000, median(sample(d, replace = TRUE)))
  cat(sprintf("%-18s paired median %+.4f [%+.4f, %+.4f]  min-diff %+.4f  q10-diff %+.4f\n",
              paste(a, "-", b), median(d), quantile(bs, 0.025), quantile(bs, 0.975),
              min(res[, a], na.rm = TRUE) - min(res[, b], na.rm = TRUE),
              quantile(res[, a], 0.1, na.rm = TRUE) - quantile(res[, b], 0.1, na.rm = TRUE)))
}
cat("\n")
cmp("IMPORTS", "BASE")
cmp("SWAP", "BASE")
cmp("ACTIVE", "BASE")
cmp("SWAP", "IMPORTS")
cmp("ACTIVE", "IMPORTS")
cmp("ACTIVE", "SWAP")
cmp("NULL_ACTIVE", "NULL_BASE")
# how much do minima over 30 move?  resample 30 of the rounds 2000 times
mins30 <- replicate(2000, {
  k <- sample(rounds, 30)
  min(res[k, "SWAP"]) - min(res[k, "BASE"])
})
cat(sprintf("\nSWAP-BASE on minima of a random 30 rounds: 5%% %+.4f  50%% %+.4f  95%% %+.4f\n",
            quantile(mins30, 0.05), quantile(mins30, 0.5), quantile(mins30, 0.95)))
cat("GENREVDONE\n")
