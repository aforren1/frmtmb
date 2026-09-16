# Check the INSTRUMENT behind "+0.0000 s [-0.0063, +0.0059]". A null
# result is only worth something if the instrument can see an effect
# of the size that would matter, so this run carries a POSITIVE
# control: the BASE arm with a known 10 ms delay added, which must
# read +0.010 s. Arms interleaved, one fresh process per measurement,
# order shuffled within each round.
#   Rscript dev/sgrev-loadcost.R <rounds>
av <- commandArgs(trailingOnly = TRUE)
R <- if (length(av)) as.integer(av[[1]]) else 30L
BASE <- "C:/Users/adf44/source/r/rellib-r3"
FIX <- "C:/Users/adf44/source/r/sgrev-lib"
PIN <- "C:/Users/adf44/source/r/pinlib"
USR <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
rbin <- file.path(R.home("bin"), "Rscript")

# the tick of the clock the child uses, measured here
tk <- diff(sort(unique(replicate(20000, as.numeric(Sys.time())))))
cat("Sys.time() tick: min ", format(min(tk), digits = 3), " s  median ",
    format(median(tk), digits = 3), " s\n", sep = "")

mk <- function(lib, pkg, delay) {
  f <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s","%s","%s"))', lib, PIN, USR),
    't0 <- as.numeric(Sys.time())',
    if (delay > 0) sprintf('Sys.sleep(%f)', delay) else '',
    sprintf('suppressMessages(library(%s))', pkg),
    't1 <- as.numeric(Sys.time())',
    'cat(t1 - t0, "\\n")'), f)
  f
}
arms <- list(
  sample_BASE = mk(BASE, "frmtmb.sample", 0),
  sample_FIX = mk(FIX, "frmtmb.sample", 0),
  core_BASE = mk(BASE, "frmtmb", 0),
  core_FIX = mk(FIX, "frmtmb", 0),
  POSCTRL_BASE_plus_10ms = mk(BASE, "frmtmb.sample", 0.010))
nm <- names(arms)
res <- matrix(NA_real_, R, length(nm), dimnames = list(NULL, nm))
set.seed(20260915)
for (r in seq_len(R)) {
  for (a in sample(nm)) {
    o <- system2(rbin, c("--vanilla", shQuote(arms[[a]])), stdout = TRUE,
                 stderr = FALSE)
    res[r, a] <- as.numeric(o[length(o)])
  }
}
cat("\n", R, " rounds, ", R * length(nm),
    " fresh processes, arms shuffled within each round\n", sep = "")
cat(sprintf("%-24s %8s %8s\n", "arm", "min", "median"))
for (a in nm) cat(sprintf("%-24s %8.4f %8.4f\n", a, min(res[, a]),
                          median(res[, a])))
boot <- function(d, B = 20000) {
  s <- replicate(B, median(sample(d, length(d), TRUE)))
  quantile(s, c(0.025, 0.975))
}
cat("\npaired per-round difference, median [95% bootstrap interval]\n")
prs <- list(c("sample_FIX", "sample_BASE"),
            c("core_FIX", "core_BASE"),
            c("POSCTRL_BASE_plus_10ms", "sample_BASE"))
for (p in prs) {
  d <- res[, p[1]] - res[, p[2]]
  ci <- boot(d)
  cat(sprintf("%-34s %+8.4f s [%+.4f, %+.4f]\n",
              paste(p[1], "-", p[2]), median(d), ci[1], ci[2]))
}
cat("\nthe positive control must read about +0.0100 s. If it does,",
    "an effect\nof that size would not have been missed.\n")
saveRDS(res, "dev/sgrev-out/loadcost.rds")
cat("DONE\n")
