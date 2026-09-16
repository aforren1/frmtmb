# What the 28 bindings cost frmtmb.sample at load.
#
# Each replicate is a FRESH process timing loadNamespace() on a
# high-resolution clock (Sys.time(), whose tick is printed; proc.time()
# ticks at 10 ms here and produced three withdrawn figures in this
# project). Arms are interleaved and shuffled within each round. The
# CONTROL arm loads frmtmb alone from each library, which the change to
# core touches only by one argument check, so it should show no
# difference; it is the instrument's own noise floor.
#
#   Rscript dev/samplegen-loadcost.R <FIXLIB> <rounds> <seed>
av <- commandArgs(trailingOnly = TRUE)
FIX <- av[1]
rounds <- as.integer(av[2])
set.seed(as.integer(av[3]))
BASE <- "C:/Users/adf44/source/r/rellib-r3"
RS <- file.path(R.home("bin"), "Rscript")
child <- function(lib, pkg) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f))
  writeLines(c(
    sprintf(".libPaths(c('%s', '%s', '%s', '%s'))", lib, BASE,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"),
    "t0 <- Sys.time()",
    sprintf("suppressMessages(loadNamespace('%s'))", pkg),
    "cat(format(as.numeric(Sys.time() - t0, units = 'secs'),",
    "           digits = 9), '\\n')"), f)
  o <- system2(RS, c("--vanilla", shQuote(f)), stdout = TRUE)
  # the timing is the last line; anything a package prints goes before it
  as.numeric(utils::tail(o, 1L))
}
ticks <- diff(unique(vapply(1:2000, function(i) as.numeric(Sys.time()), 0)))
cat(sprintf("Sys.time() tick here: min %.7f s, median %.7f s\n",
            min(ticks), median(ticks)))
arms <- list(
  sample_BASE = c(BASE, "frmtmb.sample"),
  sample_FIX = c(FIX, "frmtmb.sample"),
  control_core_BASE = c(BASE, "frmtmb"),
  control_core_FIX = c(FIX, "frmtmb"))
res <- matrix(NA_real_, rounds, length(arms), dimnames = list(NULL, names(arms)))
for (r in seq_len(rounds)) {
  for (a in sample(names(arms))) {
    res[r, a] <- child(arms[[a]][1], arms[[a]][2])
  }
}
cat(sprintf("\n%d rounds, %d fresh processes, arms shuffled per round\n",
            rounds, rounds * length(arms)))
cat(sprintf("%-20s %9s %9s\n", "arm", "min", "median"))
for (a in names(arms)) {
  cat(sprintf("%-20s %9.4f %9.4f\n", a, min(res[, a]), median(res[, a])))
}
boot <- function(d, B = 20000) {
  m <- replicate(B, median(sample(d, replace = TRUE)))
  c(median(d), quantile(m, c(0.025, 0.975)))
}
ds <- boot(res[, "sample_FIX"] - res[, "sample_BASE"])
dc <- boot(res[, "control_core_FIX"] - res[, "control_core_BASE"])
cat(sprintf("\npaired per-round difference, median [95%% bootstrap interval]\n"))
cat(sprintf("frmtmb.sample FIX - BASE   %+.4f s [%+.4f, %+.4f]\n",
            ds[1], ds[2], ds[3]))
cat(sprintf("CONTROL frmtmb FIX - BASE  %+.4f s [%+.4f, %+.4f]\n",
            dc[1], dc[2], dc[3]))
cat(sprintf("seed %s, B = 20000\n", av[3]))
