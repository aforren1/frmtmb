# genrev: re-measure the load cost with the instrument checked first.
#
# The lane used system.time()[["elapsed"]], whose tick on Windows is
# coarse; a +0.040 s claim on that clock can be two ticks.  This script
# measures the tick, then re-times with Sys.time(), which on R 4.x
# Windows reads a high-resolution counter, and reports BOTH clocks so
# the two can be compared.
av <- commandArgs(trailingOnly = TRUE)
reps <- if (length(av)) as.integer(av[1]) else 30L
RS <- file.path(R.home("bin"), "Rscript")
FIX <- "C:/Users/adf44/source/r/genrev-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r3"
COMMON <- c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6")

# --- the instrument ------------------------------------------------
tick <- function(f, n = 400000L) {
  v <- numeric(0)
  t0 <- f()
  for (i in seq_len(n)) {
    t1 <- f()
    if (!isTRUE(all.equal(t1, t0))) { v <- c(v, t1 - t0); t0 <- t1 }
    if (length(v) > 200) break
  }
  v
}
pt <- tick(function() proc.time()[["elapsed"]])
st <- tick(function() as.numeric(Sys.time()))
cat(sprintf("proc.time() elapsed tick: min %.6f s median %.6f s (n=%d)\n",
            min(pt), median(pt), length(pt)))
cat(sprintf("Sys.time()   tick:        min %.9f s median %.9f s (n=%d)\n",
            min(st), median(st), length(st)))

one <- function(lib, expr) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(
    sprintf(".libPaths(c(%s))",
            paste(deparse(c(lib, COMMON)), collapse = "")),
    "a <- Sys.time()",
    "t <- system.time(suppressMessages({",
    expr,
    "}))",
    "b <- Sys.time()",
    "cat('ELAPSED', t[['elapsed']], as.numeric(b - a, units = 'secs'), '\\n')"),
    f)
  o <- suppressWarnings(system2(RS, c("--vanilla", shQuote(f)),
                                stdout = TRUE, stderr = TRUE))
  v <- grep("^ELAPSED", o, value = TRUE)
  if (!length(v)) return(c(NA_real_, NA_real_))
  as.numeric(strsplit(trimws(sub("^ELAPSED ", "", v[1])), " +")[[1]])
}

arms <- list(
  `frmtmb FIX`  = list(FIX,  "library(frmtmb)"),
  `frmtmb BASE` = list(BASE, "library(frmtmb)"),
  # the attribution test: BASE plus exactly the namespaces the FIX
  # build adds to Imports
  `frmtmb BASE + nlme + generics` =
    list(BASE, "loadNamespace('nlme'); loadNamespace('generics'); library(frmtmb)"),
  `control Matrix FIX`  = list(FIX,  "library(Matrix)"),
  `control Matrix BASE` = list(BASE, "library(Matrix)"),
  `control nothing FIX`  = list(FIX,  "invisible(1)"),
  `control nothing BASE` = list(BASE, "invisible(1)"),
  `nlme alone` = list(FIX, "library(nlme)")
)

res <- array(NA_real_, c(reps, length(arms), 2),
             dimnames = list(NULL, names(arms), c("systime", "systime_hr")))
for (i in seq_len(reps)) {
  for (j in seq_along(arms)) {
    res[i, j, ] <- one(arms[[j]][[1]], arms[[j]][[2]])
  }
}
cat(sprintf("\n%d replicates, interleaved, one fresh process each\n", reps))
cat(sprintf("%-32s %8s %8s | %10s %10s %10s\n", "arm", "min(st)",
            "med(st)", "min(hr)", "med(hr)", "q10(hr)"))
for (j in seq_along(arms)) {
  a <- res[, j, 1]; b <- res[, j, 2]
  cat(sprintf("%-32s %8.3f %8.3f | %10.4f %10.4f %10.4f\n", names(arms)[j],
              min(a, na.rm = TRUE), median(a, na.rm = TRUE),
              min(b, na.rm = TRUE), median(b, na.rm = TRUE),
              quantile(b, 0.1, na.rm = TRUE)))
}
g <- function(nm, k) res[, nm, k]
cat(sprintf("\nFIX - BASE      minima: system.time %+.3f  hi-res %+.4f\n",
            min(g("frmtmb FIX", 1)) - min(g("frmtmb BASE", 1)),
            min(g("frmtmb FIX", 2)) - min(g("frmtmb BASE", 2))))
cat(sprintf("CONTROL Matrix  minima: system.time %+.3f  hi-res %+.4f\n",
            min(g("control Matrix FIX", 1)) - min(g("control Matrix BASE", 1)),
            min(g("control Matrix FIX", 2)) - min(g("control Matrix BASE", 2))))
cat(sprintf("CONTROL nothing minima: system.time %+.3f  hi-res %+.4f\n",
            min(g("control nothing FIX", 1)) - min(g("control nothing BASE", 1)),
            min(g("control nothing FIX", 2)) - min(g("control nothing BASE", 2))))
cat(sprintf("BASE+nlme+generics - BASE minima: hi-res %+.4f  (the attribution)\n",
            min(g("frmtmb BASE + nlme + generics", 2)) -
              min(g("frmtmb BASE", 2))))
cat(sprintf("FIX - (BASE+nlme+generics) minima: hi-res %+.4f  (residual)\n",
            min(g("frmtmb FIX", 2)) -
              min(g("frmtmb BASE + nlme + generics", 2))))
# a permutation check on the FIX-vs-BASE minima: how often does a random
# relabelling of the paired replicates produce a gap this large?
d <- cbind(g("frmtmb FIX", 2), g("frmtmb BASE", 2))
obs <- min(d[, 1]) - min(d[, 2])
set.seed(20260915)
perm <- replicate(20000, {
  s <- sample(c(TRUE, FALSE), nrow(d), TRUE)
  x <- ifelse(s, d[, 1], d[, 2]); y <- ifelse(s, d[, 2], d[, 1])
  min(x) - min(y)
})
cat(sprintf("permutation p for the FIX-BASE gap on minima: %.4f (obs %+.4f)\n",
            mean(abs(perm) >= abs(obs)), obs))
saveRDS(res, file.path(dirname(tempdir()), "genrev-loadcost.rds"))
cat("GENREVDONE\n")
