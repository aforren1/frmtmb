# What does this change cost at `library(frmtmb)`?
#
# CHECK THE INSTRUMENT BEFORE BELIEVING THE MEASUREMENT.  The first
# version of this script used `system.time()[["elapsed"]]`, which is
# `proc.time()`, and reported +0.040 s.  Review falsified it: that
# clock advances on this box in steps of min 0.010 s and median
# 0.020 s, so the 0.220 against 0.180 minima were one or two ticks and
# the CONTROL moved a whole tick, half the claimed effect.  This is the
# third timing claim to evaporate on that instrument in this project.
#
# So this version times every replicate on BOTH clocks, prints the
# measured tick of each, and reports the high-resolution number.  The
# arms are interleaved so machine load moves them together, there are
# two controls that must report about zero, and the difference on the
# paired minima carries a permutation test rather than an adjective.
#
#   Rscript dev/generics-loadcost.R <reps> <perms>
av <- commandArgs(trailingOnly = TRUE)
reps <- if (length(av) >= 1L) as.integer(av[1]) else 30L
nperm <- if (length(av) >= 2L) as.integer(av[2]) else 20000L
RS <- file.path(R.home("bin"), "Rscript")
FIX <- "C:/Users/adf44/source/r/generics-lib"
BASE <- "C:/Users/adf44/source/r/rellib-r3"
# the review's library, which holds the SWAP build of this lane.
# READ ONLY; nothing is ever installed there.
SWAP <- "C:/Users/adf44/source/r/genrev-lib"
COMMON <- c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6")

tickof <- function(f, want = 200L, budget = 5) {
  d <- numeric(0)
  t0 <- Sys.time()
  repeat {
    v <- replicate(20000, f())
    e <- diff(v)
    d <- c(d, e[e > 0])
    if (length(d) >= want) break
    if (as.numeric(difftime(Sys.time(), t0,
                            units = "secs")) > budget) break
  }
  if (!length(d)) return(c(min = NA_real_, median = NA_real_))
  c(min = min(d), median = stats::median(d))
}

one <- function(lib, expr) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(
    sprintf(".libPaths(c(%s))",
            paste(deparse(c(lib, BASE, COMMON)), collapse = "")),
    "p0 <- proc.time()[['elapsed']]; h0 <- Sys.time()",
    "suppressMessages({",
    expr,
    "})",
    "p1 <- proc.time()[['elapsed']]; h1 <- Sys.time()",
    "cat('ELAPSED', p1 - p0,",
    "    as.numeric(difftime(h1, h0, units = 'secs')), '\\n')"), f)
  o <- suppressWarnings(system2(RS, c("--vanilla", shQuote(f)),
                                stdout = TRUE, stderr = TRUE))
  v <- grep("^ELAPSED", o, value = TRUE)
  if (!length(v)) return(c(NA_real_, NA_real_))
  as.numeric(strsplit(trimws(sub("^ELAPSED ", "", v[1])), " +")[[1]])
}

arms <- list(
  `frmtmb FIX (active)` = list(FIX,  "library(frmtmb)"),
  `frmtmb SWAP (round 1)` = list(SWAP, "library(frmtmb)"),
  `frmtmb BASE` = list(BASE, "library(frmtmb)"),
  `BASE + nlme + gen` = list(
    BASE,
    "library(frmtmb); loadNamespace('nlme'); loadNamespace('generics')"),
  `control Matrix FIX`  = list(FIX,  "library(Matrix)"),
  `control Matrix BASE` = list(BASE, "library(Matrix)"),
  `control nothing FIX`  = list(FIX,  "invisible(NULL)"),
  `control nothing BASE` = list(BASE, "invisible(NULL)"),
  `lme4 (2.5c price)` = list(FIX, "library(lme4)"),
  `nlme (2.5a price)` = list(FIX, "library(nlme)"),
  `generics (2.5a price)` = list(FIX, "library(generics)")
)

st <- hr <- matrix(NA_real_, reps, length(arms),
                   dimnames = list(NULL, names(arms)))
for (i in seq_len(reps)) {
  for (j in seq_along(arms)) {
    v <- one(arms[[j]][[1]], arms[[j]][[2]])
    st[i, j] <- v[1]
    hr[i, j] <- v[2]
  }
}

pt <- tickof(function() proc.time()[["elapsed"]])
ht <- tickof(function() as.numeric(Sys.time()))

perm_p <- function(a, b, nperm) {
  obs <- min(a) - min(b)
  n <- length(a)
  cnt <- 0L
  for (k in seq_len(nperm)) {
    s <- stats::runif(n) < 0.5
    aa <- ifelse(s, b, a)
    bb <- ifelse(s, a, b)
    if (abs(min(aa) - min(bb)) >= abs(obs) - 1e-12) cnt <- cnt + 1L
  }
  c(obs = obs, p = cnt / nperm)
}

cat("```\n")
cat("== load cost, dev/generics-loadcost.R, R 4.6.1 ==\n")
cat(sprintf("%d replicates, interleaved, one fresh process each\n", reps))
cat(sprintf("proc.time() tick here: min %.4f s, median %.4f s\n",
            pt[["min"]], pt[["median"]]))
cat(sprintf("Sys.time()  tick here: min %.6f s, median %.6f s\n",
            ht[["min"]], ht[["median"]]))
cat("The first is the instrument that produced the withdrawn 0.040 s.\n\n")
cat(sprintf("%-22s %8s %8s | %9s %9s\n", "arm", "min(pt)", "med(pt)",
            "min(hr)", "med(hr)"))
for (j in seq_along(arms)) {
  cat(sprintf("%-22s %8.3f %8.3f | %9.4f %9.4f\n", names(arms)[j],
              min(st[, j]), stats::median(st[, j]),
              min(hr[, j]), stats::median(hr[, j])))
}
set.seed(20260915)
pp <- perm_p(hr[, "frmtmb FIX (active)"], hr[, "frmtmb BASE"],
             nperm)
sw <- perm_p(hr[, "frmtmb SWAP (round 1)"], hr[, "frmtmb BASE"],
             nperm)
ab <- perm_p(hr[, "frmtmb FIX (active)"],
             hr[, "frmtmb SWAP (round 1)"], nperm)
cm <- perm_p(hr[, "control Matrix FIX"], hr[, "control Matrix BASE"],
             nperm)
cn <- perm_p(hr[, "control nothing FIX"], hr[, "control nothing BASE"],
             nperm)
cat(sprintf("\nACTIVE - BASE on minima, hi-res:  %+.4f s  p = %.4f\n",
            pp[["obs"]], pp[["p"]]))
cat(sprintf("SWAP   - BASE, same:              %+.4f s  p = %.4f\n",
            sw[["obs"]], sw[["p"]]))
cat(sprintf("ACTIVE - SWAP, same:              %+.4f s  p = %.4f\n",
            ab[["obs"]], ab[["p"]]))
cat(sprintf("CONTROL Matrix, same:             %+.4f s  p = %.4f\n",
            cm[["obs"]], cm[["p"]]))
cat(sprintf("CONTROL nothing, same:            %+.4f s  p = %.4f\n",
            cn[["obs"]], cn[["p"]]))
cat(sprintf("ACTIVE - BASE on proc.time():      %+.4f s (%.1f ticks)\n",
            min(st[, "frmtmb FIX (active)"]) -
              min(st[, "frmtmb BASE"]),
            (min(st[, "frmtmb FIX (active)"]) -
               min(st[, "frmtmb BASE"])) / pt[["min"]]))
cat("\nATTRIBUTION, hi-res minima:\n")
cat(sprintf("  BASE + nlme + generics, minus BASE: %+.4f s\n",
            min(hr[, "BASE + nlme + gen"]) - min(hr[, "frmtmb BASE"])))
cat(sprintf("  ACTIVE minus that arm:              %+.4f s\n",
            min(hr[, "frmtmb FIX (active)"]) -
              min(hr[, "BASE + nlme + gen"])))
cat(sprintf("  SWAP   minus that arm:              %+.4f s\n",
            min(hr[, "frmtmb SWAP (round 1)"]) -
              min(hr[, "BASE + nlme + gen"])))
cat(sprintf("%d permutations, seed 20260915\n", nperm))
cat("```\n")
saveRDS(list(st = st, hr = hr, pt = pt, ht = ht),
        "dev/generics-out2/loadcost.rds")
